"""Read-only, bounded-memory exact public/train overlap audit on authorized Parquet files."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.overlap_audit import PublicOverlapAudit  # noqa: E402

PEAKS = ["precursor_mz", "ms2_mzs", "ms2_normalized_intensities", "adduct", "ionization_mode"]


def file_record(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def run_audit(train_path, public_path, *, batch_size=8192, progress=None):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError("Batch size must be a positive integer")
    train_path, public_path = Path(train_path), Path(public_path)
    if train_path.resolve() == public_path.resolve():
        raise ValueError("Training and public-example paths must be different")
    train = pq.ParquetFile(train_path)
    public = pq.ParquetFile(public_path)
    for parquet, required in ((train, PEAKS + ["ingest_lib"]),
                              (public, PEAKS + ["spectrum_id", "molecule_id"])):
        if not set(required) <= set(parquet.schema_arrow.names):
            raise ValueError("Missing required audit columns")
    # Public inputs should be small. Refuse accidental full-train target materialization.
    if public.metadata.num_rows > 100_000:
        raise ValueError("Public example input exceeds the 100,000-row safety ceiling")
    targets = public.read(columns=PEAKS + ["spectrum_id", "molecule_id"]).to_pylist()
    audit = PublicOverlapAudit(targets)
    precursors = pa.array(sorted(audit.precursors), type=pa.float64())
    for index, batch in enumerate(train.iter_batches(batch_size=batch_size, columns=PEAKS + ["ingest_lib"])):
        mask = pc.is_in(pc.cast(batch.column("precursor_mz"), pa.float64()), value_set=precursors)
        # Exact equality of precursors is necessary for any full-content match.
        # Arrow filters before list materialization; peak arrays of decoys stay in Arrow.
        selected = batch.filter(mask).to_pylist()
        audit.scan_batch(selected, scanned_count=batch.num_rows)
        if progress is not None and index % 25 == 0:
            progress(audit.scanned, train.metadata.num_rows)
    result = audit.finish(expected_train_rows=train.metadata.num_rows)
    result["checked_at"] = datetime.now(UTC).isoformat()
    result["inputs"] = {"train.parquet": file_record(train_path), "test.parquet": file_record(public_path)}
    result["audit_sources"] = {
        name: hashlib.sha256((ROOT / name).read_text(encoding="utf-8").encode()).hexdigest()
        for name in ("src/casmi/validate/spectrum_overlap.py", "src/casmi/validate/overlap_audit.py", "scripts/audit_public_overlap.py")
    }
    result["environment"] = {"python": platform.python_version(), "pyarrow": importlib.metadata.version("pyarrow")}
    result["batch_size"] = batch_size
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train", type=Path, required=True)
    parser.add_argument("--public", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=8192)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Refusing to overwrite an existing audit report")
    result = run_audit(args.train, args.public, batch_size=args.batch_size,
                       progress=lambda count, total: print(f"Scanned {count:,}/{total:,} training rows", flush=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as destination:
        json.dump(result, destination, indent=2, allow_nan=False)
        destination.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
