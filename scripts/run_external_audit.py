"""Audit pinned external CC0 structures and stored spectra against full CASMI train."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from datetime import UTC, datetime

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def audit(train_path, external_path, manifest_path, protocol_path):
    import pandas as pd
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from casmi.chem.canon import skeleton_key, canonical_ik14
    from casmi.validate.overlap_audit import PublicOverlapAudit
    from casmi.validate.source_manifest import sha256_file

    protocol = json.loads(Path(protocol_path).read_text(encoding="utf-8"))
    source = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    external_path, train_path = Path(external_path), Path(train_path)
    if source["source"] != protocol["external_source"] or sha256_file(external_path) != source["parquet_sha256"]:
        raise ValueError("External source manifest differs from declared inputs")
    external = pd.read_parquet(external_path)
    if set(external.license) != {"CC0"} or external.accession.duplicated().any():
        raise ValueError("External license or accession contract failed")
    targets = external.rename(columns={"accession": "spectrum_id", "metric_key": "molecule_id"}).to_dict("records")
    content = PublicOverlapAudit(targets)
    target_skeletons = set(external.skeleton)
    target_metrics = set(external.metric_key)
    training_standard, overlapping_metrics, seen_smiles = set(), set(), set()
    invalid_structures = matching_skeletons = 0
    train = pq.ParquetFile(train_path)
    precursors = pa.array(sorted(content.precursors), type=pa.float64())
    columns = ["normalized_smiles", "inchikey14", "precursor_mz", "ms2_mzs", "ms2_normalized_intensities", "adduct", "ionization_mode", "ingest_lib"]
    for index, batch in enumerate(train.iter_batches(batch_size=8192, columns=columns)):
        metadata = batch.select(["normalized_smiles", "inchikey14"]).to_pylist()
        for row in metadata:
            key, smiles = row["inchikey14"], row["normalized_smiles"]
            if not isinstance(key, str) or not key:
                raise ValueError("Training identity is missing")
            training_standard.add(key)
            if smiles in seen_smiles:
                continue
            seen_smiles.add(smiles)
            skeleton = skeleton_key(smiles)
            if skeleton is None:
                invalid_structures += 1
            elif skeleton in target_skeletons:
                matching_skeletons += 1
                metric = canonical_ik14(smiles)
                if metric is None:
                    raise ValueError("Matching training skeleton has unverifiable metric identity")
                if metric in target_metrics:
                    overlapping_metrics.add(metric)
        mask = pc.is_in(pc.cast(batch.column("precursor_mz"), pa.float64()), value_set=precursors)
        content.scan_batch(batch.filter(mask).to_pylist(), scanned_count=batch.num_rows)
        if index % 25 == 0:
            print(f"External audit scanned {content.scanned}/{train.metadata.num_rows} training rows", flush=True)
    content_summary = content.finish(expected_train_rows=train.metadata.num_rows)
    # All spectra belonging to a matched external molecule are excluded together.
    content_metrics = {targets[index]["molecule_id"] for index in content.matched}
    standard_metrics = set(external.loc[external.standard_ik14.isin(training_standard), "metric_key"])
    rejected = overlapping_metrics | content_metrics | standard_metrics
    novel = external[~external.metric_key.isin(rejected)].copy()
    novel_metrics = sorted(set(novel.metric_key))
    verified = invalid_structures == 0
    enough = len(novel_metrics) >= 30
    report = {
        "schema_version": 1,
        "status": "EXTERNAL_INDEPENDENCE_FEASIBLE" if verified and enough else "EXTERNAL_INDEPENDENCE_NOT_ESTABLISHED",
        "checked_at": datetime.now(UTC).isoformat(),
        "protocol_sha256": sha256_file(Path(protocol_path)),
        "inputs": {"train.parquet": {"sha256": sha256_file(train_path), "bytes": train_path.stat().st_size},
                   external_path.name: {"sha256": sha256_file(external_path), "bytes": external_path.stat().st_size}},
        "external_source": source,
        "identity_audit": {"training_rows_scanned": content.scanned, "training_standard_identities": len(training_standard),
                           "distinct_training_smiles": len(seen_smiles), "unparseable_training_labels": invalid_structures,
                           "training_smiles_sharing_external_skeleton": matching_skeletons,
                           "external_standard_overlap_molecules": len(standard_metrics),
                           "external_metric_overlap_molecules": len(overlapping_metrics),
                           "external_exact_content_overlap_molecules": len(content_metrics),
                           "excluded_external_molecules_union": len(rejected),
                           "novel_external_molecules": len(novel_metrics), "novel_external_spectra": len(novel),
                           "novel_metric_identities_sha256": canonical_digest(novel_metrics),
                           "novel_accessions_sha256": canonical_digest(sorted(novel.accession))},
        "stored_content_audit": content_summary,
        "model_performance_inspected": False,
        "limitations": ["External source acquisition is independently retrieved public reference data, not a new laboratory replication.",
                        "Identity checks rely on source chemical annotations and the official computational canonicalizer.",
                        "Exact stored-content comparison cannot detect arbitrary preprocessing/normalization variants of a common acquisition.",
                        "This feasibility stage does not establish supervised-family separation, model accuracy or a C3 coverage improvement."]}
    return report, novel


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True, type=Path)
    parser.add_argument("--external", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report, _ = audit(args.train, args.external, args.manifest, args.protocol)
    with args.output.open("x", encoding="utf-8") as destination:
        json.dump(report, destination, indent=2, allow_nan=False)
        destination.write("\n")
    print("External metadata/content feasibility audit completed; no model outcomes inspected", flush=True)
