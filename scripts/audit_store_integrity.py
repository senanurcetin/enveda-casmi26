"""Read-only structural audit of an existing CASMI MS/MS peak store.

Usage:
  python scripts/audit_store_integrity.py \
    --store /path/to/artifacts/store \
    --output runs/store_layout_audit.json

Requires pandas, NumPy and a Parquet reader (pyarrow). Never copies peak
arrays into memory, uploads source data or changes any store files. This
audit checks offsets/file lengths and metadata peak counts only; it is
NOT a spectrum-content or upstream data-provenance verification.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.spectra.store_integrity import inspect_store_layout, scan_store_peak_content  # noqa: E402


def inspect_store(root: Path, *, scan_peaks: bool = False,
                  chunk_peaks: int = 262144) -> dict:
    import pandas as pd

    peaks = pd.read_parquet(root / "meta.parquet", columns=["n_peaks"])
    layout = inspect_store_layout(
        root, expected_spectra=len(peaks),
        metadata_n_peaks=peaks["n_peaks"].to_numpy(),
    )
    report = layout.summary()
    if scan_peaks:
        report["peak_content"] = scan_store_peak_content(
            root, chunk_peaks=chunk_peaks,
        )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--scan-peaks", action="store_true",
                        help="opt-in full peak-content scan (bounded-memory, full disk I/O)")
    parser.add_argument("--chunk-peaks", type=int, default=262144,
                        help="peak pairs per scan chunk, default 262144")
    args = parser.parse_args()
    try:
        report = inspect_store(
            args.store, scan_peaks=args.scan_peaks, chunk_peaks=args.chunk_peaks,
        )
        serialized = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x", encoding="utf-8") as file:
                file.write(serialized)
        print(serialized, end="")
        if args.scan_peaks and report["peak_content"]["status"] != "VALID_CLEANED_PEAKS":
            raise SystemExit(2)
    except (OSError, ValueError, TypeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
