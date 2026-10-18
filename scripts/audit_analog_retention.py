"""Measure D2 index-selection bias without loading raw spectra or labels.

The historical AnalogIndex._select method retains top reference spectra per
(structure, adduct), prioritizing timsTOF acquisitions. A novelty knockout can
remove that source while leaving other-library spectra eligible in the original
data but never indexed. This read-only audit reports counts of lost cross-
library reference pairs; it never alters the model/index or claims MRR effects.

Usage:
  python scripts/audit_analog_retention.py \
    --meta artifacts/store/meta.parquet \
    --query-library enveda-np-examples \
    --json-out runs/analog_retention_audit.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_COLUMNS = (
    "train_ik14", "ingest_lib", "adduct", "test_adduct",
    "ppm", "n_peaks", "is_timstof",
)


def audit_retention(meta: pd.DataFrame, query_library: str,
                    per_structure_adduct: int = 3) -> dict:
    """Compare cross-library feasible (key, adduct) pairs before/after indexing.

    Pair counts reflect index-coverage limitations, not molecular retrieval
    performance; no spectra similarities, ranks, or competition labels used.
    """
    missing = sorted(set(EXPECTED_COLUMNS) - set(meta.columns))
    if missing:
        raise ValueError(f"Metadata columns missing: {', '.join(missing)}")
    if not isinstance(query_library, str) or not query_library.strip():
        raise ValueError("query_library must be a nonempty string")
    if type(per_structure_adduct) is not int or per_structure_adduct < 1:
        raise ValueError("per_structure_adduct must be a positive integer")
    if meta.empty:
        raise ValueError("Empty spectrum metadata")
    eligible = (
        meta.test_adduct.fillna(False).to_numpy(dtype=bool)
        & (np.abs(meta.ppm.to_numpy(dtype=float)) <= 20)
        & (meta.n_peaks.to_numpy(dtype=float) >= 5)
    )
    selected = np.flatnonzero(eligible)
    t = pd.DataFrame({
        "row": selected,
        "key": meta.train_ik14.astype("string").to_numpy()[selected],
        "adduct": meta.adduct.astype("string").to_numpy()[selected],
        "lib": meta.ingest_lib.astype("string").to_numpy()[selected],
        "prio": (
            meta.is_timstof.fillna(False).to_numpy(dtype=int)[selected] * 1000
            + meta.n_peaks.to_numpy(dtype=float)[selected].clip(0, 999)
        ),
    })
    t = t.dropna(subset=["key", "adduct", "lib"])
    t = t[(t.key != "") & (t.adduct != "")]
    if t.empty:
        raise ValueError("No eligible reference spectra")
    # Mirrors AnalogIndex._select exactly: highest priority then stable row order
    # with at most N spectra retained per standard-key/adduct.
    ranked = t.sort_values("prio", ascending=False, kind="stable")
    kept = ranked.groupby(["key", "adduct"], observed=True,
                          sort=False).head(per_structure_adduct)
    pairs = lambda frame: set(zip(frame.key, frame.adduct))
    all_other = pairs(t[t.lib != query_library])
    selected_other = pairs(kept[kept.lib != query_library])
    query_lib_pairs = pairs(t[t.lib == query_library])
    cross_possible = all_other & query_lib_pairs
    cross_indexed = selected_other & query_lib_pairs
    lost = cross_possible - cross_indexed
    return {
        "status": "GAP_FOUND" if lost else "NO_GAP_FOUND",
        "method": "metadata-only emulation of AnalogIndex._select",
        "query_library": query_library,
        "per_structure_adduct": per_structure_adduct,
        "eligible_reference_spectra": int(len(t)),
        "indexed_reference_spectra": int(len(kept)),
        "cross_library_candidate_pairs": int(len(cross_possible)),
        "cross_library_pairs_retained": int(len(cross_indexed)),
        "cross_library_pairs_lost_by_index": int(len(lost)),
        "cross_library_pair_loss_fraction": (
            round(len(lost) / len(cross_possible), 6) if cross_possible else None
        ),
        "structures_with_lost_cross_library_pairs": int(len({k for k, _ in lost})),
        "note": (
            "A pair is one (standard 2D key, adduct), not a query molecule. "
            "Instrument/mode, spectral similarity, tautomer equivalence and "
            "actual rankings are not verified by this aggregate."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meta", type=Path, required=True,
                        help="Existing SpectrumStore meta.parquet")
    parser.add_argument("--query-library", default="enveda-np-examples")
    parser.add_argument("--per-structure-adduct", type=int, default=3)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    try:
        meta = pd.read_parquet(args.meta, columns=list(EXPECTED_COLUMNS))
        report = audit_retention(meta, args.query_library, args.per_structure_adduct)
        sha = hashlib.sha256()
        with args.meta.open("rb") as file:
            for chunk in iter(lambda: file.read(1024 * 1024), b""):
                sha.update(chunk)
        report["metadata_sha256"] = sha.hexdigest()
    except (OSError, ValueError, TypeError, KeyError) as exc:
        parser.error(str(exc))
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    print(rendered)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
