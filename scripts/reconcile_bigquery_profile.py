"""Reconcile a read-only BigQuery library aggregate with a local CASMI profile.

Inputs must be metadata-only, not row-level spectra or restricted answers.
Expected BigQuery CSV columns: source_library, spectrum_count,
distinct_2d_structures (see sql/bigquery/library_quality.sql).

Example:
  python scripts/reconcile_bigquery_profile.py \
      --profile reports/data_profile.json --aggregate runs/library_quality.csv \
      --json-out runs/bigquery_reconciliation.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

REQUIRED = ("source_library", "spectrum_count", "distinct_2d_structures")


def parse_aggregate(path: Path) -> dict[str, dict[str, int]]:
    with path.open("r", newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None or any(key not in reader.fieldnames for key in REQUIRED):
            raise ValueError(f"BigQuery export requires columns: {', '.join(REQUIRED)}")
        rows = {}
        for line, row in enumerate(reader, start=2):
            source = (row["source_library"] or "").strip()
            if not source:
                raise ValueError(f"row {line}: empty source_library")
            if source in rows:
                raise ValueError(f"row {line}: duplicate source_library {source!r}")
            try:
                spectra = int(row["spectrum_count"])
                structures = int(row["distinct_2d_structures"])
            except (ValueError, TypeError) as exc:
                raise ValueError(f"row {line}: invalid numeric counts") from exc
            if spectra < 0 or structures < 0 or structures > spectra:
                raise ValueError(f"row {line}: impossible spectra/structure counts")
            rows[source] = {"spectra": spectra, "structures": structures}
    if not rows:
        raise ValueError("BigQuery aggregate contains no library rows")
    return rows


def reconcile(profile: dict, observed: dict[str, dict[str, int]]) -> dict:
    if not isinstance(profile.get("train", {}).get("by_library"), dict):
        raise ValueError("Invalid historical data profile")
    expected = profile["train"]["by_library"]
    total = profile["train"].get("spectra")
    if type(total) is not int or total < 1:
        raise ValueError("Invalid historical total spectra")
    issues = []
    for lib in sorted(set(expected) | set(observed)):
        a, b = expected.get(lib), observed.get(lib)
        if a is None:
            issues.append({"source_library": lib, "issue": "unseen_source_library",
                           "expected": None, "observed": b})
        elif b is None:
            issues.append({"source_library": lib, "issue": "missing_source_library",
                           "expected": {"spectra": a["spectra"], "structures": a["structures"]},
                           "observed": None})
        else:
            for metric in ("spectra", "structures"):
                if a[metric] != b[metric]:
                    issues.append({"source_library": lib, "issue": f"{metric}_mismatch",
                                   "expected": a[metric], "observed": b[metric]})
    measured_total = sum(row["spectra"] for row in observed.values())
    if measured_total != total:
        issues.append({"issue": "total_spectra_mismatch", "expected": total,
                       "observed": measured_total})
    # NEVER add per-library distinct structures as a global estimate: a
    # structure can be recorded in multiple source libraries.
    return {
        "status": "MATCH" if not issues else "DRIFT",
        "historical_spectra": total,
        "aggregated_spectra": measured_total,
        "historical_sources": len(expected),
        "aggregate_sources": len(observed),
        "distinct_structure_counts_are_nonadditive": True,
        "issues": issues,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, default=Path("reports/data_profile.json"))
    parser.add_argument("--aggregate", type=Path, required=True,
                        help="CSV exported from the read-only BigQuery aggregate")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    try:
        raw = args.profile.read_bytes()
        profile = json.loads(raw)
        report = reconcile(profile, parse_aggregate(args.aggregate))
        report["historical_profile_sha256"] = hashlib.sha256(raw).hexdigest()
        report["aggregate_file_sha256"] = hashlib.sha256(args.aggregate.read_bytes()).hexdigest()
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    print(rendered)
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered + "\n", encoding="utf-8")
    return 0 if report["status"] == "MATCH" else 2


if __name__ == "__main__":
    raise SystemExit(main())
