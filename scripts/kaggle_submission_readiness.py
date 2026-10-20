"""Offline, read-only CASMI26 Kaggle submission readiness audit.

Checks only local data availability, file schemas and the presence of source
scripts. Never contacts Kaggle, reads tokens, downloads data or submits.
Can run inside a prepared Kaggle Notebook or an authorized local dataset
directory. Missing prerequisites cause a nonzero exit code.

Usage:
  python scripts/kaggle_submission_readiness.py --data-dir /kaggle/input/COMPETITION
  python scripts/kaggle_submission_readiness.py --data-dir /approved/casmi --json-out runs/readiness.json
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

NEEDED = ("train.parquet", "test.parquet", "sample_submission.csv")
REQUIRED_TEST = {"molecule_id"}
REQUIRED_TRAIN = {"ingest_lib", "precursor_mz", "adduct"}


def evaluate(data_dir: Path, *, repo_root: Path | None = None) -> dict:
    base = Path(data_dir)
    repo = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]
    files = {}
    for name in NEEDED:
        path = base / name
        files[name] = {"present": path.is_file(),
                       "nonempty": path.is_file() and path.stat().st_size > 0}
    scripts = {
        path: (repo / path).is_file()
        for path in (
            "kaggle/submission_driver.py", "scripts/audit_submission.py",
            "kaggle/ours_channel.py",
        )
    }
    problems = []
    for key, value in files.items():
        if not value["present"] or not value["nonempty"]:
            problems.append("Missing or empty official competition file: " + key)
    for key, good in scripts.items():
        if not good:
            problems.append("Missing repository notebook component: " + key)
    if files["sample_submission.csv"]["nonempty"]:
        with (base / "sample_submission.csv").open(
            "r", newline="", encoding="utf-8-sig"
        ) as source:
            header = next(csv.reader(source), [])
        if header != ["molecule_id", "smiles"]:
            problems.append("Unexpected sample submission columns")
    # Do not request pyarrow merely to verify whether files exist.
    # When available, inspect Parquet footer/schema without scanning rows.
    parquet_schema = {}
    if all(files[name]["nonempty"] for name in ("train.parquet", "test.parquet")):
        try:
            import pyarrow.parquet as pq
        except ImportError:
            problems.append("pyarrow not installed; Parquet schemas unverified")
        else:
            for name, required in (
                ("train.parquet", REQUIRED_TRAIN),
                ("test.parquet", REQUIRED_TEST),
            ):
                try:
                    file = pq.ParquetFile(base / name)
                    cols = set(file.schema_arrow.names)
                    parquet_schema[name] = {
                        "rows": int(file.metadata.num_rows),
                        "required_columns_present": required <= cols,
                    }
                    if file.metadata.num_rows == 0 or not required <= cols:
                        problems.append("Empty or incompatible Parquet schema: " + name)
                except (OSError, ValueError, TypeError) as exc:
                    problems.append("Unreadable Parquet footer: " + name)
    return {
        "status": "LOCAL_INPUTS_VERIFIED" if not problems else "BLOCKED",
        "files": files, "repository_scripts": scripts,
        "parquet_schema": parquet_schema, "blockers": problems,
        "submitted": False, "scored": False,
        "disclaimer": "Local input/schema check only; no Kaggle permission, runtime, submission or MRR verification",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = evaluate(args.data_dir)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        with args.json_out.open("x", encoding="utf-8") as out:
            out.write(rendered)
    print(rendered, end="")
    if report["status"] != "LOCAL_INPUTS_VERIFIED":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
