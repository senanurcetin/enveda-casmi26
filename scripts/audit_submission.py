"""Audit a CASMI26 submission without altering the competition inference pipeline.

By default, this is a dependency-free *structural* preflight, not the official
scorer. Add --chemistry to validate RDKit SMILES and detect tautomer-equivalent
2D candidates. A valid preflight does not establish a high leaderboard score.

Examples:
  python scripts/audit_submission.py runs/submission.csv --test-parquet data/test.parquet
  python scripts/audit_submission.py runs/submission.csv --expected-ids ids.txt --chemistry
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

EXPECTED_COLUMNS = ("molecule_id", "smiles")
MAX_CANDIDATES = 25


def audit_rows(
    rows: Iterable[Mapping[str, str]],
    expected_ids: Iterable[str] | None = None,
    chemistry_key: Callable[[str], str | None] | None = None,
) -> dict:
    """Return serializable summary with molecule-granular errors/warnings.

    A repeated guess consumes a ranked slot in the official metric, so
    duplicates are warnings rather than silently being deduplicated.
    """
    issues: list[dict] = []
    seen_ids: set[str] = set()
    lengths: list[int] = []
    n_rows = 0
    empty_guesses = 0
    invalid_guesses = 0
    duplicate_guesses = 0

    def issue(level: str, code: str, molecule_id: str, detail: str) -> None:
        issues.append({"severity": level, "code": code,
                       "molecule_id": molecule_id, "detail": detail})

    for row in rows:
        n_rows += 1
        mid = str(row.get("molecule_id") or "").strip()
        raw = row.get("smiles")
        if not mid:
            issue("error", "blank_molecule_id", "", f"row {n_rows}")
            continue
        if mid in seen_ids:
            issue("error", "duplicate_molecule_id", mid, f"row {n_rows}")
        seen_ids.add(mid)
        if raw is None or not str(raw).strip():
            issue("error", "empty_prediction", mid, "no ranked candidates")
            lengths.append(0)
            continue

        guesses = [guess.strip() for guess in str(raw).split(";")]
        lengths.append(len(guesses))
        if len(guesses) > MAX_CANDIDATES:
            issue("error", "over_25", mid, f"{len(guesses)} candidates")
        if len(guesses) == 0:
            issue("error", "no_candidates", mid, "empty ranked list")
        seen_keys: set[str] = set()
        for rank, guess in enumerate(guesses, start=1):
            if not guess:
                empty_guesses += 1
                issue("error", "blank_candidate", mid, f"rank {rank}")
                continue
            key = guess
            if chemistry_key is not None:
                try:
                    key = chemistry_key(guess) or ""
                except Exception:
                    key = ""
                if not key:
                    invalid_guesses += 1
                    issue("error", "invalid_smiles", mid, f"rank {rank}")
                    continue
            if key in seen_keys:
                duplicate_guesses += 1
                issue("warning", "duplicate_2d_candidate" if chemistry_key else
                      "duplicate_smiles", mid, f"rank {rank}")
            seen_keys.add(key)

    coverage = None
    if expected_ids is not None:
        expected = set(str(mid).strip() for mid in expected_ids)
        expected.discard("")
        missing = sorted(expected - seen_ids)
        extras = sorted(seen_ids - expected)
        coverage = {
            "expected_molecules": len(expected),
            "matched_molecules": len(expected & seen_ids),
            "coverage_fraction": len(expected & seen_ids) / len(expected) if expected else None,
            "missing_ids": missing,
            "extra_ids": extras,
        }
        for mid in missing:
            issue("error", "missing_molecule", mid, "no submission row")
        for mid in extras:
            issue("error", "unexpected_molecule", mid, "not in test set")

    counts = Counter(it["severity"] for it in issues)
    return {
        "schema_version": 1,
        "status": "FAIL" if counts["error"] else ("WARN" if counts["warning"] else "PASS"),
        "chemistry_checked": chemistry_key is not None,
        "molecule_rows": n_rows,
        "unique_molecules": len(seen_ids),
        "candidate_count": {
            "min": min(lengths) if lengths else 0,
            "max": max(lengths) if lengths else 0,
            "mean": round(sum(lengths) / len(lengths), 3) if lengths else 0,
        },
        "blank_candidates": empty_guesses,
        "invalid_smiles": invalid_guesses,
        "duplicate_candidates": duplicate_guesses,
        "coverage": coverage,
        "errors": counts["error"],
        "warnings": counts["warning"],
        "issues": issues,
    }


def _chemistry_key(smiles: str) -> str | None:
    # The same canonical key used in our official-metric replica.
    from casmi.chem.canon import canonical_ik14
    return canonical_ik14(smiles)


def _expected_ids(args: argparse.Namespace) -> Sequence[str] | None:
    if args.expected_ids:
        return Path(args.expected_ids).read_text(encoding="utf-8").splitlines()
    if args.test_parquet:
        try:
            import pandas as pd
        except ImportError as exc:
            raise RuntimeError("--test-parquet requires pandas and pyarrow") from exc
        return pd.read_parquet(args.test_parquet, columns=["molecule_id"])["molecule_id"].tolist()
    return None


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("submission", type=Path, help="CASMI submission CSV")
    expected = parser.add_mutually_exclusive_group()
    expected.add_argument("--expected-ids", type=Path, help="one test molecule_id per line")
    expected.add_argument("--test-parquet", type=Path, help="official test.parquet")
    parser.add_argument("--chemistry", action="store_true",
                        help="RDKit tautomer-canonical InChIKey14 validation")
    parser.add_argument("--fail-on-warnings", action="store_true")
    parser.add_argument("--json-out", type=Path, help="save machine-readable audit")
    args = parser.parse_args(argv)

    if args.chemistry:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    try:
        with args.submission.open("r", encoding="utf-8-sig", newline="") as file:
            reader = csv.DictReader(file)
            columns = tuple(reader.fieldnames or ())
            if columns != EXPECTED_COLUMNS:
                print(f"ERROR: expected columns {EXPECTED_COLUMNS}, got {columns}", file=sys.stderr)
                return 2
            result = audit_rows(reader, _expected_ids(args),
                                _chemistry_key if args.chemistry else None)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    rendered = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)
    print(rendered)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered + "\n", encoding="utf-8")
    return 2 if result["errors"] or (args.fail_on_warnings and result["warnings"]) else 0


if __name__ == "__main__":
    raise SystemExit(main())
