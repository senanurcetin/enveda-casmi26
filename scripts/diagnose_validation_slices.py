"""Descriptive error-slice analysis on two signed CASMI validation run JSONs.

Usage:
  python scripts/diagnose_validation_slices.py runs/A.json runs/B.json \
    --min-n 15 --output runs/error_slices.json

No hypothesis test, model decision or competition score is produced.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.error_slices import diagnostic_slices  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--min-n", type=int, default=15)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        reference = json.loads(args.reference.read_text(encoding="utf-8"))
        candidate = json.loads(args.candidate.read_text(encoding="utf-8"))
        report = diagnostic_slices(reference, candidate, min_n=args.min_n)
        payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x", encoding="utf-8") as file:
                file.write(payload)
        print(payload, end="")
    except (ValueError, KeyError, OSError, TypeError, json.JSONDecodeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
