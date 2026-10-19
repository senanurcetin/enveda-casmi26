"""Describe negative-control dependence under verified original source pairing.

Usage:
  python scripts/compare_control_runs.py runs/REFERENCE.json runs/BLANK.json
  python scripts/compare_control_runs.py runs/REFERENCE.json runs/SHUFFLE.json \
      --bootstrap 4000 --seed 11 --output runs/control_diagnostic.json

Unlike compare_runs.py, this tool NEVER issues ACCEPT/REJECT model verdicts
or treats a perturbation as a valid test-time model comparison.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.control_alignment import align_control_runs  # noqa: E402
from casmi.validate.run_alignment import CLASSES  # noqa: E402
from casmi.validate.paired_bootstrap import paired_difference  # noqa: E402
from casmi.validate.simultaneous_ci import compare_mixtures  # noqa: E402


def diagnostic_summary(reference, treatment, *, n_boot=4000, seed=0):
    """Descriptive effects of deliberately removing/corrupting fragment evidence."""
    if type(n_boot) is not int or n_boot < 100:
        raise ValueError("At least 100 bootstrap draws are required")
    if type(seed) is not int or seed < 0:
        raise ValueError("Bootstrap seed must be nonnegative")
    aligned = align_control_runs(reference, treatment)
    classes = {}
    for cls in CLASSES:
        a, b = aligned["reference_ranks"][cls], aligned["treatment_ranks"][cls]
        delta, low, high = paired_difference(a, b, n_boot=n_boot, seed=seed)
        to_rr = lambda x: np.where(np.asarray(x) > 0, 1.0 / np.maximum(x, 1), 0.0)
        classes[cls] = {
            "reference_mrr": float(to_rr(a).mean()),
            "control_mrr": float(to_rr(b).mean()),
            "control_minus_reference": delta,
            "ci_low": low,
            "ci_high": high,
            "ci_type": "pointwise paired bootstrap, descriptive diagnostic only",
        }
    mixed = compare_mixtures(
        aligned["reference_ranks"], aligned["treatment_ranks"],
        n_boot=n_boot, seed=seed,
    )
    mixed.pop("verdict", None)
    mixed["interpretation"] = (
        "Approximate simultaneous band on the prespecified novelty grid. "
        "Descriptive control-minus-reference effect only; NOT an ACCEPT/REJECT "
        "model-improvement or leaderboard significance result."
    )
    return {
        "status": "DIAGNOSTIC_ONLY",
        "verified_original_source_pairing": aligned["verified_source_pairing"],
        "n_molecules": aligned["n_queries"],
        "changed_query_signatures": aligned["changed_queries"],
        "mode": aligned["control_mode"],
        "control_seed": aligned["control_seed"],
        "bootstrap_seed": seed,
        "n_bootstrap": n_boot,
        "direction": "negative_control_minus_unperturbed_reference",
        "classes": classes,
        "mixture": mixed,
        "limitations": (
            "Intentional fragment perturbations measure model reliance on MS/MS "
            "evidence within this selected population; no independence from "
            "training or competition domain shift is implied."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("control", type=Path)
    parser.add_argument("--bootstrap", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        a = json.loads(args.reference.read_text(encoding="utf-8"))
        b = json.loads(args.control.read_text(encoding="utf-8"))
        report = diagnostic_summary(a, b, n_boot=args.bootstrap, seed=args.seed)
        value = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x", encoding="utf-8") as file:
                file.write(value)
        print(value, end="")
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
