"""Evidence-safe paired comparison of two CASMI validation runs.

Usage:
  python scripts/compare_runs.py runs/A.json runs/B.json
  python scripts/compare_runs.py runs/OLD_A.json runs/OLD_B.json \
      --allow-unverified-legacy-order

Only records containing matching query IDs, spectra/knockout fingerprints
and fully covered, identical local source-inventory hashes can support an
ACCEPT/REJECT conclusion. Earlier report formats are exploratory. The novelty-mixture intervals
are approximately simultaneous across the fixed grid (not post-selection
corrected). Legacy reports are exploratory.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.paired_bootstrap import paired_difference  # noqa: E402
from casmi.validate.run_alignment import CLASSES, align_runs  # noqa: E402
from casmi.validate.simultaneous_ci import compare_mixtures  # noqa: E402

def rr(ranks):
    r = np.asarray(ranks, dtype=int)
    return np.where(r > 0, 1.0 / np.maximum(r, 1), 0.0)


def main(a_path, b_path, n_boot=4000, seed=0,
         allow_unverified_legacy_order=False):
    if not isinstance(n_boot, int) or n_boot < 100:
        raise ValueError("Use at least 100 bootstrap resamples")
    a, b = (json.loads(Path(p).read_text(encoding="utf-8"))
            for p in (a_path, b_path))
    aligned = align_runs(
        a, b, allow_unverified_legacy_order=allow_unverified_legacy_order,
    )
    if aligned["decision_eligible"]:
        label = "VERIFIED query identities, spectra AND complete local input snapshot"
    elif aligned["verified"]:
        label = ("QUERY-PAIRED ONLY; missing verified complete source snapshot. "
                 "Results are exploratory, never ACCEPT/REJECT")
    else:
        label = "UNVERIFIED legacy ordering; numbers cannot support a scientific verdict"
    print(f"A = {a['tag']}   B = {b['tag']}   n = {aligned['n']}   {label}\n")
    for cls in CLASSES:
        ar, br = aligned["a"][cls], aligned["b"][cls]
        m, lo, hi = paired_difference(ar, br, n_boot=n_boot, seed=seed)
        verdict = ("better" if lo > 0 else "worse" if hi < 0 else "within noise")
        if not aligned["decision_eligible"]:
            verdict = "exploratory only"
        print(f"{cls}: MRR {rr(ar).mean():.4f} -> {rr(br).mean():.4f}   "
              f"delta {m:+.4f} [{lo:+.4f}, {hi:+.4f}]  {verdict}")

    result = compare_mixtures(
        aligned["a"], aligned["b"], n_boot=n_boot, seed=seed,
    )
    verdict = result["verdict"] if aligned["decision_eligible"] else "UNVERIFIED"
    means = [item["delta_mrr"] for item in result["scenarios"]]
    lows = [item["ci_low"] for item in result["scenarios"]]
    highs = [item["ci_high"] for item in result["scenarios"]]
    print(f"\nNovelty-mixture effect range: {min(means):+.4f} to {max(means):+.4f}")
    print(
        f"Approximate simultaneous {result['confidence']:.0%} confidence "
        f"bounds across {len(means)} fixed mixture scenarios: "
        f"lower {min(lows):+.4f} to {max(lows):+.4f}; "
        f"upper {min(highs):+.4f} to {max(highs):+.4f}"
    )
    print("CAUTION: the interval does not correct for repeated model selection "
          "or Kaggle distribution shift.")
    print(f"VERDICT: {verdict}"
          + (" (source/query evidence insufficient; diagnostic numerical pattern: "
             + result["verdict"] + ")" if not aligned["decision_eligible"] else ""))
    return verdict


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_a", type=Path)
    parser.add_argument("run_b", type=Path)
    parser.add_argument("--bootstrap", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--allow-unverified-legacy-order", action="store_true",
        help="Explore old, non-manifested runs; never yields ACCEPT/REJECT",
    )
    args = parser.parse_args()
    try:
        main(args.run_a, args.run_b, n_boot=args.bootstrap, seed=args.seed,
             allow_unverified_legacy_order=args.allow_unverified_legacy_order)
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    cli()
