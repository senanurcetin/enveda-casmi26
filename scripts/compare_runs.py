"""Evidence-safe paired comparison of two CASMI validation runs.

Usage:
  python scripts/compare_runs.py runs/A.json runs/B.json
  python scripts/compare_runs.py runs/OLD_A.json runs/OLD_B.json \
      --allow-unverified-legacy-order

Only records containing matching query IDs and spectra/knockout fingerprints
can support an ACCEPT/REJECT conclusion. Legacy reports are exploratory.
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

GRID = [(f1, f3) for f1 in (0.15, 0.20, 0.25, 0.30)
        for f3 in (0.2, 0.3, 0.4)]


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
    label = "VERIFIED paired identities + input spectra" if aligned["verified"] else (
        "UNVERIFIED legacy ordering; numbers cannot support a scientific verdict"
    )
    print(f"A = {a['tag']}   B = {b['tag']}   n = {aligned['n']}   {label}\n")
    d = {}
    for cls in CLASSES:
        ar, br = aligned["a"][cls], aligned["b"][cls]
        d[cls] = rr(br) - rr(ar)
        m, lo, hi = paired_difference(ar, br, n_boot=n_boot, seed=seed)
        verdict = "better" if lo > 0 else "worse" if hi < 0 else "within noise"
        print(f"{cls}: MRR {rr(ar).mean():.4f} -> {rr(br).mean():.4f}   "
              f"delta {m:+.4f} [{lo:+.4f}, {hi:+.4f}]  {verdict}")

    rng = np.random.default_rng(seed)
    n = aligned["n"]
    idx = rng.integers(0, n, size=(n_boot, n))
    lows, highs, means = [], [], []
    for f1, f3 in GRID:
        per_query = f1 * d["C1"] + (1 - f1 - f3) * d["C2"] + f3 * d["C3"]
        boots = per_query[idx].mean(axis=1)
        means.append(float(per_query.mean()))
        lows.append(float(np.quantile(boots, 0.025)))
        highs.append(float(np.quantile(boots, 0.975)))

    statistical = ("ACCEPT" if min(lows) > 0 else
                   "REJECT" if max(highs) < 0 else "NOT PROVEN")
    verdict = statistical if aligned["verified"] else "UNVERIFIED"
    print(f"\nExpected mixture delta (C1/C2/C3 grid): "
          f"{min(means):+.4f} to {max(means):+.4f}; "
          f"95% CI lower bounds {min(lows):+.4f} to {max(lows):+.4f}.")
    print(f"VERDICT: {verdict}"
          + (" (legacy numerical result: " + statistical + ")" if not aligned["verified"] else ""))
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
