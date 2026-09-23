"""Paired comparison of two validation runs on the same queries.

Usage: python scripts/compare_runs.py runs/A.json runs/B.json
Prints ΔMRR (B − A) per class with a bootstrap 95% CI, and the Δ of the mixture score over the
class-share grid. A change is accepted only if the mixture Δ is positive across the whole grid.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.knockout import mixture  # noqa: E402
from casmi.validate.metric import paired_difference  # noqa: E402


def main(a_path, b_path):
    a, b = (json.loads(Path(p).read_text()) for p in (a_path, b_path))
    assert a["n_queries"] == b["n_queries"], "runs are not on the same queries"
    print(f"A = {a['tag']}   B = {b['tag']}   n = {a['n_queries']}\n")
    deltas = {}
    for cls in ("C1", "C2", "C3"):
        d, lo, hi = paired_difference(a["classes"][cls]["ranks"], b["classes"][cls]["ranks"])
        deltas[cls] = d
        verdict = "better" if lo > 0 else "worse" if hi < 0 else "within noise"
        print(f"{cls}: MRR {a['classes'][cls]['mrr']:.4f} -> {b['classes'][cls]['mrr']:.4f}   "
              f"Δ {d:+.4f} [{lo:+.4f}, {hi:+.4f}]  {verdict}")
    grid = [(f1, f3, mixture(deltas, f1, f3)) for f1 in (0.15, 0.2, 0.25, 0.3) for f3 in (0.2, 0.3, 0.4)]
    lo, hi = min(g[2] for g in grid), max(g[2] for g in grid)
    print(f"\nexpected-LB Δ over the class-share grid: {lo:+.4f} … {hi:+.4f}  "
          f"-> {'ACCEPT' if lo > 0 else 'REJECT' if hi < 0 else 'MIXED — depends on the class mix'}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
