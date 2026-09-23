"""Paired comparison of two validation runs on the same queries.

Usage: python scripts/compare_runs.py runs/A.json runs/B.json

Per class: ΔMRR (B − A) with a paired bootstrap 95% CI.
Mixture: for every class-share point (f1, f3) the expected-LB difference is bootstrapped *jointly*
over queries (the same molecules appear in every simulated class, so their deltas are resampled
together). Decision rule:
  ACCEPT        the lower CI bound is > 0 at every grid point
  REJECT        the upper CI bound is < 0 at every grid point
  NOT PROVEN    anything else — the change may be real but this validation cannot show it
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.metric import paired_difference  # noqa: E402

GRID = [(f1, f3) for f1 in (0.15, 0.20, 0.25, 0.30) for f3 in (0.2, 0.3, 0.4)]


def rr(ranks):
    r = np.asarray(ranks)
    return np.where(r > 0, 1.0 / np.maximum(r, 1), 0.0)


def main(a_path, b_path, n_boot=4000, seed=0):
    a, b = (json.loads(Path(p).read_text()) for p in (a_path, b_path))
    assert a["n_queries"] == b["n_queries"], "runs are not on the same queries"
    print(f"A = {a['tag']}   B = {b['tag']}   n = {a['n_queries']}\n")
    d = {}
    for cls in ("C1", "C2", "C3"):
        d[cls] = rr(b["classes"][cls]["ranks"]) - rr(a["classes"][cls]["ranks"])
        m, lo, hi = paired_difference(a["classes"][cls]["ranks"], b["classes"][cls]["ranks"])
        verdict = "better" if lo > 0 else "worse" if hi < 0 else "within noise"
        print(f"{cls}: MRR {a['classes'][cls]['mrr']:.4f} -> {b['classes'][cls]['mrr']:.4f}   "
              f"Δ {m:+.4f} [{lo:+.4f}, {hi:+.4f}]  {verdict}")

    rng = np.random.default_rng(seed)
    n = len(d["C1"])
    idx = rng.integers(0, n, size=(n_boot, n))
    lows, highs, means = [], [], []
    for f1, f3 in GRID:
        per_query = f1 * d["C1"] + (1 - f1 - f3) * d["C2"] + f3 * d["C3"]
        boots = per_query[idx].mean(axis=1)
        means.append(per_query.mean())
        lows.append(np.quantile(boots, .025))
        highs.append(np.quantile(boots, .975))
    verdict = "ACCEPT" if min(lows) > 0 else "REJECT" if max(highs) < 0 else "NOT PROVEN"
    print(f"\nexpected-LB Δ over the class-share grid: point {min(means):+.4f} … {max(means):+.4f} | "
          f"95% CI lower bounds {min(lows):+.4f} … {max(lows):+.4f}  ->  {verdict}")
    return verdict


if __name__ == "__main__":
    main(*sys.argv[1:3])
