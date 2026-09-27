"""Choose LightGBM ranker hyper-parameters by cross-fitting on the training queries.

V-real and V-in are not touched here: each configuration is trained on one skeleton fold and scored
on the other (both directions), on the same table and features as the current ranker. Score =
class-mix-weighted MRR (C1 0.2, C2 0.5, C3 0.3) of the truth's rank within its candidate list.
One seed per fit keeps the sweep cheap; the chosen setting is then trained with all seeds by
train_ranker.py --params and checked once on V-real.

  python scripts/tune_ranker.py --table rank_train_gen2_d3_d4.parquet --like ranker-v16-gen2
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, PROJECT_ROOT, SEED  # noqa: E402
from casmi.rank.ranker import Ranker, RankerConfig  # noqa: E402

MIX = {"C1": 0.2, "C2": 0.5, "C3": 0.3}
BASE = RankerConfig().params
GRID = [
    {},
    {"num_leaves": 63},
    {"num_leaves": 15},
    {"min_data_in_leaf": 200},
    {"learning_rate": 0.03},
    {"feature_fraction": 0.6},
    {"num_leaves": 63, "min_data_in_leaf": 200},
    {"num_leaves": 63, "learning_rate": 0.03, "feature_fraction": 0.6},
    {"lambda_l2": 10.0, "min_data_in_leaf": 100},
    {"lambdarank_truncation_level": 10},
]


def group_mrr(df, score):
    t = df[["qid", "cls", "label"]].assign(s=score)
    t["r"] = t.groupby(["qid", "cls"], sort=False).s.rank(ascending=False, method="first")
    first = t[t.label == 1].groupby(["qid", "cls"]).r.min()
    rr = (1.0 / first).where(first <= 25, 0.0)
    return rr.groupby(level="cls").mean().to_dict()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="rank_train_gen2_d3_d4.parquet")
    ap.add_argument("--like", default="ranker-v16-gen2", help="ranker whose feature list to use")
    args = ap.parse_args()
    import lightgbm as lgb
    t0 = time.time()
    cols = lgb.Booster(model_file=str(ARTIFACT_DIR / "models" / f"{args.like}_seed0.txt")).feature_name()
    tbl = pd.read_parquet(ARTIFACT_DIR / args.table, columns=cols + ["qid", "cls", "label", "fold"])
    tbl = tbl[tbl.groupby(["qid", "cls"]).label.transform("max") > 0].reset_index(drop=True)
    print(f"{len(tbl):,} rows, {tbl.qid.nunique()} queries, {len(cols)} features ({time.time() - t0:.0f}s)", flush=True)
    rng = np.random.default_rng(SEED)
    results = []
    for over in GRID:
        cfg = RankerConfig(seeds=(0,), params={**BASE, **over})
        score = np.zeros(len(tbl))
        iters = []
        for f in (0, 1):
            tr = tbl[tbl.fold != f]
            q = tr.qid.unique()
            es = set(rng.choice(q, len(q) // 5, replace=False))
            r = Ranker(cfg).fit(tr[~tr.qid.isin(es)], tr[tr.qid.isin(es)], cols)
            te = (tbl.fold == f).to_numpy()
            score[te] = r.predict(tbl[te])
            iters += r.best_iters
        m = group_mrr(tbl, score)
        mix = sum(MIX[c] * m.get(c, 0.0) for c in MIX)
        results.append({"params": json.dumps(over), **{c: round(m.get(c, 0), 4) for c in MIX}, "mix": round(mix, 4),
                        "iters": iters})
        print(f"{json.dumps(over):<70} C1 {m.get('C1', 0):.4f}  C2 {m.get('C2', 0):.4f}  C3 {m.get('C3', 0):.4f}  "
              f"mix {mix:.4f}  iters {iters}  ({time.time() - t0:.0f}s)", flush=True)
    res = pd.DataFrame(results).sort_values("mix", ascending=False)
    res.to_csv(PROJECT_ROOT / "runs" / "tune_ranker.csv", index=False)
    print(f"\nbest: {res.iloc[0].params}  (baseline mix {results[0]['mix']:.4f})")


if __name__ == "__main__":
    main()
