"""Train the LambdaRank ranker on simulated training queries, evaluate on V-in, log the run.

  python scripts/train_ranker.py [--no-d3] [--tag name]

Early stopping uses a query-disjoint 20% of the *training* queries — V-in is touched only once,
for the final report, so its numbers stay an honest estimate. The V-in ranks go through the same
metric-key deduplication as a submission, and the run is logged in the format compare_runs.py
reads, so the ranker is compared pairwise against the interim rule on identical queries.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, PROJECT_ROOT, SEED  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.rank.features import feature_columns  # noqa: E402
from casmi.rank.ranker import Ranker  # noqa: E402
from casmi.validate.knockout import mixture_grid  # noqa: E402
from casmi.validate.metric import rank_of_truth, summarize  # noqa: E402


def load(name, d3):
    if d3:
        return pd.read_parquet(ARTIFACT_DIR / f"rank_{name}_d3.parquet")
    return pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob(f"rank_{name}_[0-9].parquet"))], ignore_index=True)


def vin_ranks(tbl, scores, pool, keys, vin_order):
    """1-based rank of the truth after metric-key dedup, per (query, class), in V-in query order."""
    tbl = tbl.assign(score=scores)
    ranks = {"C1": {}, "C2": {}}
    for (qid, cls), g in tbl.groupby(["qid", "cls"], sort=False):
        g = g.sort_values("score", ascending=False, kind="stable")
        truth_keys = g.pool_key[g.label == 1].tolist()
        truth = keys.get(truth_keys[0], pool.smiles_of(truth_keys[0])) if truth_keys else None
        ranked = keys.dedupe_ranked(g.pool_key.tolist()[:60], pool.smiles_of, 25)
        ranks[cls][qid] = rank_of_truth([mk for _, mk in ranked], truth) if truth else 0
    return {cls: [ranks[cls].get(q, 0) for q in vin_order] for cls in ranks}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-d3", action="store_true")
    ap.add_argument("--train-files", default="rank_train_d3.parquet")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--drop-d3", action="store_true", help="train without D3 columns on the same table")
    ap.add_argument("--no-membership", action="store_true",
                    help="drop features that encode how the pool was assembled (train/COCONUT/enveda-180 "
                         "membership, spectrum counts): their class-conditional distribution is set by the "
                         "simulation design, not by chemistry, and the ranker learns them as shortcuts")
    ap.add_argument("--c2-weight", type=float, default=1.0, help="relative weight of Class-2 groups")
    args = ap.parse_args()
    use_d3 = not args.no_d3
    tag = args.tag or ("ranker-d1d2d3" if use_d3 else "ranker-d1d2")
    t0 = time.time()

    train = pd.read_parquet(ARTIFACT_DIR / args.train_files) if use_d3 else load("train", False)
    vin = load("vin", use_d3)
    # np_likeness / has_sugar are only populated for COCONUT structures: NaN itself encodes membership
    MEMBERSHIP = {"prior_in_train", "prior_in_coconut", "prior_enveda180_only", "prior_log_spectra", "d1_nref",
                  "prior_np_likeness", "prior_has_sugar"}
    cols = [c for c in feature_columns(train) if c != "fold" and not (args.drop_d3 and c.startswith("d3_"))
            and not (args.no_membership and c in MEMBERSHIP)]
    qids = train.qid.unique()
    rng = np.random.default_rng(SEED)
    es = set(rng.choice(qids, len(qids) // 5, replace=False))
    tr, va = train[~train.qid.isin(es)], train[train.qid.isin(es)]
    print(f"train {tr.qid.nunique()} queries / {len(tr):,} rows | early-stop {va.qid.nunique()} queries | "
          f"{len(cols)} features", flush=True)

    ranker = Ranker().fit(tr, va, cols, class_weight={"C1": 1.0, "C2": args.c2_weight})
    print(f"best iterations per seed: {ranker.best_iters}")
    print("top features by gain:\n" + ranker.importance().head(12).round(0).to_string())

    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool.parquet"))
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    # V-in query order must match run_validation's, so compare_runs pairs the same molecules
    base = json.loads(sorted((PROJECT_ROOT / "runs").glob("*_d1-d2-analog_*.json"))[-1].read_text())
    order = pd.unique(vin.qid)
    ranks = vin_ranks(vin, ranker.predict(vin), pool, keys, order)
    keys.save()
    ranks["C3"] = [0] * len(order)

    report = {"tag": tag, "config": {"features": cols, "d3": use_d3, "seeds": list(ranker.cfg.seeds),
                                     "c2_weight": args.c2_weight, "train_files": args.train_files, "drop_d3": args.drop_d3, "no_membership": args.no_membership},
              "n_queries": len(order), "classes": {}}
    for cls in ("C1", "C2", "C3"):
        s = summarize(ranks[cls])
        report["classes"][cls] = {**s.__dict__, "ranks": ranks[cls]}
        print(f"{cls}: {s}")
    grid = mixture_grid({c: report["classes"][c]["mrr"] for c in ("C1", "C2", "C3")})
    report["mixture"] = grid.to_dict("records")
    assert len(order) == base["n_queries"], "V-in query set differs from the baseline run"
    h = hashlib.sha1(json.dumps(report["config"], sort_keys=True).encode()).hexdigest()[:8]
    out = PROJECT_ROOT / "runs" / f"{time.strftime('%Y%m%d-%H%M')}_{tag}_{h}.json"
    out.write_text(json.dumps(report, indent=2))
    (ARTIFACT_DIR / "models").mkdir(exist_ok=True)
    for i, m in enumerate(ranker.models):
        m.save_model(str(ARTIFACT_DIR / "models" / f"{tag}_seed{i}.txt"))
    print(f"logged runs/{out.name} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
