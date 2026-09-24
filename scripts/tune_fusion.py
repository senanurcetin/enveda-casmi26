"""Combine the learned ranker with the interim rule (max of direct and propagated evidence).

The ranker wins Class 2 (0.662 vs 0.252) but loses Class 1 (0.839 vs 0.908) to the simple rule.
Two ways to keep both:
  fusion   score = (1 - w) · RRF(ranker rank) + w · RRF(rule rank)
  gate     a candidate whose direct-match similarity is >= tau goes first; the rest by ranker

w and tau are chosen on OUT-OF-FOLD ranker scores of the *training* queries (the ranker is
cross-fitted on the skeleton folds), maximising the worst-case mixture score over the class-share
grid. V-in is scored once, with the chosen setting, for the paired comparison.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE, PROJECT_ROOT, SEED  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.rank.ranker import Ranker  # noqa: E402
from casmi.validate.knockout import mixture_grid  # noqa: E402
from casmi.validate.metric import rank_of_truth, summarize  # noqa: E402

TRAIN_FILE = "rank_train_v13_d3.parquet"
GRID = [(f1, f3) for f1 in (0.15, 0.20, 0.25, 0.30) for f3 in (0.2, 0.3, 0.4)]


def columns_of(tag):
    import lightgbm as lgb
    return lgb.Booster(model_file=str(ARTIFACT_DIR / "models" / f"{tag}_seed0.txt")).feature_name()


def rule_score(t):
    return np.maximum(t.d1_max.to_numpy(), t.an_p4.to_numpy()) + 1e-3 * t.an_p2.to_numpy() - 1e-6 * t.abs_ppm.to_numpy()


def combined(t, ranker_score, mode, param):
    g = t.assign(rs=ranker_score, ru=rule_score(t)).groupby(["qid", "cls"], sort=False)
    r_rank = g.rs.rank(ascending=False, method="first").to_numpy()
    u_rank = g.ru.rank(ascending=False, method="first").to_numpy()
    if mode == "fusion":
        return (1 - param) / (60 + r_rank) + param / (60 + u_rank)
    # gate: candidates with direct similarity >= tau first (ordered by it), then the ranker order
    gated = t.d1_max.to_numpy() >= param
    return np.where(gated, 10.0 + t.d1_max.to_numpy(), 1.0 / (1.0 + r_rank))


def mrr_by_class(t, score):
    """Fast, dedup-free reciprocal rank on the training table (tuning only)."""
    t = t.assign(s=score)
    t["r"] = t.groupby(["qid", "cls"], sort=False).s.rank(ascending=False, method="first")
    first = t[t.label == 1].groupby(["qid", "cls"]).r.min()
    rr = (1.0 / first).where(first <= 25, 0.0)
    return rr.groupby(level="cls").mean().to_dict()


def worst_case(mrr):
    return min(f1 * mrr.get("C1", 0) + (1 - f1 - f3) * mrr.get("C2", 0) for f1, f3 in GRID)


def main():
    t0 = time.time()
    train = pd.read_parquet(ARTIFACT_DIR / TRAIN_FILE)
    train = train[train.groupby(["qid", "cls"]).label.transform("max") > 0].reset_index(drop=True)
    cols = columns_of("ranker-v13-clean")

    # cross-fitted (out-of-fold) ranker scores on the training queries
    oof = np.zeros(len(train))
    rng = np.random.default_rng(SEED)
    for f in (0, 1):
        tr, te = train[train.fold != f], train.fold == f
        qids = tr.qid.unique()
        es = set(rng.choice(qids, len(qids) // 5, replace=False))
        r = Ranker().fit(tr[~tr.qid.isin(es)], tr[tr.qid.isin(es)], cols)
        oof[te.to_numpy()] = r.predict(train[te])
    print(f"cross-fitted ranker in {time.time() - t0:.0f}s", flush=True)

    rows = [("ranker only", "fusion", 0.0), ("rule only", "fusion", 1.0)]
    rows += [(f"fusion w={w}", "fusion", w) for w in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6)]
    rows += [(f"gate tau={t}", "gate", t) for t in (0.7, 0.75, 0.8, 0.85, 0.9, 0.95)]
    results = []
    for name, mode, p in rows:
        m = mrr_by_class(train, combined(train, oof, mode, p))
        results.append((name, mode, p, m.get("C1", 0), m.get("C2", 0), worst_case(m)))
    res = pd.DataFrame(results, columns=["setting", "mode", "param", "C1", "C2", "worst_case_mix"])
    print("\nout-of-fold, training queries:\n" + res.round(4).to_string(index=False))
    best = res.iloc[res.worst_case_mix.idxmax()]
    print(f"\nchosen on training OOF: {best.setting}")

    # V-in: once, with the full ranker already trained on all training queries
    import lightgbm as lgb
    vin = pd.read_parquet(ARTIFACT_DIR / "rank_vin_d3.parquet")
    models = [lgb.Booster(model_file=str(p)) for p in sorted((ARTIFACT_DIR / "models").glob("ranker-v13-clean_seed*.txt"))]
    rs = np.mean([m.predict(vin[cols]) for m in models], axis=0)
    score = combined(vin, rs, best["mode"], best.param)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    order = pd.unique(vin.qid)
    ranks = {"C1": {}, "C2": {}}
    for (qid, cls), g in vin.assign(s=score).groupby(["qid", "cls"], sort=False):
        g = g.sort_values("s", ascending=False, kind="stable")
        tk = g.pool_key[g.label == 1].tolist()
        truth = keys.get(tk[0], pool.smiles_of(tk[0])) if tk else None
        ranked = keys.dedupe_ranked(g.pool_key.tolist()[:60], pool.smiles_of, 25)
        ranks[cls][qid] = rank_of_truth([mk for _, mk in ranked], truth) if truth else 0
    keys.save()
    tag = "fusion-v13-clean"
    report = {"tag": tag, "config": {"base": "ranker-v13-clean", "mode": best["mode"], "param": float(best.param),
                                     "chosen_on": "training OOF worst-case mixture"},
              "n_queries": len(order), "classes": {}, "tuning": res.to_dict("records")}
    for cls in ("C1", "C2"):
        rk = [ranks[cls].get(q, 0) for q in order]
        s = summarize(rk)
        report["classes"][cls] = {**s.__dict__, "ranks": rk}
        print(f"V-in {cls}: {s}")
    report["classes"]["C3"] = {**summarize([0] * len(order)).__dict__, "ranks": [0] * len(order)}
    report["mixture"] = mixture_grid({c: report["classes"][c]["mrr"] for c in ("C1", "C2", "C3")}).to_dict("records")
    h = hashlib.sha1(json.dumps(report["config"], sort_keys=True).encode()).hexdigest()[:8]
    out = PROJECT_ROOT / "runs" / f"{time.strftime('%Y%m%d-%H%M')}_{tag}_{h}.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"logged runs/{out.name} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
