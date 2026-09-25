"""Score the current ranker on V-in built against another candidate-pool variant (CASMI_POOL_TAG).

V-in answers are all COCONUT structures, so a bigger pool can only *cost* here (more decoys per
molecule); its benefit — answers that exist only in the new source — is invisible to V-in and has to
be measured on the leaderboard. This script measures the cost, paired against the default pool.

  CASMI_POOL_TAG=_pcl python scripts/eval_pool_variant.py --ranker ranker-v13-clean
"""
import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_TAG, PROJECT_ROOT  # noqa: E402
from casmi.channels.fpmodel import load  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries, mixture_grid  # noqa: E402
from casmi.validate.metric import rank_of_truth, summarize  # noqa: E402

spec = importlib.util.spec_from_file_location("ad3", ROOT / "scripts" / "add_d3_features.py")
ad3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ad3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranker", default="ranker-v13-clean")
    args = ap.parse_args()
    t0 = time.time()
    import lightgbm as lgb
    device = "cuda" if torch.cuda.is_available() else "cpu"
    models = {f: load(ARTIFACT_DIR / "models" / f"fp_fold{f}.pt", device) for f in (0, 1)}
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE)
    store = SpectrumStore(ARTIFACT_DIR / "store")
    vin_q = build_queries(store, store.meta.train_ik14.to_numpy(dtype=object))

    tbl = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob(f"rank_vin{POOL_TAG}_[0-9].parquet"))], ignore_index=True)
    tbl = ad3.add_d3(tbl, [(q, -1) for q in vin_q], models, pool, fps, device)
    rk = [lgb.Booster(model_file=str(p)) for p in sorted((ARTIFACT_DIR / "models").glob(f"{args.ranker}_seed*.txt"))]
    cols = rk[0].feature_name()
    if any(c.startswith("d4_") for c in cols):
        from casmi.channels.fragments import add_d4_table
        tbl = add_d4_table(tbl, {q.qid: q.spectra for q in vin_q}, pool.smiles_of, ARTIFACT_DIR / "frag_cache.pkl")
    score = np.mean([m.predict(tbl[cols]) for m in rk], axis=0)

    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    order = [q.qid for q in vin_q]
    ranks, n_cand = {"C1": {}, "C2": {}}, []
    for (qid, cls), g in tbl.assign(s=score).groupby(["qid", "cls"], sort=False):
        g = g.sort_values("s", ascending=False, kind="stable")
        tk = g.pool_key[g.label == 1].tolist()
        truth = keys.get(tk[0], pool.smiles_of(tk[0])) if tk else None
        ranked = keys.dedupe_ranked(g.pool_key.tolist()[:60], pool.smiles_of, 25)
        ranks[cls][qid] = rank_of_truth([mk for _, mk in ranked], truth) if truth else 0
        if cls == "C2":
            n_cand.append(len(g))
    keys.save()
    tag = f"{args.ranker}{POOL_TAG}"
    report = {"tag": tag, "config": {"ranker": args.ranker, "pool": POOL_FILE}, "n_queries": len(order),
              "median_candidates": float(np.median(n_cand)), "classes": {}}
    print(f"pool {POOL_FILE}: {len(pool):,} structures, median {np.median(n_cand):.0f} candidates per molecule")
    for cls in ("C1", "C2"):
        r = [ranks[cls].get(q, 0) for q in order]
        s = summarize(r)
        report["classes"][cls] = {**s.__dict__, "ranks": r}
        print(f"V-in {cls}: {s}")
    report["classes"]["C3"] = {**summarize([0] * len(order)).__dict__, "ranks": [0] * len(order)}
    report["mixture"] = mixture_grid({c: report["classes"][c]["mrr"] for c in ("C1", "C2", "C3")}).to_dict("records")
    h = hashlib.sha1(json.dumps(report["config"], sort_keys=True).encode()).hexdigest()[:8]
    out = PROJECT_ROOT / "runs" / f"{time.strftime('%Y%m%d-%H%M')}_{tag}_{h}.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"logged runs/{out.name} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
