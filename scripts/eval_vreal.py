"""Score (pool variant × scorer) combinations on V-real and check them against the leaderboard.

V-real tables are built per pool with `build_rank_data.py --set vreal` (CASMI_POOL_TAG selects the
pool). Here D3 is joined out-of-fold, every scorer ranks every (query, class) group, and the
reciprocal ranks are logged per stratum (pcl / coco) and class (C1 / C2).

  python scripts/eval_vreal.py                         # every leaderboard-known configuration
  python scripts/eval_vreal.py --only _pcl:ranker-x    # one configuration

The leaderboard-known configurations are the calibration targets: a validation that deserves to
drive decisions must order them as the leaderboard did.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, PROJECT_ROOT  # noqa: E402
from casmi.channels.fpmodel import load  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402
from casmi.validate.metric import rank_of_truth  # noqa: E402


def _module(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / file)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


brd = _module("brd", "build_rank_data.py")
ad3 = _module("ad3", "add_d3_features.py")

# (pool tag, scorer) -> public LB, from EXPERIMENTS.md
LEADERBOARD = {
    ("", "rule"): 0.190,
    ("", "ranker-d1d2d3"): 0.245,
    ("", "ranker-v13-clean"): 0.258,
    ("_pcl", "ranker-v13-clean"): 0.267,
    ("_pcl", "ranker-v13-clean-pcl"): 0.273,
    ("_pclx", "ranker-v13-clean-pcl"): 0.268,
    ("_pclp", "ranker-v13-clean-pcl"): 0.252,
}


def rule_score(t):
    """The interim rule of BaselineEngine.rank as one sortable score."""
    return (np.maximum(t.d1_max.to_numpy(), t.an_p4.to_numpy()) + 1e-3 * t.an_p2.to_numpy()
            - 1e-6 * np.abs(t.ppm.to_numpy()))


def table_for(tag, items, models, device):
    cached = ARTIFACT_DIR / f"rank_vreal{tag}_d3.parquet"
    if cached.exists():
        return pd.read_parquet(cached)
    tbl = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob(f"rank_vreal{tag}_[0-9].parquet"))],
                    ignore_index=True)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / f"pool{tag}.parquet"))
    fps = np.load(ARTIFACT_DIR / f"pool_fp2048{tag}.npy")
    tbl = ad3.add_d3(tbl, items, models, pool, fps, device)
    tbl.to_parquet(cached, index=False)
    return tbl


def with_d4(tag, tbl, items):
    cached = ARTIFACT_DIR / f"rank_vreal{tag}_d3_d4.parquet"
    if cached.exists():
        return pd.read_parquet(cached)
    from casmi.channels.fragments import add_d4_table
    out = add_d4_table(tbl, {q.qid: q.spectra for q, _ in items}, smiles_any, ARTIFACT_DIR / "frag_cache.pkl")
    out.to_parquet(cached, index=False)
    return out


def score(tbl, scorer):
    if scorer == "rule":
        return rule_score(tbl)
    import lightgbm as lgb
    rk = [lgb.Booster(model_file=str(p)) for p in sorted((ARTIFACT_DIR / "models").glob(f"{scorer}_seed*.txt"))]
    cols = rk[0].feature_name()
    return np.mean([m.predict(tbl[cols]) for m in rk], axis=0)


def ranks_of(tbl, s, items, stratum_of, keys, smiles):
    """Reciprocal-rank inputs per (stratum, class); a group with no candidates or no truth scores 0."""
    got = {}
    for (qid, cls), g in tbl.assign(s=s).groupby(["qid", "cls"], sort=False):
        g = g.sort_values("s", ascending=False, kind="stable")
        key = qid.split(":", 1)[1]
        truth = keys.get(key, smiles[key])
        cand_smiles = dict(zip(g.pool_key, g.smiles)) if "smiles" in g else None
        ranked = keys.dedupe_ranked(g.pool_key.tolist()[:60], lambda k: cand_smiles[k] if cand_smiles else smiles_any(k), 25)
        got[(qid, cls)] = rank_of_truth([mk for _, mk in ranked], truth)
    out = {}
    for q, _ in items:
        st = stratum_of[q.pool_key]
        for cls in ("C1", "C2"):
            out.setdefault(f"{st}-{cls}", []).append(got.get((q.qid, cls), 0))
    return out


_POOL_SMILES = {}


def uses_d4(scorer):
    import lightgbm as lgb
    p = sorted((ARTIFACT_DIR / "models").glob(f"{scorer}_seed*.txt"))[0]
    return any(c.startswith("d4_") for c in lgb.Booster(model_file=str(p)).feature_name())


def smiles_any(k):
    return _POOL_SMILES[k]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="tag:scorer[,tag:scorer...]")
    args = ap.parse_args()
    t0 = time.time()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    models = {f: load(ARTIFACT_DIR / "models" / f"fp_fold{f}.pt", device) for f in (0, 1)}
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    base_pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool.parquet"))
    vin_skel = {skeleton_key(base_pool.smiles_of(q.pool_key)) for q in build_queries(store, key_of_row)}
    items, stratum_of = brd.vreal_queries(store, key_of_row, vin_skel)
    smiles = brd.structure_smiles()
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    configs = ([tuple(c.split(":")) for c in args.only.split(",")] if args.only else list(LEADERBOARD))
    print(f"V-real: {len(items)} queries ({sum(stratum_of[q.pool_key] == 'pcl' for q, _ in items)} pcl) "
          f"loaded in {time.time() - t0:.0f}s", flush=True)
    rows = []
    for tag in dict.fromkeys(t for t, _ in configs):
        tbl = table_for(tag, items, models, device)
        pool_tbl = pd.read_parquet(ARTIFACT_DIR / f"pool{tag}.parquet", columns=["pool_key", "smiles"])
        _POOL_SMILES.clear()
        _POOL_SMILES.update(zip(pool_tbl.pool_key, pool_tbl.smiles))
        for t_, scorer in configs:
            if t_ != tag:
                continue
            t_use = tbl
            if scorer != "rule" and uses_d4(scorer):
                t_use = with_d4(tag, tbl, items)
            r = ranks_of(t_use, score(t_use, scorer), items, stratum_of, keys, smiles)
            mrr = {k: float(np.mean([1 / x if x else 0 for x in v])) for k, v in r.items()}
            lb = LEADERBOARD.get((tag, scorer))
            rows.append({"pool": tag or "(default)", "scorer": scorer, **{k: round(v, 4) for k, v in mrr.items()},
                         "LB": lb})
            report = {"tag": f"vreal_{scorer}{tag}", "config": {"set": "vreal", "pool": tag, "scorer": scorer},
                      "n_queries": len(items), "strata": {k: {"mrr": mrr[k], "ranks": v} for k, v in r.items()},
                      "leaderboard": lb}
            h = hashlib.sha1(json.dumps(report["config"], sort_keys=True).encode()).hexdigest()[:8]
            (PROJECT_ROOT / "runs" / f"{time.strftime('%Y%m%d-%H%M')}_vreal_{scorer}{tag}_{h}.json").write_text(
                json.dumps(report, indent=1))
            print(f"  {tag or '(default)':<10} {scorer:<22} " + "  ".join(f"{k} {v:.3f}" for k, v in mrr.items())
                  + (f"   LB {lb:.3f}" if lb else ""), flush=True)
    keys.save()
    res = pd.DataFrame(rows)
    res.to_csv(PROJECT_ROOT / "runs" / "vreal_summary.csv", index=False)
    print(f"\n{res.to_string(index=False)}\n({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
