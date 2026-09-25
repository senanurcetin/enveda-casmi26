"""Score rankers on the V-real / V-in tables that carry Class-3 hedge candidates.

  python scripts/eval_gen.py --set vreal --scorers ranker-v13-clean-pcl:nogen,ranker-v15-gen

A scorer suffixed ':nogen' ranks the table with the generated rows removed — exactly what the
submission does today — so every comparison is against the production behaviour on the same
queries. Ranks go through metric-key deduplication; Class 1/2/3 ranks are logged per query in the
compare_runs.py format (a query without a group in some class scores 0 there under every scorer).
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
from casmi.config import ARTIFACT_DIR, PROJECT_ROOT  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.validate.knockout import mixture_grid  # noqa: E402
from casmi.validate.metric import rank_of_truth, summarize  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="vreal", choices=["vreal", "vin"])
    ap.add_argument("--scorers", required=True)
    ap.add_argument("--stratum", default="", help="restrict V-real to one stratum (coco / pcl)")
    args = ap.parse_args()
    t0 = time.time()
    import lightgbm as lgb
    tbl = pd.read_parquet(ARTIFACT_DIR / f"rank_{args.set}_gen_d3_d4.parquet")
    if args.stratum:
        tbl = tbl[tbl.stratum == args.stratum]
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet"))
    g = tbl[tbl.gen == 1]
    pool.register_generated(g.pool_key.to_numpy(), g.gen_smiles.to_numpy())
    smiles = pd.read_parquet(ARTIFACT_DIR / "structures.parquet", columns=["pool_key", "smiles"])
    smiles = dict(zip(smiles.pool_key, smiles.smiles))
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    order = list(pd.unique(tbl.qid))
    for spec in args.scorers.split(","):
        name, _, mode = spec.partition(":")
        t = tbl[tbl.gen == 0] if mode == "nogen" else tbl
        rk = [lgb.Booster(model_file=str(p)) for p in sorted((ARTIFACT_DIR / "models").glob(f"{name}_seed*.txt"))]
        cols = rk[0].feature_name()
        s = np.mean([m.predict(t[cols]) for m in rk], axis=0)
        got = {}
        for (qid, cls), grp in t.assign(s=s).groupby(["qid", "cls"], sort=False):
            grp = grp.sort_values("s", ascending=False, kind="stable")
            key = qid.split(":", 1)[1]
            truth = keys.get(key, smiles[key])
            ranked = keys.dedupe_ranked(grp.pool_key.tolist()[:60], pool.smiles_of, 25)
            got[(qid, cls)] = rank_of_truth([mk for _, mk in ranked], truth)
        report = {"tag": f"{args.set}{args.stratum}_{spec.replace(':', '-')}", "config": {"set": args.set, "scorer": spec, "stratum": args.stratum},
                  "n_queries": len(order), "classes": {}}
        line = []
        for cls in ("C1", "C2", "C3"):
            r = [got.get((q, cls), 0) for q in order]
            sm = summarize(r)
            report["classes"][cls] = {**sm.__dict__, "ranks": r}
            line.append(f"{cls} {sm.mrr:.4f}")
        report["mixture"] = mixture_grid({c: report["classes"][c]["mrr"] for c in ("C1", "C2", "C3")}).to_dict("records")
        h = hashlib.sha1(json.dumps(report["config"], sort_keys=True).encode()).hexdigest()[:8]
        out = PROJECT_ROOT / "runs" / f"{time.strftime('%Y%m%d-%H%M')}_{report['tag']}_{h}.json"
        out.write_text(json.dumps(report))
        print(f"{spec:<32} " + "  ".join(line) + f"   -> runs/{out.name}", flush=True)
    keys.save()
    print(f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
