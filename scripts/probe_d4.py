"""D4 variants, judged alone: MRR of the truth when a Class-2 candidate list (pool candidates) is
ordered by one D4 statistic. Cheap screen before rebuilding ranker tables.

  python scripts/probe_d4.py [--n 250]
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.channels.fragments import FragmentScorer  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

spec = importlib.util.spec_from_file_location("brd", ROOT / "scripts" / "build_rank_data.py")
brd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(brd)

VARIANTS = {
    "base": {},
    "tol 0.005/5ppm": {"tol_da": 0.005, "tol_ppm": 5.0},
    "tol 0.02/20ppm": {"tol_da": 0.02, "tol_ppm": 20.0},
    "H ±1": {"h_max": 1},
    "H ±3": {"h_max": 3},
    "weight raw": {"weight": "raw"},
    "weight log": {"weight": "log"},
    "weight count": {"weight": "count"},
    "depth 3": {"depth": 3},
}


def mrr(v, lab):
    r = (v > v[lab].max()).sum() + 1
    return 1.0 / r if r <= 25 else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=250)
    args = ap.parse_args()
    t0 = time.time()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet"))
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    spectra = {q.qid: q.spectra for q in vin}
    spectra.update({q.qid: q.spectra for q, _ in brd.vreal_queries(store, key_of_row, vin_skel)[0]})
    cols = ["qid", "cls", "pool_key", "label", "gen"]
    sets = {}
    for name in ("vin", "vreal"):
        t = pd.read_parquet(ARTIFACT_DIR / f"rank_{name}_gen2_0.parquet", columns=cols + (["stratum"] if name == "vreal" else []))
        t = t[(t.cls == "C2") & (t.gen == 0)]
        if name == "vreal":
            t = t[t.stratum == "coco"]
        t = t[t.groupby("qid").label.transform("max") > 0]
        q = t.qid.unique()[: args.n]
        sets[name] = t[t.qid.isin(q)]
    print({k: v.qid.nunique() for k, v in sets.items()}, f"({time.time() - t0:.0f}s)", flush=True)
    for vname, kw in VARIANTS.items():
        sc = FragmentScorer(pool.smiles_of, **kw)
        out = []
        for name, t in sets.items():
            r = {"int_max": [], "pk_max": [], "int": []}
            for qid, g in t.groupby("qid", sort=False):
                d4 = sc.score(spectra[qid], g.pool_key.to_numpy())
                lab = g.label.to_numpy() == 1
                r["int"].append(mrr(d4[:, 0], lab))
                r["int_max"].append(mrr(d4[:, 1], lab))
                r["pk_max"].append(mrr(d4[:, 3], lab))
            out.append(f"{name}: " + " ".join(f"{k} {np.mean(v):.3f}" for k, v in r.items()))
        print(f"{vname:<16} " + " | ".join(out) + f"  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
