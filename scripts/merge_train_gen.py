"""Tables with Class-3 hedge candidates (built with --generate on pool _pcl) -> D3 + D4 joined.

  train: rank_train_gen_* (C1 own-library, C2, C3) + rank_e180_gen_* (C1)  -> rank_train_gen_d3_d4.parquet
  vreal: rank_vreal_gen_0                                                 -> rank_vreal_gen_d3_d4.parquet
  vin:   rank_vin_gen_0                                                   -> rank_vin_gen_d3_d4.parquet
"""
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_TAG  # noqa: E402
from casmi.channels.fpmodel import load  # noqa: E402
from casmi.channels.fragments import add_d4_table  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

assert POOL_TAG == "_pcl", "run with CASMI_POOL_TAG=_pcl"
spec = importlib.util.spec_from_file_location("ad3", ROOT / "scripts" / "add_d3_features.py")
ad3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ad3)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="_gen", help="table suffix: rank_{set}{suffix}_*.parquet")
    ap.add_argument("--n-train", type=int, default=1500)
    args = ap.parse_args()
    sfx = args.suffix
    t0 = time.time()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    models = {f: load(ARTIFACT_DIR / "models" / f"fp_fold{f}.pt", device) for f in (0, 1)}
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE)
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    brd = ad3.brd
    excl = vin_skel | brd.vreal_skeletons()
    sets = {
        "vin": ([(q, -1) for q in vin], [f"rank_vin{sfx}_0.parquet"]),
        "vreal": (brd.vreal_queries(store, key_of_row, vin_skel)[0], [f"rank_vreal{sfx}_0.parquet"]),
        "train": (brd.training_queries(store, key_of_row, pool, excl, args.n_train)
                  + brd.e180_class1_queries(store, key_of_row, pool, excl, 1200),
                  sorted(f.name for f in ARTIFACT_DIR.glob(f"rank_train{sfx}_[0-9].parquet")) + [f"rank_e180{sfx}_0.parquet"]),
    }
    for name, (items, files) in sets.items():
        out = ARTIFACT_DIR / f"rank_{name}{sfx}_d3_d4.parquet"
        if out.exists():
            print(f"{out.name} exists")
            continue
        tbl = pd.concat([pd.read_parquet(ARTIFACT_DIR / f) for f in files], ignore_index=True)
        if name == "train":
            tbl = tbl[~((tbl.qid.str.startswith("enveda-180:")) & (tbl.cls != "C1"))]
        g = tbl[tbl.gen == 1]
        pool.register_generated(g.pool_key.to_numpy(), g.gen_smiles.to_numpy())
        tbl = ad3.add_d3(tbl, items, models, pool, fps, device)
        tbl = add_d4_table(tbl, {q.qid: q.spectra for q, _ in items}, pool.smiles_of, ARTIFACT_DIR / "frag_cache.pkl")
        tbl.to_parquet(out, index=False)
        print(f"{out.name}: {len(tbl):,} rows, groups {tbl.groupby('cls').qid.nunique().to_dict()}, "
              f"generated truths {int(tbl[(tbl.gen == 1)].label.sum())} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
