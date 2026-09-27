"""Recompute the D3 columns of the _gen2 tables with another fingerprint model pair
(models/fp_fold{0,1}{suffix}.pt), writing rank_{set}_gen2{suffix}_d3_d4.parquet.

  python scripts/redo_d3.py --suffix _L
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.channels.fpmodel import load  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

spec = importlib.util.spec_from_file_location("ad3", ROOT / "scripts" / "add_d3_features.py")
ad3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ad3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="_L")
    args = ap.parse_args()
    t0 = time.time()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    models = {f: load(ARTIFACT_DIR / "models" / f"fp_fold{f}{args.suffix}.pt", device) for f in (0, 1)}
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet"))
    fps = np.load(ARTIFACT_DIR / "pool_fp2048_pcl.npy")
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    brd = ad3.brd
    excl = vin_skel | brd.vreal_skeletons()
    items = {"vin": [(q, -1) for q in vin],
             "vreal": brd.vreal_queries(store, key_of_row, vin_skel)[0],
             "train": brd.training_queries(store, key_of_row, pool, excl, 3000)
             + brd.e180_class1_queries(store, key_of_row, pool, excl, 1200)}
    for name, its in items.items():
        out = ARTIFACT_DIR / f"rank_{name}_gen2{args.suffix}_d3_d4.parquet"
        if out.exists():
            continue
        tbl = pd.read_parquet(ARTIFACT_DIR / f"rank_{name}_gen2_d3_d4.parquet")
        tbl = tbl[[c for c in tbl.columns if not c.startswith("d3_")]]
        g = tbl[tbl.gen == 1]
        pool.register_generated(g.pool_key.to_numpy(), g.gen_smiles.to_numpy())
        tbl = ad3.add_d3(tbl, its, models, pool, fps, device)
        tbl.to_parquet(out, index=False)
        print(f"{out.name}: {len(tbl):,} rows ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
