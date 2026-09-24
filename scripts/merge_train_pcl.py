"""v13 training table rebuilt on the PubChemLite pool (CASMI_POOL_TAG=_pcl):
Class 2 and v1-style Class 1 from the NP libraries + instrument-matched Class 1 from enveda-180,
with out-of-fold D3 scored against the _pcl pool fingerprints.
Output: artifacts/rank_train_v13_pcl_d3.parquet
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_TAG  # noqa: E402
from casmi.channels.fpmodel import load  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

assert POOL_TAG == "_pcl", "run with CASMI_POOL_TAG=_pcl"
spec = importlib.util.spec_from_file_location("ad3", ROOT / "scripts" / "add_d3_features.py")
ad3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ad3)

device = "cuda" if torch.cuda.is_available() else "cpu"
tr = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob("rank_train_pcl_[0-9].parquet"))], ignore_index=True)
e180 = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob("rank_e180_pcl_[0-9].parquet"))], ignore_index=True)
tbl = pd.concat([tr, e180[e180.cls == "C1"]], ignore_index=True)

models = {f: load(ARTIFACT_DIR / "models" / f"fp_fold{f}.pt", device) for f in (0, 1)}
pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
fps = np.load(ARTIFACT_DIR / POOL_FP_FILE)
store = SpectrumStore(ARTIFACT_DIR / "store")
key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in build_queries(store, key_of_row)}
items = ad3.brd.training_queries(store, key_of_row, pool, vin_skel, 1500) + \
        ad3.brd.e180_class1_queries(store, key_of_row, pool, vin_skel, 1200)
tbl = ad3.add_d3(tbl, items, models, pool, fps, device)
tbl.to_parquet(ARTIFACT_DIR / "rank_train_v13_pcl_d3.parquet", index=False)
print(f"v13 on _pcl: {len(tbl):,} rows |", tbl.groupby("cls").qid.nunique().to_dict())
