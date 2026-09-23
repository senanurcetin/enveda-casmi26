"""Training table v3: v1 Class-2 groups + instrument-matched Class-1 groups from enveda-180.

v2 showed cross-library Class-1 matches are still near-duplicates (truth similarity median 1.000,
aggregator libraries re-host the same spectra), while real Class-1 matches are cross-instrument
(V-in truth similarity median 0.848). enveda-180 is timsTOF — the test instrument — so its
structures with references elsewhere reproduce V-in's situation exactly.
Output: artifacts/rank_train_v3_d3.parquet
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

spec = importlib.util.spec_from_file_location("ad3", ROOT / "scripts" / "add_d3_features.py")
ad3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ad3)

device = "cuda" if torch.cuda.is_available() else "cpu"
v1 = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob("rank_train_[0-9].parquet"))], ignore_index=True)
e180 = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob("rank_e180_[0-9].parquet"))], ignore_index=True)
tbl = pd.concat([v1[v1.cls == "C2"], e180[e180.cls == "C1"]], ignore_index=True)

models, pool, fps, _, train_items = ad3.load_everything(device)
store = SpectrumStore(ARTIFACT_DIR / "store")
key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in build_queries(store, key_of_row)}
e180_items = ad3.brd.e180_class1_queries(store, key_of_row, pool, vin_skel, 1200)
tbl = ad3.add_d3(tbl, train_items + e180_items, models, pool, fps, device)
tbl.to_parquet(ARTIFACT_DIR / "rank_train_v3_d3.parquet", index=False)
g = tbl.groupby("cls").qid.nunique()
print(f"v3: {len(tbl):,} rows | C1 groups {g.get('C1', 0)} (enveda-180, cross-instrument) | C2 groups {g.get('C2', 0)}")
truth = tbl[(tbl.cls == "C1") & (tbl.label == 1)].d1_max
print(f"C1 truth d1_max: median {truth.median():.3f}, p25 {truth.quantile(.25):.3f}  (V-in: 0.848 / 0.715)")
