"""Training table v2: v1 Class-2 groups + Class-1 groups rebuilt with cross-library knockout.

v1 simulated Class 1 by removing only the query's own spectra, so (a) >= 25% of "Class-1" groups
had no reference left at all (they were Class 2 in disguise) and (b) the rest matched near-duplicate
spectra from the same library and protocol. Test Class-1 matches are cross-instrument.
Output: artifacts/rank_train_v2_d3.parquet
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402

spec = importlib.util.spec_from_file_location("ad3", ROOT / "scripts" / "add_d3_features.py")
ad3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ad3)

device = "cuda" if torch.cuda.is_available() else "cpu"
v1 = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob("rank_train_[0-9].parquet"))], ignore_index=True)
c1x = pd.concat([pd.read_parquet(f) for f in sorted(ARTIFACT_DIR.glob("rank_train_c1x_[0-9].parquet"))], ignore_index=True)
tbl = pd.concat([v1[v1.cls == "C2"], c1x[c1x.cls == "C1"]], ignore_index=True)
models, pool, fps, _, train_items = ad3.load_everything(device)
tbl = ad3.add_d3(tbl, train_items, models, pool, fps, device)
tbl.to_parquet(ARTIFACT_DIR / "rank_train_v2_d3.parquet", index=False)
g = tbl.groupby("cls").qid.nunique()
print(f"v2: {len(tbl):,} rows | C1 groups {g.get('C1', 0)} (cross-library) | C2 groups {g.get('C2', 0)}")
truth = tbl[(tbl.cls == "C1") & (tbl.label == 1)].d1_max
print(f"C1 truth d1_max: median {truth.median():.3f}, p25 {truth.quantile(.25):.3f}  (V-in: 0.848 / 0.715)")
