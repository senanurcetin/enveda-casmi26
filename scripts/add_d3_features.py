"""Join out-of-fold D3 (spectrum → fingerprint) features onto the ranker tables.

The D3 score of a candidate depends only on the query spectra and the candidate fingerprint —
not on the library knockout — so it can be added after the tables are built. What matters is
that no query is scored by a model that trained on its structure:
  train queries in fold f  ->  model trained on fold 1-f
  V-in queries             ->  mean of both fold models (neither saw V-in skeletons)

Query spectra are rebuilt with the same deterministic code that built the tables.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE  # noqa: E402
from casmi.channels.fpinfer import molecule_logp  # noqa: E402
from casmi.channels.fpmodel import candidate_scores, load  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

spec = importlib.util.spec_from_file_location("brd", ROOT / "scripts" / "build_rank_data.py")
brd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(brd)


def add_context(tbl, cols):
    g = tbl.groupby(["qid", "cls"], sort=False)
    for c in cols:
        tbl[f"{c}_rank"] = g[c].rank(ascending=False, method="min")
        tbl[f"{c}_gap"] = g[c].transform("max") - tbl[c]
        sd = g[c].transform("std").replace(0, np.nan)
        tbl[f"{c}_z"] = ((tbl[c] - g[c].transform("mean")) / sd).fillna(0.0)
    return tbl


def add_d3(tbl, items, models, pool, fps, device):
    """D3 columns for every (query, class) group of `tbl`; `items` = [(Query, fold)]."""
    by_qid = {q.qid: (q, fold) for q, fold in items}
    if "gen_smiles" in tbl:
        g = tbl[tbl.gen == 1]
        pool.register_generated(g.pool_key.to_numpy(), g.gen_smiles.to_numpy())
    ll, cos = np.zeros(len(tbl)), np.zeros(len(tbl))
    for qid, idx in tbl.groupby("qid", sort=False).indices.items():
        q, fold = by_qid[qid]
        use = [models[0], models[1]] if fold == -1 else [models[1 - fold]]
        lp_on, lp_off = molecule_logp(use, q.spectra, device)
        bits = np.unpackbits(pool.packed_fps(tbl.pool_key.to_numpy()[idx], fps), axis=1).astype(np.float32)
        ll[idx], cos[idx] = candidate_scores(lp_on, lp_off, bits)
    tbl = tbl.assign(d3_ll=ll, d3_cos=cos)
    return add_context(tbl, ["d3_ll", "d3_cos"])


def load_everything(device):
    models = {f: load(ARTIFACT_DIR / "models" / f"fp_fold{f}.pt", device) for f in (0, 1)}
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE)
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    train_items = brd.training_queries(store, key_of_row, pool, vin_skel, 1500)
    return models, pool, fps, [(q, -1) for q in vin], train_items


def main():
    """v1 tables: rank_{vin,train}_[0-9].parquet -> rank_{vin,train}_d3.parquet"""
    device = "cuda" if torch.cuda.is_available() else "cpu"
    models, pool, fps, vin_items, train_items = load_everything(device)
    for name, items in (("vin", vin_items), ("train", train_items)):
        files = sorted(ARTIFACT_DIR.glob(f"rank_{name}_[0-9].parquet"))
        tbl = add_d3(pd.concat([pd.read_parquet(f) for f in files], ignore_index=True), items, models, pool, fps, device)
        tbl.to_parquet(ARTIFACT_DIR / f"rank_{name}_d3.parquet", index=False)
        print(f"{name}: {tbl.qid.nunique()} queries, {len(tbl):,} rows -> rank_{name}_d3.parquet", flush=True)


if __name__ == "__main__":
    main()
