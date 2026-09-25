"""Join D4 (fragment explanation) features onto ranker / validation tables.

Like D3, a D4 score depends only on the query spectra and the candidate structure, not on the
knockout, so it is added after the tables are built. Query spectra are rebuilt with the same
deterministic builders that made the tables. Fragment masses are cached across tables
(artifacts/frag_cache.pkl): the same candidates recur in every pool variant.

  python scripts/add_d4_features.py --table rank_train_v13_pcl_d3 --queries train,e180 --pool _pcl
  python scripts/add_d4_features.py --table rank_vreal_pcl_d3 --queries vreal --pool _pcl
  python scripts/add_d4_features.py --table rank_vin_pcl_d3 --queries vin --pool _pcl
Writes <table>_d4.parquet next to the input.
"""
import argparse
import importlib.util
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.channels.fragments import FragmentScorer, add_d4  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

spec = importlib.util.spec_from_file_location("brd", ROOT / "scripts" / "build_rank_data.py")
brd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(brd)
CACHE = ARTIFACT_DIR / "frag_cache.pkl"


def queries(kinds, store, key_of_row, pool):
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    out = []
    for kind in kinds:
        if kind == "vin":
            out += [(q, -1) for q in vin]
        elif kind == "vreal":
            out += brd.vreal_queries(store, key_of_row, vin_skel)[0]
        elif kind == "train":
            out += brd.training_queries(store, key_of_row, pool, vin_skel, 1500)
        elif kind == "e180":
            out += brd.e180_class1_queries(store, key_of_row, pool, vin_skel, 1200)
    return {q.qid: q for q, _ in out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--pool", default="_pcl", help="pool tag the table was built on (for candidate SMILES)")
    args = ap.parse_args()
    t0 = time.time()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / f"pool{args.pool}.parquet"))
    by_qid = queries(args.queries.split(","), store, key_of_row, pool)
    tbl = pd.read_parquet(ARTIFACT_DIR / f"{args.table}.parquet")
    if "gen_smiles" in tbl:
        g = tbl[tbl.gen == 1]
        pool.register_generated(g.pool_key.to_numpy(), g.gen_smiles.to_numpy())
    scorer = FragmentScorer(pool.smiles_of)
    if CACHE.exists():
        scorer.cache.update(pickle.loads(CACHE.read_bytes()))
    n0 = len(scorer.cache)
    print(f"{args.table}: {len(tbl):,} rows, {tbl.qid.nunique()} queries, {n0:,} cached structures "
          f"({time.time() - t0:.0f}s)", flush=True)
    parts = []
    groups = tbl.groupby(["qid", "cls"], sort=False)
    for i, ((qid, _), g) in enumerate(groups):
        parts.append(add_d4(g.reset_index(drop=True), scorer, by_qid[qid].spectra))
        if (i + 1) % 500 == 0:
            print(f"  {i + 1}/{groups.ngroups} groups, {len(scorer.cache) - n0:,} new structures "
                  f"({time.time() - t0:.0f}s)", flush=True)
    out = pd.concat(parts, ignore_index=True)
    out.to_parquet(ARTIFACT_DIR / f"{args.table}_d4.parquet", index=False)
    CACHE.write_bytes(pickle.dumps(scorer.cache, protocol=pickle.HIGHEST_PROTOCOL))
    print(f"wrote {args.table}_d4.parquet ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
