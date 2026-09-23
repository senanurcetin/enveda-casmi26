"""Knockout validation: score an engine per simulated class, with CIs, and log the run.

Usage: python scripts/run_validation.py [--max-queries N] [--tag name] [--offset-ppm X]
"""
import argparse
import hashlib
import json
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, PROJECT_ROOT  # noqa: E402
from casmi.channels.library import LibraryIndex  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.engine import BaselineEngine, EngineConfig, molecule_mass  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import CLASSES, build_queries, knockout_for, mixture_grid  # noqa: E402
from casmi.validate.metric import rank_of_truth, summarize  # noqa: E402


def equivalent_keys(query, pool, keys: MetricKeys, cfg):
    """Pool keys in the query's mass window that are the same structure under the metric."""
    truth = keys.get(query.pool_key, pool.smiles_of(query.pool_key))
    win = pool.window(molecule_mass(query.spectra, cfg.mass_offset_ppm), cfg.tol_ppm * 3)
    same = {query.pool_key}
    for k in pool.key[win.idx]:
        if k != query.pool_key and keys.get(k, pool.smiles_of(k)) == truth:
            same.add(k)
    return truth, same


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-queries", type=int, default=None)
    ap.add_argument("--tag", default="d1-baseline")
    ap.add_argument("--offset-ppm", type=float, default=0.0)
    args = ap.parse_args()

    cfg = replace(EngineConfig(), mass_offset_ppm=args.offset_ppm)
    t0 = time.time()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool.parquet"))
    library = LibraryIndex(store, key_of_row, cfg.library)
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    queries = build_queries(store, key_of_row, max_queries=args.max_queries)
    engine = BaselineEngine(pool, library, cfg)
    print(f"loaded in {time.time() - t0:.0f}s | pool {len(pool):,} | library structures {len(library.rows_of):,} "
          f"| queries {len(queries)}", flush=True)

    ranks = {c: [] for c in CLASSES}
    recall = {c: [] for c in CLASSES}
    n_cand, n_equiv = [], []
    for qi, q in enumerate(queries):
        truth, same = equivalent_keys(q, pool, keys, cfg)
        n_equiv.append(len(same))
        for cls in CLASSES:
            ko = knockout_for(q, cls, library.rows_of, same)
            feats = engine.features([dict(s) for s in q.spectra], ko.excluded_rows, ko.excluded_keys)
            ranked = keys.dedupe_ranked(engine.rank(feats), pool.smiles_of, cfg.top_n)
            ranks[cls].append(rank_of_truth([mk for _, mk in ranked], truth))
            recall[cls].append(bool(set(feats.pool_key) & same))
            if cls == "C2":
                n_cand.append(len(feats))
        if (qi + 1) % 25 == 0:
            keys.save()
            print(f"  {qi + 1}/{len(queries)} queries, {time.time() - t0:.0f}s", flush=True)
    keys.save()

    report = {"tag": args.tag, "config": json.loads(json.dumps(asdict(cfg), default=str)),
              "n_queries": len(queries), "median_candidates": float(np.median(n_cand)),
              "queries_with_tautomer_duplicates": int(np.sum(np.array(n_equiv) > 1)), "classes": {}}
    print(f"\nmedian candidates per molecule: {np.median(n_cand):.0f} | "
          f"queries whose structure also sits in the pool under another key: {report['queries_with_tautomer_duplicates']}")
    for cls in CLASSES:
        s = summarize(ranks[cls])
        report["classes"][cls] = {**s.__dict__, "pool_recall": float(np.mean(recall[cls])), "ranks": ranks[cls]}
        print(f"{cls}: {s}  pool recall {np.mean(recall[cls]):.3f}")
    grid = mixture_grid({c: report["classes"][c]["mrr"] for c in CLASSES})
    report["mixture"] = grid.to_dict("records")
    print("\nexpected LB over plausible class mixes (rows f1, columns f3):\n"
          + grid.pivot(index="f1", columns="f3", values="score").round(4).to_string())

    h = hashlib.sha1(json.dumps(report["config"], sort_keys=True).encode()).hexdigest()[:8]
    out = PROJECT_ROOT / "runs" / f"{time.strftime('%Y%m%d-%H%M')}_{args.tag}_{h}.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"\nlogged runs/{out.name} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
