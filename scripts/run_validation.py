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
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE, POOL_TAG, PROJECT_ROOT  # noqa: E402
from casmi.channels.library import LibraryIndex  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.engine import BaselineEngine, EngineConfig, molecule_mass  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.pool.fingerprint_provenance import verify_fingerprint_provenance  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import (  # noqa: E402
    CLASSES, build_queries, knockout_for, mixture_grid, include_equivalent_source_rows,
)
from casmi.validate.metric import rank_of_truth, summarize  # noqa: E402
from casmi.validate.run_alignment import query_signature, SIGNATURE_VERSION  # noqa: E402
from casmi.validate.spectral_controls import perturb_spectra, MODES  # noqa: E402
from casmi.validate.error_slices import describe_query  # noqa: E402
from casmi.validate.source_manifest import (  # noqa: E402
    MANIFEST_SCHEMA, verify_manifest, require_validation_coverage,
)


def equivalent_keys(query, pool, keys: MetricKeys, cfg):
    """Pool keys in the query's mass window that are the same structure under the metric.
    The cheap skeleton key filters out everything that cannot match before canonicalising."""
    truth_smiles = pool.smiles_of(query.pool_key)
    truth, skel = keys.get(query.pool_key, truth_smiles), skeleton_key(truth_smiles)
    win = pool.window(molecule_mass(query.spectra, cfg.mass_offset_ppm), cfg.tol_ppm * 3)
    same = {query.pool_key}
    for k in pool.key[win.idx]:
        if k == query.pool_key or skeleton_key(pool.smiles_of(k)) != skel:
            continue
        if keys.get(k, pool.smiles_of(k)) == truth:
            same.add(k)
    return truth, same


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-queries", type=int, default=None)
    ap.add_argument("--library", default="enveda-np-examples", help="query library")
    ap.add_argument("--max-spectra", type=int, default=None, help="spectra sampled per molecule")
    ap.add_argument("--query-seed", type=int, default=0, help="query/spectrum sampling seed")
    ap.add_argument("--tag", default="d1-baseline")
    ap.add_argument("--offset-ppm", type=float, default=0.0)
    ap.add_argument("--analog", action="store_true", help="enable D2 analog propagation")
    ap.add_argument("--require-fingerprint-provenance", action="store_true",
                    help="fail closed if the D2 pool/fingerprint source snapshot is absent")
    ap.add_argument("--fresh-analog-index", action="store_true",
                    help="rebuild D2 reference index from spectra, bypass all index caches")
    ap.add_argument("--analog-selection-policy", choices=("legacy", "source_balanced"),
                    default="legacy", help="D2 index reference selection; legacy keeps historical baseline")
    ap.add_argument("--spectral-control", choices=MODES, default="none",
                    help="diagnostic only: suppress or perturb fragments, not a model setting")
    ap.add_argument("--source-manifest", type=Path, default=None,
                    help="optional SHA-256 artifact inventory; fail closed on mismatch")
    ap.add_argument("--control-seed", type=int, default=0,
                    help="per-molecule deterministic negative-control seed")
    args = ap.parse_args()

    if any(n is not None and n < 1 for n in (args.max_queries, args.max_spectra)):
        ap.error("--max-queries and --max-spectra must be positive")
    if args.control_seed < 0:
        ap.error("--control-seed must be nonnegative")
    if args.require_fingerprint_provenance and not args.analog:
        ap.error("--require-fingerprint-provenance requires --analog")
    if args.fresh_analog_index and not args.analog:
        ap.error("--fresh-analog-index requires --analog")
    if args.analog_selection_policy != "legacy" and not args.analog:
        ap.error("--analog-selection-policy requires --analog")
    cfg = EngineConfig()
    cfg = replace(cfg, mass_offset_ppm=args.offset_ppm,
                  analog=replace(cfg.analog, reference_policy=args.analog_selection_policy))
    t0 = time.time()
    source_status = None
    source_required_inputs = None
    if args.source_manifest is not None:
        with args.source_manifest.open("r", encoding="utf-8") as stream:
            recorded_manifest = json.load(stream)
        source_required_inputs = require_validation_coverage(
            recorded_manifest, pool_tag=POOL_TAG, analog=args.analog,
        )
        source_status = verify_manifest(ARTIFACT_DIR, recorded_manifest)
        print(f"Verified {source_status['files_checked']} local source files "
              f"(including {len(source_required_inputs)} required inputs)", flush=True)
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    library = LibraryIndex(store, key_of_row, cfg.library)
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    queries = build_queries(store, key_of_row, library=args.library,
                            max_queries=args.max_queries, max_spectra=args.max_spectra,
                            seed=args.query_seed)
    if not queries:
        raise ValueError('No eligible query molecules for validation')
    # Persist identities and input fingerprints after C1 exclusion metadata is
    # finalized. Otherwise runs with different knockout rules might be paired.
    query_ids = [str(q.qid) for q in queries]
    query_signatures = []
    source_query_signatures = []
    query_context = []
    # Some tautomers share the official metric key while having distinct
    # standard pool keys. C1 must also hide these structures in the query's
    # own library. Do not use the filtered D1 index to infer all raw rows.
    source_idx = np.flatnonzero((store.meta.ingest_lib == args.library).to_numpy())
    analog = fps = None
    fingerprint_status = None
    analog_index_status = None
    if args.analog:
        from casmi.channels.analog import AnalogIndex
        fp_keys = pd.read_parquet(ARTIFACT_DIR / POOL_FP_KEYS_FILE).pool_key.to_numpy()
        if not np.array_equal(fp_keys, pool.key):
            raise ValueError("D2 fingerprint keys are not aligned with the candidate pool")
        fps = np.load(ARTIFACT_DIR / POOL_FP_FILE, allow_pickle=False, mmap_mode="r")
        provenance_file = ARTIFACT_DIR / f"pool_fp_provenance{POOL_TAG}.json"
        if provenance_file.exists():
            with provenance_file.open("r", encoding="utf-8") as stream:
                recorded_fingerprint = json.load(stream)
            fingerprint_status = verify_fingerprint_provenance(
                pool.table, fp_keys, fps, recorded_fingerprint,
            )
            print("Verified D2 candidate SMILES/mass/fingerprint provenance",
                  flush=True)
        elif args.require_fingerprint_provenance:
            raise ValueError(
                "No D2 fingerprint provenance manifest; rebuild a separate "
                "CASMI_POOL_TAG variant instead of reusing frozen arrays"
            )
        else:
            print("Warning: legacy D2 fingerprint source linkage not verified; "
                  "new experiments should use --require-fingerprint-provenance",
                  flush=True)
        index_cache = None if args.fresh_analog_index else ARTIFACT_DIR / "analog_index.npz"
        analog = AnalogIndex(store, key_of_row, cfg.analog, cache=index_cache)
        if args.fresh_analog_index:
            analog_index_status = analog.fresh_index_provenance()
            print("D2 analog reference matrix rebuilt from source (cache bypassed)",
                  flush=True)
        else:
            print("Warning: D2 index cache is not decision-grade evidence; "
                  "use --fresh-analog-index with a signed input manifest",
                  flush=True)
        print(f"analog index: {analog.matrix.shape[0]:,} spectra, {analog.matrix.nnz:,} nonzeros", flush=True)
    engine = BaselineEngine(pool, library, cfg, analog, fps)
    print(f"loaded in {time.time() - t0:.0f}s | pool {len(pool):,} | library structures {len(library.rows_of):,} "
          f"| queries {len(queries)}", flush=True)

    ranks = {c: [] for c in CLASSES}
    recall = {c: [] for c in CLASSES}
    n_cand, n_equiv = [], []
    for qi, q in enumerate(queries):
        query_context.append(describe_query(q.spectra))
        truth, same = equivalent_keys(q, pool, keys, cfg)
        include_equivalent_source_rows(q, same, source_idx, key_of_row)
        source_query_signatures.append(query_signature(q))
        if args.spectral_control != "none":
            q = replace(q, spectra=perturb_spectra(
                q.spectra, mode=args.spectral_control, seed=args.control_seed,
                query_id=str(q.qid),
            ))
        query_signatures.append(query_signature(q))
        n_equiv.append(len(same))
        per_class = {}
        for cls in ("C1", "C2"):
            ko = knockout_for(q, cls, library.rows_of, same)
            per_class[cls] = engine.features([dict(s) for s in q.spectra], ko.excluded_rows, ko.excluded_keys)
        # C3 = C2 with the structure (and its tautomer variants) removed from the pool
        per_class["C3"] = per_class["C2"][~per_class["C2"].pool_key.isin(same)].reset_index(drop=True)
        for cls, feats in per_class.items():
            ranked = keys.dedupe_ranked(engine.rank(feats), pool.smiles_of, cfg.top_n)
            ranks[cls].append(rank_of_truth([mk for _, mk in ranked], truth))
            recall[cls].append(bool(set(feats.pool_key) & same))
        n_cand.append(len(per_class["C2"]))
        if (qi + 1) % 25 == 0:
            keys.save()
            print(f"  {qi + 1}/{len(queries)} queries, {time.time() - t0:.0f}s", flush=True)
    keys.save()

    config = {**json.loads(json.dumps(asdict(cfg), default=str)), "channels": ["d1"] + (["d2"] if args.analog else [])}
    report = {"tag": args.tag, "config": config,
              "fingerprint_provenance": fingerprint_status,
              "analog_index_provenance": analog_index_status,
              "source_manifest_sha256": (
                  source_status["files_digest_sha256"] if source_status is not None else None
              ),
              "source_manifest_file_count": (
                  source_status["files_checked"] if source_status is not None else None
              ),
              "source_manifest_schema": MANIFEST_SCHEMA if source_status is not None else None,
              "source_manifest_pool_tag": POOL_TAG if source_status is not None else None,
              "source_manifest_coverage": (
                  "required-validation-inputs-v1" if source_status is not None else None
              ),
              "source_manifest_required_inputs": (
                  list(source_required_inputs) if source_required_inputs is not None else None
              ),
              "spectral_control": {
                  "mode": args.spectral_control, "seed": args.control_seed,
                  "interpretation": "negative-control diagnostic, not submission/model improvement",
              },
              "n_queries": len(queries), "query_ids": query_ids,
              "query_signatures": query_signatures,
              "query_context": query_context,
              "query_signature_version": SIGNATURE_VERSION,
              "query_sampling": {"library": args.library,
                                 "max_queries": args.max_queries,
                                 "max_spectra": args.max_spectra,
                                 "seed": args.query_seed,
                                 "spectral_control": args.spectral_control,
                                 "control_seed": args.control_seed},
              "median_candidates": float(np.median(n_cand)),
              "queries_with_tautomer_duplicates": int(np.sum(np.array(n_equiv) > 1)), "classes": {}}
    print(f"\nmedian candidates per molecule: {np.median(n_cand):.0f} | "
          f"queries whose structure also sits in the pool under another key: {report['queries_with_tautomer_duplicates']}")
    for cls in CLASSES:
        s = summarize(ranks[cls])
        report["classes"][cls] = {**s.__dict__, "pool_recall": float(np.mean(recall[cls])), "ranks": ranks[cls]}
        print(f"{cls}: {s}  pool recall {np.mean(recall[cls]):.3f}")
    grid = mixture_grid({c: report["classes"][c]["mrr"] for c in CLASSES})
    report["mixture"] = grid.to_dict("records")
    if args.spectral_control != "none":
        report["spectral_control"]["source_query_signatures"] = source_query_signatures
        print("NEGATIVE CONTROL: transformed fragments; diagnostic scores are NOT predicted leaderboard performance.")
    print("\n" + ("DIAGNOSTIC control mixture" if args.spectral_control != "none"
                   else "expected LB over plausible class mixes")
          + " (rows f1, columns f3):\n"
          + grid.pivot(index="f1", columns="f3", values="score").round(4).to_string())

    manifest = {"config": report["config"], "query_sampling": report["query_sampling"],
                "query_ids": query_ids, "query_signatures": query_signatures}
    h = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:12]
    out = PROJECT_ROOT / "runs" / f"{time.strftime('%Y%m%d-%H%M%S')}_{args.tag}_{h}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x", encoding="utf-8") as file:
        json.dump(report, file, indent=2)
    print(f"\nlogged runs/{out.name} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
