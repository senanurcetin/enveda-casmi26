"""Aggregate-only retrospective study of frozen models under corrected knockouts.

All source artifacts must be rebuilt in a fresh scratch directory. This runner
does not train models, create predictions for test.parquet, or submit to Kaggle.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import sys
import time
from dataclasses import asdict, replace
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def main():
    import numpy as np
    import pandas as pd
    import rdkit
    import torch

    from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE, POOL_TAG
    from casmi.channels.analog import AnalogIndex
    from casmi.channels.fpmodel import load
    from casmi.channels.library import LibraryIndex
    from casmi.chem.keys import MetricKeys
    from casmi.engine import BaselineEngine, EngineConfig
    from casmi.engine import molecule_mass
    from casmi.pool.pool import CandidatePool
    from casmi.pool.fingerprint_provenance import verify_fingerprint_provenance
    from casmi.predict import Predictor
    from casmi.spectra.store import SpectrumStore
    from casmi.validate.knockout import build_queries, include_equivalent_source_rows, knockout_for, CLASSES
    from casmi.validate.knockout import equivalent_keys as pool_equivalents
    from casmi.validate.metric import rank_of_truth, summarize
    from casmi.validate.run_alignment import query_signature, SIGNATURE_VERSION
    from casmi.validate.source_manifest import create_manifest, required_validation_inputs, verify_manifest
    from casmi.validate.spectral_controls import perturb_spectra
    from casmi.validate.paired_bootstrap import paired_difference
    from casmi.validate.simultaneous_ci import compare_mixtures
    from scripts.make_submission import LoadedRanker
    from scripts.run_validation import equivalent_keys
    from scripts.compare_control_runs import diagnostic_summary

    if rdkit.__version__ != "2026.03.3":
        raise RuntimeError("The official metric RDKit version is required")
    started = time.time()
    np.random.seed(0)
    torch.manual_seed(0)
    torch.set_num_threads(2)
    torch.use_deterministic_algorithms(True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    cfg = EngineConfig()
    full_ranker = LoadedRanker("ranker-v16-gen2")
    if len(full_ranker.models) != 5:
        raise ValueError("Exactly five frozen ranker seed models are required")
    columns = set(full_ranker.columns)
    full_cfg = replace(cfg, generate="gen" in columns, gen_consensus="gen_n_seeds" in columns,
                       np_score="np_score" in columns)
    static = list(required_validation_inputs(pool_tag=POOL_TAG, analog=True))
    static += ["structures.parquet", f"pool_source_quality{POOL_TAG}.json"]
    static += ["models/fp_fold0.pt", "models/fp_fold1.pt"]
    static += [f"models/ranker-v16-gen2_seed{seed}.txt" for seed in range(5)]
    manifest = create_manifest(ARTIFACT_DIR, static)
    verify_manifest(ARTIFACT_DIR, manifest)
    store = SpectrumStore(ARTIFACT_DIR / "store")
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    truth_pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "structures.parquet"))
    source_quality = json.loads((ARTIFACT_DIR / f"pool_source_quality{POOL_TAG}.json").read_text())
    fp_keys = pd.read_parquet(ARTIFACT_DIR / POOL_FP_KEYS_FILE).pool_key.to_numpy()
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE, mmap_mode="r", allow_pickle=False)
    provenance = json.loads((ARTIFACT_DIR / f"pool_fp_provenance{POOL_TAG}.json").read_text())
    fingerprint_status = verify_fingerprint_provenance(pool.table, fp_keys, fps, provenance)
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    library = LibraryIndex(store, key_of_row, cfg.library)
    analog = AnalogIndex(store, key_of_row, cfg.analog, cache=None)
    analog_status = analog.fresh_index_provenance()
    fp_models = [load(ARTIFACT_DIR / "models" / f"fp_fold{fold}.pt", device) for fold in (0, 1)]
    predictors = {
        "d1": Predictor(BaselineEngine(pool, library, cfg), pool, None),
        "d1_d2": Predictor(BaselineEngine(pool, library, cfg, analog, fps), pool, None),
        "frozen_v16": Predictor(BaselineEngine(pool, library, full_cfg, analog, fps), pool, None,
                                ranker=full_ranker, fp_models=fp_models, device=device),
    }
    keys = MetricKeys(ARTIFACT_DIR / "study_metric_keys.parquet")
    for predictor in predictors.values():
        predictor.keys = keys
    queries = build_queries(store, key_of_row, library="enveda-np-examples",
                            max_spectra=5, max_queries=None, seed=0)
    if not queries:
        raise ValueError("No eligible designated V-in queries")
    source_rows = np.flatnonzero((store.meta.ingest_lib == "enveda-np-examples").to_numpy())
    modes = ("none", "blank", "intensity_shuffle")
    query_ids = [str(query.qid) for query in queries]
    source_signatures = []
    transformed_signatures = {mode: [] for mode in modes}
    ranks = {model: {mode: {cls: [] for cls in CLASSES} for mode in modes} for model in predictors}
    extra_exclusions = []
    equivalent_counts = []
    print(f"Prepared {len(queries)} eligible V-in molecules; full store has {len(store):,} spectra", flush=True)
    for index, query in enumerate(queries):
        # Labels and reference exclusions survive candidate-source quarantine.
        # Never remove a query because its truth is absent from the candidate pool.
        truth, same = equivalent_keys(query, truth_pool, keys, cfg)
        candidate_truth, candidate_same = pool_equivalents(
            query.pool_key, molecule_mass(query.spectra, cfg.mass_offset_ppm), pool, keys,
            window_ppm=cfg.tol_ppm * 3, truth_smiles=truth_pool.smiles_of(query.pool_key),
        )
        if candidate_truth != truth:
            raise ValueError("Independent label canonicalization mismatch")
        same.update(candidate_same)
        if not truth:
            raise ValueError("A query truth could not be canonicalized")
        include_equivalent_source_rows(query, same, source_rows, key_of_row)
        source_signatures.append(query_signature(query))
        extra_exclusions.append(int(len(np.setdiff1d(query.source_library_rows, query.spectrum_rows))))
        equivalent_counts.append(len(same))
        for mode in modes:
            spectra = perturb_spectra(query.spectra, mode=mode, seed=0, query_id=str(query.qid))
            transformed_signatures[mode].append(query_signature(replace(query, spectra=spectra)))
            for cls in CLASSES:
                ko = knockout_for(query, cls, library.rows_of, same)
                for model, predictor in predictors.items():
                    # C3 generation must run AFTER pool exclusion; filtering a C2
                    # list afterward would miss truth-regeneration opportunities.
                    ranked = predictor.rank(spectra, excluded_rows=ko.excluded_rows,
                                            excluded_keys=ko.excluded_keys,
                                            knocked_keys=same if cls != "C1" else frozenset())
                    ranks[model][mode][cls].append(int(rank_of_truth([key for _, key in ranked], truth)))
        if (index + 1) % 10 == 0:
            print(f"Completed {index + 1}/{len(queries)} paired molecules; {time.time() - started:.0f}s", flush=True)
            keys.save()
    keys.save()
    # No source inputs may change while caches/temporary generated candidates evolve.
    verify_manifest(ARTIFACT_DIR, manifest)
    models = {}
    private_reports = {}
    for model in predictors:
        config = {"engine": asdict(full_cfg if model == "frozen_v16" else cfg), "model": model}
        private_reports[model] = {}
        models[model] = {"conditions": {}, "controls": {}}
        for mode in modes:
            blocks = {cls: {**asdict(summarize(ranks[model][mode][cls])), "ranks": ranks[model][mode][cls]} for cls in CLASSES}
            private_reports[model][mode] = {
                "n_queries": len(queries), "query_ids": query_ids,
                "query_signatures": transformed_signatures[mode], "query_signature_version": SIGNATURE_VERSION,
                "source_manifest_sha256": manifest["files_digest_sha256"], "config": config,
                "spectral_control": {"mode": mode, "seed": 0, "source_query_signatures": source_signatures},
                "query_sampling": {"library": "enveda-np-examples", "max_spectra": 5,
                                   "max_queries": None, "seed": 0, "spectral_control": mode, "control_seed": 0},
                "classes": blocks,
            }
            models[model]["conditions"][mode] = {cls: {key: value for key, value in blocks[cls].items() if key != "ranks"} for cls in CLASSES}
        for mode in modes[1:]:
            models[model]["controls"][mode] = diagnostic_summary(private_reports[model]["none"], private_reports[model][mode], n_boot=4000, seed=0)
    comparisons = {}
    if transformed_signatures["none"] != source_signatures:
        raise ValueError("Unperturbed query signatures changed")
    for left, right in (("d1", "d1_d2"), ("d1_d2", "frozen_v16")):
        class_effects = {}
        for cls in CLASSES:
            delta, low, high = paired_difference(ranks[left]["none"][cls], ranks[right]["none"][cls], n_boot=4000, seed=0)
            class_effects[cls] = {"delta_mrr": delta, "ci_low": low, "ci_high": high}
        mixture = compare_mixtures(ranks[left]["none"], ranks[right]["none"], n_boot=4000, seed=0)
        mixture.pop("verdict", None)
        comparisons[f"{right}_minus_{left}"] = {"status": "RETROSPECTIVE_PAIRED_DIAGNOSTIC", "classes": class_effects, "mixture": mixture}
    result = {
        "schema_version": 1, "status": "COMPLETE_RETROSPECTIVE_REVALIDATION",
        "checked_at": datetime.now(UTC).isoformat(), "elapsed_seconds": round(time.time() - started, 1),
        "cohort": {"library": "enveda-np-examples", "selection": "all eligible molecules; no max-query subsample",
                   "n_molecules": len(queries), "max_spectra": 5, "query_seed": 0,
                   "identity_sha256": digest_json(query_ids), "source_signatures_sha256": digest_json(source_signatures)},
        "source_artifact_manifest": manifest, "fingerprint_provenance": fingerprint_status,
        "candidate_source_quality": source_quality,
        "queries_with_truth_absent_from_candidate_pool": sum(query.pool_key not in pool.row_of for query in queries),
        "analog_index_provenance": analog_status,
        "knockout_audit": {"queries_with_additional_same_library_exclusions": sum(value > 0 for value in extra_exclusions),
                            "additional_same_library_rows_excluded": sum(extra_exclusions),
                            "queries_with_metric_equivalent_pool_variants": sum(value > 1 for value in equivalent_counts)},
        "model_configurations": {"d1": asdict(cfg), "d1_d2": asdict(cfg), "frozen_v16": asdict(full_cfg)},
        "models": models, "paired_comparisons": comparisons,
        "private_detail_sha256": digest_json(private_reports),
        "environment": {"python": sys.version.split()[0], "device": device, **{name: importlib.metadata.version(name) for name in ("rdkit", "numpy", "pandas", "pyarrow", "scipy", "numba", "torch", "lightgbm")}},
        "limitations": [
            "Retrospective designated V-in evaluation, not a new untouched external holdout or competition score.",
            "Frozen learned weights were not retrained; their complete training-data lineage is not independently attested.",
            "Cohort was used in historical research; intervals do not remove adaptive model-selection or domain-shift bias.",
            "Controls diagnose fragment reliance; intensity shuffling retains peak positions and precursor information.",
            "Mixture intervals cover only the predeclared grid; no ACCEPT/REJECT model-promotion claim.",
            "No raw spectra, molecular structures, query identifiers, ranks or model files are exported.",
        ],
    }
    out = Path(os.environ.get("CASMI_STUDY_OUTPUT", "/kaggle/working/revalidation_study.json"))
    with out.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Completed retrospective study of {len(queries)} molecules; aggregate-only report written", flush=True)


if __name__ == "__main__":
    main()
