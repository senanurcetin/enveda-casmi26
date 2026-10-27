"""Execute fixed C3 budgets with fresh weights and independently audited external data.

Fail closed on insufficient external novelty. Test features are computed only after
both policies, all learned artifacts and the development choice have been sealed.
Only aggregate outcomes and private-detail hashes may leave scratch storage.
"""
from __future__ import annotations

from dataclasses import asdict
from datetime import UTC, datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
POLICIES = ("legacy", "expanded")


def primary_coverage_interval(legacy, expanded, *, seed, n_boot=4000):
    import numpy as np
    from casmi.validate.rank_inputs import validated_bootstrap
    validated_bootstrap(n_boot, seed)
    a, b = np.asarray(legacy), np.asarray(expanded)
    if a.ndim != 1 or a.shape != b.shape or len(a) < 2 or not np.isin(a, [0, 1]).all() or not np.isin(b, [0, 1]).all():
        raise ValueError("Paired binary molecule coverage required")
    delta = b.astype(float) - a.astype(float)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for start in range(0, n_boot, 100):
        ids = rng.integers(len(delta), size=(min(100, n_boot - start), len(delta)))
        draws[start:start + len(ids)] = delta[ids].mean(axis=1)
    low, high = np.quantile(draws, [.025, .975])
    return {"n_molecules": len(a), "legacy_covered": int(a.sum()), "expanded_covered": int(b.sum()),
            "gained": int(((a == 0) & (b == 1)).sum()), "lost": int(((a == 1) & (b == 0)).sum()),
            "delta_coverage": float(delta.mean()), "ci_low": float(low), "ci_high": float(high),
            "n_bootstrap": n_boot, "seed": seed, "method": "Molecule-paired binary coverage percentile bootstrap"}


def development_choice(legacy, expanded):
    import numpy as np
    a, b = np.asarray(legacy), np.asarray(expanded)
    if a.ndim != 1 or a.shape != b.shape or not len(a) or not np.isin(a, [0, 1]).all() or not np.isin(b, [0, 1]).all():
        raise ValueError("Paired development coverage required")
    gained, lost = int(((a == 0) & (b == 1)).sum()), int(((a == 1) & (b == 0)).sum())
    return {"selected_policy": "expanded" if lost == 0 and gained >= 1 else "legacy",
            "n_molecules": len(a), "legacy_covered": int(a.sum()), "expanded_covered": int(b.sum()),
            "gained": gained, "lost": lost, "policy": "Prespecified development-only C3 coverage gate; both test policies reported"}


def main():
    import numpy as np
    import pandas as pd
    import rdkit
    import torch
    from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE, POOL_TAG
    from casmi.channels.library import LibraryIndex
    from casmi.channels.analog import AnalogIndex
    from casmi.chem.keys import MetricKeys
    from casmi.engine import BaselineEngine, EngineConfig, molecule_mass
    from casmi.pool.pool import CandidatePool
    from casmi.pool.fingerprint_provenance import verify_fingerprint_provenance
    from casmi.spectra.store import SpectrumStore
    from casmi.predict import Predictor
    from casmi.rank.features import feature_columns
    from casmi.rank.ranker import Ranker, RankerConfig
    from casmi.validate.budget_features import budget_features
    from casmi.validate.external_population import prepare_internal, select_external
    from casmi.validate.cohort_bootstrap import joint_cohort_intervals
    from casmi.validate.knockout import Query, CLASSES, knockout_for, equivalent_keys, include_equivalent_source_rows
    from casmi.validate.spectral_controls import perturb_spectra
    from casmi.validate.metric import rank_of_truth, summarize
    from casmi.validate.scaffold_study import digest
    from casmi.validate.source_manifest import create_manifest, required_validation_inputs, verify_manifest, sha256_file
    from scripts.run_external_audit import audit
    from scripts.run_advanced_study import SearchMemo
    from scripts.external_study_training import fit_fingerprints

    started = time.time()
    protocol_path = ROOT / "protocols/external-c3-v1.json"
    execution_path = ROOT / "protocols/external-c3-v1-execution.json"
    source_protocol_path = ROOT / os.environ.get("CASMI_EXTERNAL_PROTOCOL", "protocols/external-source-v2.json")
    protocol = json.loads(protocol_path.read_text())
    execution = json.loads(execution_path.read_text())
    seed, device = protocol["seed"], "cpu"
    if rdkit.__version__ != "2026.03.3":
        raise ValueError("Official metric RDKit required")
    np.random.seed(seed); torch.manual_seed(seed); torch.set_num_threads(2)
    torch.use_deterministic_algorithms(True)
    independence, novel = audit(Path(os.environ["CASMI_DATA"]) / "train.parquet",
        Path(os.environ["CASMI_EXTERNAL_INPUT"]), Path(os.environ["CASMI_EXTERNAL_MANIFEST"]), source_protocol_path)
    if independence["status"] != "EXTERNAL_INDEPENDENCE_FEASIBLE":
        raise ValueError("External independence gate failed; no supervised fitting or predictions permitted")
    store = SpectrumStore(ARTIFACT_DIR / "store")
    meta = store.meta
    row_keys = meta.train_ik14.astype(str).to_numpy()
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    truth_pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "structures.parquet"))
    previous_protocol = json.loads((ROOT / "protocols/advanced-internal-v1.json").read_text())
    previous_audit = json.loads((ROOT / "reports/advanced_study.json").read_text())["partitions"]
    splits, groups, labels, split_audit = prepare_internal(meta, truth_pool, previous_protocol, previous_audit, protocol)
    external_items, external_population = select_external(novel, splits, protocol["split_salt"])
    static_names = list(required_validation_inputs(pool_tag=POOL_TAG, analog=True)) + ["structures.parquet", f"pool_source_quality{POOL_TAG}.json"]
    manifest = create_manifest(ARTIFACT_DIR, static_names)
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE, mmap_mode="r", allow_pickle=False)
    fp_keys = pd.read_parquet(ARTIFACT_DIR / POOL_FP_KEYS_FILE).pool_key.to_numpy()
    fp_status = verify_fingerprint_provenance(pool.table, fp_keys, fps,
        json.loads((ARTIFACT_DIR / f"pool_fp_provenance{POOL_TAG}.json").read_text()))
    fp_models, fp_lineage = fit_fingerprints(store, labels, splits, {**protocol, **execution}, protocol_path, ARTIFACT_DIR, device)
    cfg = EngineConfig(generate=True, gen_consensus=True, gen_expand=True, gen_max_analogs=30, gen_cap=900, np_score=True)
    library = LibraryIndex(store, row_keys, cfg.library)
    analog = AnalogIndex(store, row_keys, cfg.analog, cache=None)
    analog_status = analog.fresh_index_provenance()
    memo = SearchMemo(analog)
    keys = MetricKeys(ARTIFACT_DIR / "external_metric_keys.parquet")
    predictor = Predictor(BaselineEngine(pool, library, cfg, memo, fps), pool, keys,
                          fp_models=list(fp_models.values()), device=device)
    predictor.ranker = type("FeatureDeclaration", (), {"columns": ["d4_int"]})()
    source_rows = {name: np.flatnonzero((meta.ingest_lib == name).to_numpy()) for name in previous_protocol["query_libraries"]}

    def prepare(item, external=False):
        if external:
            spectra = [{"mz": np.asarray(r.ms2_mzs), "it": np.asarray(r.ms2_normalized_intensities),
                        "precursor_mz": float(r.precursor_mz), "adduct": str(r.adduct), "mode": str(r.ionization_mode), "ce_ev": 0.0}
                       for r in item["records"].itertuples()]
            query = Query("external:" + item["metric"], item["key"], spectra, np.empty(0, dtype=np.int64))
            smiles = item["smiles"]
        else:
            available = groups[item["key"]]
            rng = np.random.default_rng(int(hashlib.sha256(item["key"].encode()).hexdigest()[:16], 16))
            chosen = np.sort(rng.choice(available, 5, replace=False)) if len(available) > 5 else available
            spectra = []
            for index in chosen:
                mz, intensity = store.peaks(index); row = meta.iloc[index]
                spectra.append({"mz": mz, "it": intensity, "precursor_mz": float(row.precursor_mz),
                                "adduct": str(row.adduct), "mode": str(row.ionization_mode), "ce_ev": float(row.ce_ev)})
            query = Query(item["library"] + ":" + item["key"], item["key"], spectra, chosen)
            smiles = labels[item["key"]]
        truth, same = equivalent_keys(item["key"], molecule_mass(spectra), truth_pool, keys, truth_smiles=smiles)
        candidate_truth, candidates = equivalent_keys(item["key"], molecule_mass(spectra), pool, keys, truth_smiles=smiles)
        if truth != item["metric"] or candidate_truth != truth:
            raise ValueError("Population metric and query truth disagree")
        same.update(candidates); same.add(item["key"])
        if not external:
            include_equivalent_source_rows(query, same, source_rows[item["library"]], row_keys)
        elif any(key in library.rows_of for key in same):
            raise ValueError("A verified external truth unexpectedly has training reference spectra")
        return query, same, truth

    def frames(query, same, truth, classes, mode):
        spectra = perturb_spectra(query.spectra, mode=mode, seed=seed, query_id=query.qid)
        for cls in classes:
            knockout = knockout_for(query, cls, library.rows_of, same)
            union = predictor.features(spectra, excluded_rows=knockout.excluded_rows,
                excluded_keys=knockout.excluded_keys, knocked_keys=same if cls != "C1" else frozenset())
            if not union.empty:
                union = union.assign(qid=query.qid, cls=cls)
                union["label"] = [int(keys.get(key, pool.smiles_of(key)) == truth) for key in union.pool_key]
            yield cls, {policy: budget_features(union, policy) for policy in POLICIES}

    tables, development = {}, {p: [] for p in POLICIES}
    for split in ("train", "dev"):
        blocks = {p: [] for p in POLICIES}
        for index, item in enumerate(splits[split]):
            memo.reset(); query, same, truth = prepare(item)
            predictor.fp_models = [fp_models[1 - item["fold"]]] if split == "train" else list(fp_models.values())
            for cls, budgets in frames(query, same, truth, CLASSES, "none"):
                for policy, frame in budgets.items():
                    if len(frame): blocks[policy].append(frame)
                    if split == "dev" and cls == "C3":
                        development[policy].append(int(frame.label.max()) if len(frame) else 0)
            if (index + 1) % 10 == 0:
                print(f"Prepared {split} budget features {index + 1}/{len(splits[split])}", flush=True)
        tables[split] = {p: pd.concat(blocks[p], ignore_index=True) for p in POLICIES}
    choice = development_choice(development["legacy"], development["expanded"])
    if choice["lost"]:
        raise ValueError("Legacy-preservation invariant failed on development")
    models, training_audit = {}, {}
    for policy in POLICIES:
        columns = [c for c in feature_columns(tables["train"][policy]) if c in tables["dev"][policy] and not c.startswith(("prior_", "audit_"))]
        base = RankerConfig()
        config = RankerConfig(seeds=tuple(execution["ranker"]["seeds"]), num_boost_round=execution["ranker"]["rounds"],
            early_stopping=execution["ranker"]["early_stopping"],
            params={**base.params, "num_threads": 2, "deterministic": True, "force_col_wise": True})
        model = Ranker(config).fit(tables["train"][policy], tables["dev"][policy], columns)
        for index, fitted in enumerate(model.models):
            fitted.save_model(str(ARTIFACT_DIR / "models" / f"external_{policy}_{index}.txt"))
        models[policy] = model
        training_audit[policy] = {"columns": columns, "best_iterations": model.best_iters, "config": asdict(config),
            "positive_training_groups": int(tables["train"][policy].groupby(["qid", "cls"]).label.max().sum())}
        print(f"Fitted prespecified budget ranker {policy}", flush=True)
    model_hashes = {p.name: sha256_file(p) for p in sorted((ARTIFACT_DIR / "models").iterdir())}
    if len(model_hashes) != 8:
        raise ValueError("Unexpected fresh learned artifact inventory")
    seal = {"status": "ALL_MODELS_AND_DEVELOPMENT_CHOICE_SEALED_BEFORE_TEST_FEATURES", "sealed_at": datetime.now(UTC).isoformat(),
            "model_files": model_hashes, "models": training_audit, "fingerprint_training": fp_lineage,
            "development_choice": choice, "partitions": split_audit, "external_population": external_population,
            "protocol_sha256": sha256_file(protocol_path), "execution_protocol_sha256": sha256_file(execution_path),
            "source_protocol_sha256": sha256_file(source_protocol_path)}
    seal_hash = digest(seal)
    seal_path = ARTIFACT_DIR / "external_training_seal.json"
    seal_path.write_text(json.dumps(seal, indent=2, allow_nan=False))
    del tables
    print("Sealed all weights, populations and development choice; starting test features", flush=True)
    predictor.fp_models = list(fp_models.values())
    ranks, coverage, sizes, private_ids = {}, {}, {}, {}
    for cohort, items, classes, modes in (("internal", splits["test"], CLASSES, ("none",)),
                                          ("external", external_items, ("C2", "C3"), ("none", "blank"))):
        names = (*POLICIES, "mass_only", "d1", "d1_d2")
        ranks[cohort] = {p: {m: {c: [] for c in classes} for m in modes} for p in names}
        coverage[cohort] = {p: {c: [] for c in classes} for p in POLICIES}
        sizes[cohort] = {p: {c: [] for c in classes} for p in POLICIES}
        private_ids[cohort] = []
        for index, item in enumerate(items):
            memo.reset(); query, same, truth = prepare(item, external=cohort == "external")
            private_ids[cohort].append(query.qid)
            def ranked(frame, order):
                guesses = keys.dedupe_ranked(frame.pool_key.to_numpy()[order].tolist()[:60], pool.smiles_of, 25)
                return int(rank_of_truth([metric for _, metric in guesses], truth))
            for mode in modes:
                for cls, budgets in frames(query, same, truth, classes, mode):
                    for policy, frame in budgets.items():
                        value = ranked(frame, np.argsort(-models[policy].predict(frame), kind="stable")) if len(frame) else 0
                        ranks[cohort][policy][mode][cls].append(value)
                        if mode == "none":
                            coverage[cohort][policy][cls].append(int(frame.label.max()) if len(frame) else 0)
                            sizes[cohort][policy][cls].append(len(frame))
                    legacy = budgets["legacy"]
                    fixed = legacy[legacy.gen == 0].reset_index(drop=True) if len(legacy) else legacy
                    for name in ("mass_only", "d1", "d1_d2"):
                        value = 0
                        if len(fixed):
                            if name == "mass_only": order = np.argsort(np.abs(fixed.ppm), kind="stable")
                            else:
                                primary, secondary = fixed.d1_max.to_numpy(), fixed.d1_mean.to_numpy()
                                if name == "d1_d2":
                                    primary = np.maximum(primary, fixed.an_p4.fillna(0).to_numpy()) if "an_p4" in fixed else primary
                                    secondary = fixed.an_p2.fillna(0).to_numpy() if "an_p2" in fixed else secondary
                                order = np.lexsort((np.abs(fixed.ppm.to_numpy()), -secondary, -primary))
                            value = ranked(fixed, order)
                        ranks[cohort][name][mode][cls].append(value)
            if (index + 1) % 10 == 0:
                print(f"Completed sealed {cohort} test {index + 1}/{len(items)}", flush=True)
    if digest(json.loads(seal_path.read_text())) != seal_hash or model_hashes != {p.name: sha256_file(p) for p in sorted((ARTIFACT_DIR / "models").iterdir())}:
        raise ValueError("Learned artifacts or seal changed during test")
    verify_manifest(ARTIFACT_DIR, manifest)
    rr = lambda v: np.where(np.asarray(v) > 0, 1 / np.maximum(np.asarray(v), 1), 0.)
    effects = {"internal": {}, "external": {}}
    for cohort in effects:
        for cls in coverage[cohort]["legacy"]:
            effects[cohort][f"{cohort}_mrr_expanded_minus_legacy:{cls}"] = rr(ranks[cohort]["expanded"]["none"][cls]) - rr(ranks[cohort]["legacy"]["none"][cls])
            if cohort == "external":
                effects[cohort][f"external_coverage_expanded_minus_legacy:{cls}"] = np.asarray(coverage[cohort]["expanded"][cls]) - np.asarray(coverage[cohort]["legacy"][cls])
                for policy in POLICIES:
                    effects[cohort][f"external_{policy}_blank_minus_none:{cls}"] = rr(ranks[cohort][policy]["blank"][cls]) - rr(ranks[cohort][policy]["none"][cls])
    if {name for cohort in effects.values() for name in cohort} != set(execution["secondary_family"]):
        raise ValueError("Prespecified secondary family differs")
    summaries = {cohort: {p: {mode: {cls: asdict(summarize(values, n_boot=4000, seed=seed)) for cls, values in classes.items()}
        for mode, classes in conditions.items()} for p, conditions in policies.items()} for cohort, policies in ranks.items()}
    primary = primary_coverage_interval(coverage["internal"]["legacy"]["C3"], coverage["internal"]["expanded"]["C3"], seed=seed)
    if primary["lost"]:
        raise ValueError("Legacy-preserving expansion lost internal test truths")
    result = {"schema_version": 1, "status": "COMPLETE_C3_BUDGET_AND_EXTERNAL_VALIDATION", "checked_at": datetime.now(UTC).isoformat(),
        "elapsed_seconds": round(time.time() - started, 1), "protocol": protocol, "execution_protocol": execution,
        "independence_audit": independence, "partitions": split_audit, "external_population": external_population,
        "training_seal": seal, "training_seal_sha256": seal_hash, "model_bytes_unchanged_after_test": True,
        "source_manifest": manifest, "fingerprint_provenance": fp_status, "analog_index_provenance": analog_status,
        "engine_config": asdict(cfg), "models": summaries, "primary": primary,
        "secondary_joint_intervals": joint_cohort_intervals(effects, n_boot=4000, seed=seed),
        "candidate_coverage": {cohort: {p: {cls: {"n_molecules": len(v), "covered": int(sum(v)),
            "oracle_mrr_ceiling": float(np.mean(v)), "mean_candidates": float(np.mean(sizes[cohort][p][cls])),
            "maximum_candidates": int(max(sizes[cohort][p][cls]))} for cls, v in classes.items()} for p, classes in policies.items()} for cohort, policies in coverage.items()},
        "private_details_sha256": digest({"ranks": ranks, "coverage": coverage, "ids": private_ids, "development": development}),
        "environment": {"python": sys.version.split()[0], "device": device,
            **{name: importlib.metadata.version(name) for name in ("rdkit", "numpy", "pandas", "pyarrow", "torch", "lightgbm")}},
        "limitations": ["External standard/metric novelty is computationally audited against the pinned training snapshot; source annotations are not independently authenticated laboratory identity.",
            "External collision energy is not harmonized to eV; unavailable numeric eV is encoded zero. This is a distribution-shift limitation.",
            "External libraries have a restricted chemical/instrument distribution and are not a representative clinical or environmental population.",
            "Generic-family supervised separation does not remove task-defined shared candidate knowledge or fixed NP-score priors.",
            "Both predeclared budgets and every adverse outcome are retained. Positive generator improvement is not guaranteed.",
            "Historical study populations, fitting and architecture differ; direct score comparison does not isolate improvement.",
            "No competition submission, paid infrastructure, service deployment or raw competition-record export."]}
    with Path(os.environ["CASMI_EXTERNAL_STUDY_OUTPUT"]).open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write("\n")
    print("Completed sealed C3 and external study with all declared outcomes", flush=True)


if __name__ == "__main__": main()
