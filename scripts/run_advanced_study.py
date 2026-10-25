"""Execute the predeclared internal family-disjoint study; export aggregates only.

All supervised weights are fitted from scratch. No archived model is loaded.
The test phase begins only after every ranker and its feature list is sealed.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import sys
import time
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]


def ablation_columns(columns, name):
    """Remove a channel and its derived context together; retain fixed candidates."""
    membership = ("prior_",)
    if name == "membership_only":
        selected = [c for c in columns if c.startswith(membership)]
    else:
        prefixes = {"full": (), "no_d1": ("d1_",), "no_d2": ("an_",),
                    "no_d3": ("d3_",), "no_d4": ("d4_", "d4x_"),
                    "no_gen_features": ("gen",), "no_np_score": ("np_score",)}
        if name not in prefixes:
            raise ValueError("Unknown prespecified ablation")
        selected = [c for c in columns if not c.startswith(membership + prefixes[name])]
    if not selected:
        raise ValueError("Ablation has no usable features")
    return selected


class SearchMemo:
    """Per-query memoization keyed by full cleaned evidence and reference exclusions."""
    def __init__(self, index):
        self.index, self.cache = index, {}
        self.hits = self.misses = 0

    def reset(self):
        self.cache.clear()

    def search(self, spectra, mass, excluded_rows=None):
        import numpy as np
        key = hashlib.sha256()
        key.update(float(mass).hex().encode())
        for s in spectra:
            key.update(str((s["adduct"], s["mode"], float(s["precursor_mz"]).hex())).encode())
            for field in ("mz_c", "it_c"):
                data = np.asarray(s[field], dtype=np.float64)
                key.update(len(data).to_bytes(8, "little"))
                key.update(data.tobytes())
        key.update(np.unique(np.asarray(excluded_rows if excluded_rows is not None else [], dtype=np.int64)).tobytes())
        token = key.hexdigest()
        if token not in self.cache:
            self.cache[token] = self.index.search(spectra, mass, excluded_rows)
            self.misses += 1
        else:
            self.hits += 1
        return dict(self.cache[token])


def main():
    import numpy as np
    import pandas as pd
    import rdkit
    import torch
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE, POOL_TAG
    from casmi.channels.library import LibraryIndex
    from casmi.channels.analog import AnalogIndex
    from casmi.channels.fpmodel import FPModel, FPModelConfig, save
    from casmi.channels.fpinfer import FP_CLEAN, MAX_PEAKS, ADDUCT_ID
    from casmi.chem.keys import MetricKeys
    from casmi.chem.canon import canonical_ik14
    from casmi.spectra.clean import clean_peaks_fast
    from casmi.spectra.store import SpectrumStore
    from casmi.pool.pool import CandidatePool
    from casmi.pool.fingerprint_provenance import verify_fingerprint_provenance
    from casmi.engine import BaselineEngine, EngineConfig, molecule_mass
    from casmi.predict import Predictor
    from casmi.rank.features import feature_columns
    from casmi.rank.ranker import Ranker, RankerConfig
    from casmi.validate.knockout import Query, CLASSES, knockout_for, include_equivalent_source_rows, equivalent_keys
    from casmi.validate.metric import rank_of_truth, summarize
    from casmi.validate.paired_bootstrap import paired_difference
    from casmi.validate.scaffold_study import family_key, partition, digest, joint_intervals
    from casmi.validate.spectral_controls import perturb_spectra
    from casmi.validate.run_alignment import query_signature
    from casmi.validate.source_manifest import create_manifest, required_validation_inputs, verify_manifest, sha256_file

    started = time.time()
    protocol_path = ROOT / "protocols/advanced-internal-v1.json"
    protocol = json.loads(protocol_path.read_text())
    if rdkit.__version__ != "2026.03.3":
        raise RuntimeError("Official metric RDKit is required")
    seed = protocol["seed"]
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(2)
    torch.use_deterministic_algorithms(True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    store = SpectrumStore(ARTIFACT_DIR / "store")
    meta = store.meta
    key_of_row = meta.train_ik14.astype(str).to_numpy()
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    truth_pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "structures.parquet"))
    labels = dict(zip(truth_pool.key, truth_pool.smiles))
    source_quality = json.loads((ARTIFACT_DIR / f"pool_source_quality{POOL_TAG}.json").read_text())
    static_names = list(required_validation_inputs(pool_tag=POOL_TAG, analog=True))
    static_names += ["structures.parquet", f"pool_source_quality{POOL_TAG}.json"]
    manifest = create_manifest(ARTIFACT_DIR, static_names)
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE, mmap_mode="r", allow_pickle=False)
    fp_keys = pd.read_parquet(ARTIFACT_DIR / POOL_FP_KEYS_FILE).pool_key.to_numpy()
    fp_status = verify_fingerprint_provenance(pool.table, fp_keys, fps,
        json.loads((ARTIFACT_DIR / f"pool_fp_provenance{POOL_TAG}.json").read_text()))
    cfg = EngineConfig(generate=True, gen_consensus=True, np_score=True)
    library = LibraryIndex(store, key_of_row, cfg.library)
    raw_analog = AnalogIndex(store, key_of_row, cfg.analog, cache=None)
    analog_status = raw_analog.fresh_index_provenance()
    memo = SearchMemo(raw_analog)
    keys = MetricKeys(ARTIFACT_DIR / "advanced_metric_keys.parquet")
    print("Prepared fresh artifact inventory and D2 index", flush=True)

    # Population construction is label/metadata-only and cannot consult scores.
    eligible = (meta.ingest_lib.isin(protocol["query_libraries"]).to_numpy()
                & meta.test_adduct.to_numpy() & (meta.n_peaks.to_numpy() >= 5)
                & (np.abs(meta.ppm.to_numpy()) <= 20))
    rows = np.flatnonzero(eligible)
    choices = pd.DataFrame({"row": rows, "key": key_of_row[rows],
                            "library": meta.ingest_lib.astype(str).to_numpy()[rows]})
    mass_of = dict(zip(truth_pool.key, truth_pool.mass))
    choices = choices[choices.key.map(mass_of).between(*protocol["query_mass_range"])]
    counts = choices.groupby(["key", "library"]).size().reset_index(name="count")
    best = counts.sort_values(["count", "key", "library"], ascending=[False, True, True]).drop_duplicates("key")
    groups = {key: group.row.to_numpy() for key, group in
              choices.merge(best[["key", "library"]], on=["key", "library"]).groupby("key")}
    excluded_families = {family_key(labels[key]) for key in pd.unique(key_of_row[
        (meta.ingest_lib == "enveda-np-examples").to_numpy()]) if key in labels}
    representatives = {}
    invalid_families = blocked_families = 0
    for row in best.itertuples():
        family = family_key(labels[row.key])
        if family is None:
            invalid_families += 1
        elif family in excluded_families:
            blocked_families += 1
        elif family not in representatives or row.key < representatives[family]["key"]:
            representatives[family] = {"key": row.key, "family": family, "library": row.library}
    records = []
    # Only family representatives need expensive tautomer canonicalization.
    for value in representatives.values():
        metric = canonical_ik14(labels[value["key"]])
        if metric:
            records.append({**value, "metric": metric})
    splits = partition(records, protocol["split_counts"], protocol["split_salt"])
    split_audit = {name: {"n_molecules": len(values), "identity_sha256": digest([r["key"] for r in values]),
                          "family_sha256": digest([r["group"] for r in values]),
                          "metric_sha256": digest([r["metric"] for r in values]),
                          "libraries": pd.Series([r["library"] for r in values]).value_counts().to_dict()}
                   for name, values in splits.items()}
    split_audit["cross_split_identity_overlaps"] = 0
    split_audit["invalid_family_labels_excluded"] = invalid_families
    split_audit["historical_vin_family_rows_excluded"] = blocked_families
    split_audit["eligible_representative_families"] = len(records)
    print("Prepared declared 600/150/250 family-disjoint population", flush=True)

    # Both D3 models train on selected training identities only, on opposite folds.
    train_lookup = {row["key"]: row for row in splits["train"]}
    fp_eligible = (meta.test_adduct.to_numpy() & (meta.n_peaks.to_numpy() >= 5)
                   & (np.abs(meta.ppm.to_numpy()) <= 20))
    training_rows = {}
    for index in np.flatnonzero(fp_eligible & np.isin(key_of_row, list(train_lookup))):
        training_rows.setdefault(key_of_row[index], []).append(int(index))
    target_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_models, fp_lineage = {}, []
    settings = protocol["fp_model"]
    (ARTIFACT_DIR / "models").mkdir(exist_ok=True)
    for fold in (0, 1):
        rng = np.random.default_rng(seed + fold)
        selected_rows, selected_keys = [], []
        for key, item in sorted(train_lookup.items()):
            if item["fold"] != fold:
                continue
            available = np.asarray(training_rows.get(key, []), dtype=np.int64)
            cap = protocol["max_training_spectra_per_structure"]
            if len(available) > cap:
                available = np.sort(rng.choice(available, cap, replace=False))
            selected_rows.extend(available.tolist())
            selected_keys.extend([key] * len(available))
        if len(selected_rows) < 100:
            raise ValueError("Insufficient fingerprint training spectra in declared fold")
        if not set(selected_keys) <= {r["key"] for r in splits["train"]}:
            raise ValueError("Fingerprint training crossed the declared training partition")
        n = len(selected_rows)
        mz = np.zeros((n, MAX_PEAKS), np.float32)
        it = np.zeros_like(mz)
        npk = np.zeros(n, np.int64)
        unique_keys = sorted(set(selected_keys))
        targets = {}
        for key in unique_keys:
            mol = Chem.MolFromSmiles(labels[key])
            bits = np.zeros(2048, np.uint8)
            DataStructs.ConvertToNumpyArray(target_gen.GetFingerprint(mol), bits)
            if not bits.any():
                raise ValueError("Invalid fingerprint training target")
            targets[key] = bits
        for j, index in enumerate(selected_rows):
            cm, ci = clean_peaks_fast(*store.peaks(index), float(meta.precursor_mz.iat[index]), FP_CLEAN)
            mz[j, :len(cm)], it[j, :len(ci)], npk[j] = cm, ci, len(cm)
        side = meta.iloc[selected_rows]
        inputs = {"mz": mz, "it": it, "npk": npk, "prec": side.precursor_mz.to_numpy(np.float32),
                  "adduct": side.adduct.map(ADDUCT_ID).to_numpy(np.int64),
                  "ce": side.ce_ev.to_numpy(np.float32),
                  "positive": (side.ionization_mode == "positive").to_numpy(np.float32)}
        y = np.asarray([targets[key] for key in selected_keys], dtype=np.float32)
        torch.manual_seed(seed + fold)
        model = FPModel(FPModelConfig(d_model=settings["d_model"], n_layers=settings["n_layers"],
                                      d_ff=settings["d_ff"])).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=settings["lr"], weight_decay=.01)
        losses = []
        model.train()
        for epoch in range(settings["epochs"]):
            total = 0.0
            order = rng.permutation(n)
            for start in range(0, n, settings["batch_size"]):
                ids = order[start:start + settings["batch_size"]]
                x = {name: torch.as_tensor(np.ascontiguousarray(value[ids]), device=device,
                     dtype=torch.long if name in ("npk", "adduct") else torch.float32) for name, value in inputs.items()}
                target = torch.as_tensor(y[ids], device=device)
                loss = torch.nn.functional.binary_cross_entropy_with_logits(model(**x), target)
                if not torch.isfinite(loss):
                    raise ValueError("Nonfinite fingerprint training loss")
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                total += float(loss.detach().cpu()) * len(ids)
            losses.append(total / n)
            print(f"Trained fresh fingerprint fold {fold}, epoch {epoch + 1}/{settings['epochs']}", flush=True)
        model.eval()
        model_path = ARTIFACT_DIR / "models" / f"advanced_fp_fold{fold}.pt"
        lineage = {"fold": fold, "n_spectra": n, "n_structures": len(unique_keys),
                   "training_rows_sha256": digest(selected_rows), "training_identities_sha256": digest(unique_keys),
                   "training_families_sha256": digest(sorted({train_lookup[k]["group"] for k in unique_keys})),
                   "dev_test_identities_used": 0, "epochs_fixed_before_execution": settings["epochs"],
                   "training_bce_by_epoch": losses, "initialization": "fresh seeded random; no archived weights"}
        save(model, model_path, {"lineage": lineage, "protocol_sha256": sha256_file(protocol_path)})
        fp_lineage.append({**lineage, "model_sha256": sha256_file(model_path)})
        fp_models[fold] = model
        del mz, it, npk, inputs, y

    predictor = Predictor(BaselineEngine(pool, library, cfg, memo, fps), pool, keys,
                          fp_models=list(fp_models.values()), device=device)
    # Trigger the actual inference D4 implementation without fitting any archived model.
    predictor.ranker = type("FeatureDeclaration", (), {"columns": ["d4_int"]})()
    source_rows = {name: np.flatnonzero((meta.ingest_lib == name).to_numpy()) for name in protocol["query_libraries"]}
    private_query_hashes = {name: [] for name in splits}
    truth_by_qid, query_context = {}, {}

    def prepare_query(item):
        available = groups[item["key"]]
        rng = np.random.default_rng(int(hashlib.sha256(item["key"].encode()).hexdigest()[:16], 16))
        chosen = np.sort(rng.choice(available, protocol["max_spectra"], replace=False)) if len(available) > protocol["max_spectra"] else available
        spectra = []
        for index in chosen:
            mz, intensity = store.peaks(index)
            row = meta.iloc[index]
            spectra.append({"mz": mz, "it": intensity, "precursor_mz": float(row.precursor_mz),
                            "adduct": str(row.adduct), "ce_ev": float(row.ce_ev), "mode": str(row.ionization_mode)})
        query = Query(f"{item['library']}:{item['key']}", item["key"], spectra, chosen)
        truth, same = equivalent_keys(item["key"], molecule_mass(spectra), truth_pool, keys,
                                       truth_smiles=labels[item["key"]])
        candidate_truth, candidates = equivalent_keys(item["key"], molecule_mass(spectra), pool, keys,
                                                       truth_smiles=labels[item["key"]])
        if truth != item["metric"] or candidate_truth != truth:
            raise ValueError("Split labels and inference canonical truth disagree")
        same.update(candidates)
        same.add(item["key"])
        include_equivalent_source_rows(query, same, source_rows[item["library"]], key_of_row)
        truth_by_qid[query.qid] = truth
        m = molecule_mass(spectra)
        unique = lambda field: next(iter(values)) if len(values := {s[field] for s in spectra}) == 1 else "mixed"
        query_context[query.qid] = {"library": item["library"], "adduct": unique("adduct"), "mode": unique("mode"),
                                  "mass": "150-300" if m < 300 else "300-600" if m < 600 else "600-1170",
                                  "spectra": "1" if len(spectra) == 1 else "2-3" if len(spectra) < 4 else "4-5",
                                  "instrument": "timsTOF" if meta.is_timstof.iloc[chosen].all() else "other_or_mixed"}
        return query, same

    def frames(query, same, mode):
        spectra = perturb_spectra(query.spectra, mode=mode, seed=seed, query_id=query.qid)
        for cls in CLASSES:
            knockout = knockout_for(query, cls, library.rows_of, same)
            frame = predictor.features(spectra, excluded_rows=knockout.excluded_rows,
                                       excluded_keys=knockout.excluded_keys,
                                       knocked_keys=same if cls != "C1" else frozenset())
            if not frame.empty:
                frame = frame.assign(qid=query.qid, cls=cls)
                frame["label"] = [int(keys.get(key, pool.smiles_of(key)) == truth_by_qid[query.qid]) for key in frame.pool_key]
            yield cls, frame

    tables = {}
    for split in ("train", "dev"):
        blocks = []
        for index, item in enumerate(splits[split]):
            memo.reset()
            query, same = prepare_query(item)
            private_query_hashes[split].append(query_signature(query))
            # No training spectrum may be scored by its own fitted D3 fold.
            predictor.fp_models = [fp_models[1 - item["fold"]]] if split == "train" else list(fp_models.values())
            blocks.extend(frame for _, frame in frames(query, same, "none") if not frame.empty)
            if (index + 1) % 25 == 0:
                print(f"Prepared {split} features {index + 1}/{len(splits[split])}", flush=True)
        tables[split] = pd.concat(blocks, ignore_index=True)
    columns = [name for name in feature_columns(tables["train"]) if name in tables["dev"]]
    models, training_audit = {}, {}
    for name in protocol["ablations"]:
        selected = ablation_columns(columns, name)
        base = RankerConfig()
        config = RankerConfig(seeds=tuple(protocol["ranker"]["seeds"]),
            num_boost_round=protocol["ranker"]["rounds"], early_stopping=protocol["ranker"]["early_stopping"],
            params={**base.params, "num_threads": 2, "deterministic": True, "force_col_wise": True})
        model = Ranker(config)
        model.fit(tables["train"], tables["dev"], selected)
        model_hashes = []
        for index, fitted in enumerate(model.models):
            path = ARTIFACT_DIR / "models" / f"advanced_{name}_{index}.txt"
            fitted.save_model(str(path))
            model_hashes.append(sha256_file(path))
        models[name] = model
        training_audit[name] = {"columns": selected, "best_iterations": model.best_iters,
                                "model_sha256": model_hashes, "config": asdict(config),
                                "positive_training_groups": int(tables["train"].groupby(["qid", "cls"]).label.max().sum()),
                                "zero_positive_groups_retained_in_test_denominator": True}
        print(f"Fitted prespecified ranker {name}", flush=True)
    model_files = {path.name: sha256_file(path) for path in sorted((ARTIFACT_DIR / "models").iterdir())}
    if len(model_files) != 2 + len(protocol["ablations"]) * len(protocol["ranker"]["seeds"]):
        raise ValueError("Unexpected model artifact inventory; fresh training required")
    seal = {"status": "ALL_MODELS_SEALED_BEFORE_TEST_FEATURES", "sealed_at": datetime.now(UTC).isoformat(),
            "model_files": model_files,
            "models": training_audit, "fingerprint_training": fp_lineage,
            "query_input_sha256": {name: digest(values) for name, values in private_query_hashes.items() if name != "test"},
            "test_identities_sha256": split_audit["test"]["identity_sha256"], "protocol_sha256": sha256_file(protocol_path)}
    seal_sha256 = digest(seal)
    (ARTIFACT_DIR / "training_seal.json").write_text(json.dumps(seal, indent=2, allow_nan=False))
    del tables
    print("Sealed all learned models; starting final test phase", flush=True)

    names = ["mass_only", "d1", "d1_d2"] + list(models)
    ranks = {name: {mode: {cls: [] for cls in CLASSES} for mode in protocol["controls"]} for name in names}
    coverage = {cls: {"pool_truth_present": 0, "with_generation_truth_present": 0} for cls in CLASSES}
    test_ids = []
    predictor.fp_models = list(fp_models.values())

    def ranked(frame, order, truth):
        guesses = keys.dedupe_ranked(frame.pool_key.to_numpy()[order].tolist()[:60], pool.smiles_of, 25)
        return int(rank_of_truth([metric for _, metric in guesses], truth))

    for index, item in enumerate(splits["test"]):
        memo.reset()
        query, same = prepare_query(item)
        test_ids.append(query.qid)
        private_query_hashes["test"].append(query_signature(query))
        for mode in protocol["controls"]:
            for cls, frame in frames(query, same, mode):
                if frame.empty:
                    for name in names:
                        ranks[name][mode][cls].append(0)
                    continue
                truth = truth_by_qid[query.qid]
                fixed = frame[frame.gen == 0].reset_index(drop=True)
                for name in ("mass_only", "d1", "d1_d2"):
                    if fixed.empty:
                        value = 0
                    elif name == "mass_only":
                        value = ranked(fixed, np.argsort(np.abs(fixed.ppm), kind="stable"), truth)
                    else:
                        primary = fixed.d1_max.to_numpy()
                        secondary = fixed.d1_mean.to_numpy()
                        if name == "d1_d2":
                            primary = np.maximum(primary, fixed.an_p4.fillna(0).to_numpy()) if "an_p4" in fixed else primary
                            secondary = fixed.an_p2.fillna(0).to_numpy() if "an_p2" in fixed else secondary
                        value = ranked(fixed, np.lexsort((np.abs(fixed.ppm.to_numpy()), -secondary, -primary)), truth)
                    ranks[name][mode][cls].append(value)
                for name, model in models.items():
                    ranks[name][mode][cls].append(ranked(frame, np.argsort(-model.predict(frame), kind="stable"), truth))
                if mode == "none":
                    coverage[cls]["pool_truth_present"] += int(fixed.label.max()) if len(fixed) else 0
                    coverage[cls]["with_generation_truth_present"] += int(frame.label.max())
        if (index + 1) % 10 == 0:
            print(f"Completed sealed test {index + 1}/{len(splits['test'])} with controls", flush=True)
    if digest(json.loads((ARTIFACT_DIR / "training_seal.json").read_text())) != seal_sha256:
        raise ValueError("Model seal changed during test")
    if {path.name: sha256_file(path) for path in sorted((ARTIFACT_DIR / "models").iterdir())} != model_files:
        raise ValueError("Learned artifact bytes changed during test")
    verify_manifest(ARTIFACT_DIR, manifest)
    summaries = {name: {mode: {cls: asdict(summarize(values, n_boot=protocol["bootstrap_draws"], seed=seed))
                               for cls, values in classes.items()} for mode, classes in conditions.items()}
                 for name, conditions in ranks.items()}
    rr = lambda values: np.where(np.asarray(values) > 0, 1.0 / np.maximum(np.asarray(values), 1), 0.0)
    effects = {}
    for name in protocol["ablations"]:
        if name != "full":
            for cls in CLASSES:
                effects[f"{name}_minus_full:{cls}"] = rr(ranks[name]["none"][cls]) - rr(ranks["full"]["none"][cls])
    for mode in protocol["controls"][1:]:
        for cls in CLASSES:
            effects[f"full_{mode}_minus_none:{cls}"] = rr(ranks["full"][mode][cls]) - rr(ranks["full"]["none"][cls])
    primary = paired_difference(ranks["d1_d2"]["none"]["C2"], ranks["full"]["none"]["C2"],
                                n_boot=protocol["bootstrap_draws"], seed=seed)
    slices = {}
    for dimension in ("library", "adduct", "mode", "mass", "spectra", "instrument"):
        slices[dimension] = {}
        for value in sorted({query_context[qid][dimension] for qid in test_ids}):
            ids = [i for i, qid in enumerate(test_ids) if query_context[qid][dimension] == value]
            if len(ids) < protocol["slice_min_n"]:
                slices[dimension][value] = {"n": len(ids), "status": "SUPPRESSED_SMALL_CELL"}
            else:
                cells = {}
                for cls in CLASSES:
                    a = np.asarray(ranks["d1_d2"]["none"][cls])[ids]
                    b = np.asarray(ranks["full"]["none"][cls])[ids]
                    delta, low, high = paired_difference(a, b, n_boot=protocol["bootstrap_draws"], seed=seed)
                    cells[cls] = {"baseline_mrr": float(rr(a).mean()), "full_mrr": float(rr(b).mean()),
                                  "delta_mrr": delta, "ci_low": low, "ci_high": high}
                slices[dimension][value] = {"n": len(ids), "status": "DESCRIPTIVE_ONLY", "classes": cells}
    result = {"schema_version": 1, "status": "COMPLETE_INTERNAL_FAMILY_HOLDOUT",
              "checked_at": datetime.now(UTC).isoformat(), "elapsed_seconds": round(time.time() - started, 1),
              "protocol": protocol, "protocol_sha256": sha256_file(protocol_path), "partitions": split_audit,
              "training_seal": seal, "training_seal_sha256": seal_sha256,
              "source_manifest": manifest, "candidate_source_quality": source_quality,
              "fingerprint_provenance": fp_status, "analog_index_provenance": analog_status,
              "search_memo": {"scope": "per-query full-cleaned-content plus exact exclusion set; no index file cache",
                              "hits": memo.hits, "misses": memo.misses},
              "engine_config": asdict(cfg), "models": summaries, "candidate_coverage": coverage,
              "primary": {"contrast": protocol["primary"], "delta_mrr": primary[0], "ci_low": primary[1], "ci_high": primary[2],
                          "interpretation": "Single predeclared paired diagnostic; no automatic promotion"},
              "secondary_joint_intervals": joint_intervals(effects, n_boot=protocol["bootstrap_draws"], seed=seed),
              "error_slices": slices, "private_details_sha256": digest({"ranks": ranks, "test_ids": test_ids, "signatures": private_query_hashes}),
              "environment": {"python": sys.version.split()[0], "device": device,
                              **{name: importlib.metadata.version(name) for name in ("rdkit", "numpy", "pandas", "pyarrow", "scipy", "numba", "torch", "lightgbm")}},
              "limitations": ["Internal source data, not newly acquired external spectra or an independent laboratory replication.",
                              "Family-disjoint supervised query fitting does not remove the task-defined retrieval reference corpus.",
                              "Generic Murcko topology and acyclic connectivity are computational grouping rules, not a guarantee of distant chemistry.",
                              "Fixed training sample and small fresh fingerprint architecture trade capacity for auditable zero-cost execution; no state-of-the-art claim.",
                              "Feature ablations reuse a fixed generated candidate universe; no_gen_features does not remove candidate generation.",
                              "Slices are descriptive; subgroup intervals do not establish multiplicity-adjusted subgroup effects.",
                              "No raw spectra, structures, query identifiers, per-query ranks, model weights or submitted predictions are exported."]}
    output = Path(os.environ.get("CASMI_ADVANCED_OUTPUT", "/kaggle/working/advanced_study.json"))
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print("Completed sealed internal family-holdout study with fresh weights and all declared controls", flush=True)


if __name__ == "__main__":
    main()
