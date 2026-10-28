"""Authenticate the purged-training study without promoting original novelty failures."""
from __future__ import annotations

from datetime import UTC, datetime
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts.build_purged_external_notebook import KERNEL, DATASETS, REGISTERED_COMMIT
from scripts.check_kaggle_publications import source_hash
from scripts.research_snapshot import snapshot_hashes
from casmi.validate.source_manifest import required_validation_inputs, _canonical_files_digest


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def text_hash(path):
    return hashlib.sha256((ROOT / path).read_text(encoding="utf-8").encode()).hexdigest()


def require(condition, message):
    if not condition: raise ValueError(message)


def valid_hash(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def close(a, b):
    return type(a) in (int, float) and math.isfinite(a) and math.isclose(a, b, abs_tol=1e-12)


def validate_result(result):
    expected = {"schema_version", "status", "checked_at", "elapsed_seconds", "protocol", "execution_protocol",
        "independence_audit", "purged_training", "partitions", "external_population", "training_seal", "training_seal_sha256",
        "model_bytes_unchanged_after_test", "source_manifest", "fingerprint_provenance", "analog_index_provenance",
        "engine_config", "models", "primary", "secondary_joint_intervals", "candidate_coverage", "private_details_sha256",
        "environment", "limitations", "raw_inputs", "source_code_sha256", "execution", "installed_packages"}
    require(set(result) == expected and result["schema_version"] == 1
        and result["status"] == "COMPLETE_C3_AND_PURGED_REFERENCE_EXTERNAL_HOLDOUT", "Incomplete or unexpected purged-study schema")
    protocol = json.loads((ROOT / "protocols/purged-external-v1.json").read_text())
    execution = json.loads((ROOT / "protocols/purged-external-v1-execution.json").read_text())
    original = json.loads((ROOT / "reports/advanced_study.json").read_text())
    prepared = json.loads((ROOT / "reports/gnps_expanded_source_preparation.json").read_text())
    require(result["protocol"] == protocol and result["execution_protocol"] == execution, "Registered design differs")
    require(result["source_code_sha256"] == snapshot_hashes(ROOT / "kaggle/audits/purged-external/notebook.ipynb"), "Executed source differs")
    require(valid_hash(result["private_details_sha256"]) and result["model_bytes_unchanged_after_test"] is True, "Private details or immutable model seal missing")
    seal, purge, partitions, external = (result[name] for name in ("training_seal", "purged_training", "partitions", "external_population"))
    require(digest(seal) == result["training_seal_sha256"]
        and seal["status"] == "ALL_MODELS_AND_DEVELOPMENT_CHOICE_SEALED_BEFORE_TEST_FEATURES", "Training seal altered")
    require(seal["protocol_sha256"] == text_hash("protocols/purged-external-v1.json")
        and seal["execution_protocol_sha256"] == text_hash("protocols/purged-external-v1-execution.json")
        and seal["source_protocol_sha256"] == text_hash("protocols/external-source-v3.json"), "Sealed protocol identities differ")
    require(seal["partitions"] == partitions and seal["external_population"] == external
        and datetime.fromisoformat(seal["sealed_at"]) <= datetime.fromisoformat(result["checked_at"]), "Sealed populations or temporal ordering differ")
    for name, count in protocol["internal_counts"].items():
        block = partitions[name]
        require(block["n_molecules"] == count and sum(block["libraries"].values()) == count, "Internal population denominator differs")
        require(all(valid_hash(block[field]) for field in ("identity_sha256", "family_sha256", "metric_sha256")), "Invalid internal population hash")
    require(partitions["previous_population_reconstructed"] is True and all(partitions[k] == 0 for k in
        ("previous_selected_metric_overlap", "previous_selected_family_overlap", "cross_split_identity_overlaps")), "Historical/current population leakage")
    n = external["n_molecules"]
    require(type(n) is int and 30 <= n <= 150 and external["supervised_family_overlap"] == 0
        and all(valid_hash(external[k]) for k in ("metric_sha256", "family_sha256", "accessions_sha256")), "External population contract differs")
    require(purge["status"] == "EXTERNAL_FAMILIES_PURGED_AND_INDEPENDENCE_VERIFIED"
        and purge["model_performance_inspected"] is False and purge["historical_population_reconstructed_before_purge"] is True,
        "Purge/historical reconstruction was not completed before fitting")
    require(purge["original_training"] == original["raw_inputs"]["train.parquet"]
        and purge["original_source_preparation"] == prepared and purge["protocol_sha256"] == seal["protocol_sha256"], "Original source lineage differs")
    require(purge["source_amendment_sha256"] == text_hash("protocols/purged-external-v1-source-amendment.json")
        and all(type(v) is int and v >= 0 for v in purge["family_grouping_quarantine"].values()), "Source grouping amendment or quarantine counts differ")
    require(purge["original_training_rows"] == 2539608 and type(purge["retained_training_rows"]) is int
        and purge["retained_training_rows"] >= 100000 and purge["removed_training_rows"] > 0
        and purge["retained_training_rows"] + purge["removed_training_rows"] == purge["original_training_rows"], "Purge row accounting differs")
    require(purge["selected_external_molecules"] == n and purge["selected_external_families"] >= n
        and purge["retained_external_family_overlap"] == 0 and valid_hash(purge["selected_metric_sha256"]), "External family purge failed")
    derived = purge["selected_source_preparation"]
    require(derived["source"] == prepared["source"] and derived["counts"]["eligible_metric_molecules"] == n
        and derived["counts"]["eligible_spectra"] == purge["selected_external_spectra_available"], "Selected public-source accounting differs")
    require(seal["purged_source_sha256"]["protected_population"] == purge["protected_population_sha256"], "Protected historical population hash differs")
    serialized_purge = json.dumps(purge, indent=2, allow_nan=False) + "\n"
    require(seal["purged_source_sha256"]["purge_report"] == hashlib.sha256(serialized_purge.encode()).hexdigest(), "Purge report differs from pre-test seal")
    for audit in (result["independence_audit"], purge["post_purge_independence"]):
        a = audit["identity_audit"]
        require(audit["status"] == "EXTERNAL_INDEPENDENCE_FEASIBLE" and audit["model_performance_inspected"] is False
            and audit["protocol_sha256"] == text_hash("protocols/external-source-v3.json"), "Purged snapshot independence not established")
        require(a["training_rows_scanned"] == purge["retained_training_rows"] and a["novel_external_molecules"] == n
            and a["novel_external_spectra"] == purge["selected_external_spectra_available"]
            and all(a[k] == 0 for k in ("unparseable_training_labels", "external_standard_overlap_molecules",
                "external_metric_overlap_molecules", "external_exact_content_overlap_molecules", "excluded_external_molecules_union")), "External identities/content remain in purged training")
        require(audit["inputs"]["selected_external.parquet"] == {"bytes": derived["parquet_bytes"], "sha256": derived["parquet_sha256"]}, "Selected source bytes differ")
    require(result["independence_audit"]["inputs"] == purge["post_purge_independence"]["inputs"], "Training data changed between purge and fitting")
    raw = result["raw_inputs"]
    require(set(raw) == {"original_train.parquet", "purged_train.parquet", "coconut.csv", "pubchemlite.csv", "gnps_cc0.parquet"}, "Raw source inventory differs")
    require(raw["original_train.parquet"] == original["raw_inputs"]["train.parquet"]
        and raw["purged_train.parquet"] == result["independence_audit"]["inputs"]["train.parquet"]
        and all(raw[k] == original["raw_inputs"][k] for k in ("coconut.csv", "pubchemlite.csv"))
        and raw["gnps_cc0.parquet"] == {"sha256": prepared["parquet_sha256"], "bytes": prepared["parquet_bytes"]}, "Raw input byte identities differ")
    files = seal["model_files"]
    expected_files = {f"external_fp_fold{i}.pt" for i in (0, 1)} | {f"external_{p}_{i}.txt" for p in ("legacy", "expanded") for i in range(3)}
    require(set(files) == expected_files and all(valid_hash(v) for v in files.values()), "Fresh eight-model inventory differs")
    require(len(seal["fingerprint_training"]) == 2
        and sum(v["n_structures"] for v in seal["fingerprint_training"]) == protocol["internal_counts"]["train"], "Fingerprint training identities omitted")
    for fold, lineage in enumerate(seal["fingerprint_training"]):
        require(lineage["fold"] == fold and lineage["dev_test_identities_used"] == 0
            and lineage["epochs_fixed_before_execution"] == 6 and lineage["model_sha256"] == files[f"external_fp_fold{fold}.pt"]
            and lineage["initialization"] == "fresh seeded random; no archived weights", "Fingerprint fitting crossed declared boundaries")
        require(lineage["n_spectra"] >= 100 and 0 < lineage["n_structures"] <= 300
            and lineage["n_spectra"] <= 128 * lineage["n_structures"]
            and len(lineage["training_bce_by_epoch"]) == 6
            and all(type(v) in (int, float) and math.isfinite(v) and v >= 0 for v in lineage["training_bce_by_epoch"]), "Invalid fixed-epoch fitting lineage")
    require(set(seal["models"]) == {"legacy", "expanded"}, "Missing fixed policy model")
    for training in seal["models"].values():
        require(training["columns"] and len(set(training["columns"])) == len(training["columns"])
            and not any(c.startswith(("prior_", "audit_")) for c in training["columns"]), "Membership shortcut entered scientific features")
        require(len(training["best_iterations"]) == 3 and all(type(v) is int and 1 <= v <= 400 for v in training["best_iterations"]), "Invalid ranker stopping record")
        config = training["config"]
        require(config["seeds"] == [0, 1, 2] and config["num_boost_round"] == 400 and config["early_stopping"] == 40
            and config["params"]["num_threads"] == 2, "Registered ranker fitting configuration differs")
        require(all(any(c.startswith(prefix) for c in training["columns"]) for prefix in ("d1_", "an_", "d3_", "d4_", "gen"))
            and "np_score" in training["columns"], "Declared model channel omitted")
    choice = seal["development_choice"]
    require(choice["n_molecules"] == 100 and choice["lost"] == 0
        and choice["expanded_covered"] - choice["legacy_covered"] == choice["gained"], "Development preservation accounting differs")
    require(all(type(choice[k]) is int and 0 <= choice[k] <= 100 for k in ("legacy_covered", "expanded_covered", "gained", "lost")), "Invalid development coverage counts")
    require(choice["selected_policy"] == ("expanded" if choice["gained"] >= 1 else "legacy"), "Policy choice contradicts registered development gate")
    require(set(result["models"]) == {"internal", "external"}, "Missing evaluation cohort")
    for cohort, count, classes, modes in (("internal", 150, {"C1", "C2", "C3"}, {"none"}), ("external", n, {"C2", "C3"}, {"none", "blank"})):
        summaries = result["models"][cohort]
        require(set(summaries) == {"legacy", "expanded", "mass_only", "d1", "d1_d2"}, "Missing policy or declared diagnostic baseline")
        for conditions in summaries.values():
            require(set(conditions) == modes, "Missing spectral condition")
            for cells in conditions.values():
                require(set(cells) == classes, "Missing novelty class")
                for cell in cells.values():
                    require(set(cell) == {"n", "mrr", "ci_low", "ci_high", "top1", "top10", "recall_at_25"} and cell["n"] == count, "Test denominator differs")
                    require(all(type(v) in (int, float) and math.isfinite(v) and 0 <= v <= 1 for k, v in cell.items() if k != "n"), "Invalid metric")
                    require(cell["top1"] <= cell["top10"] <= cell["recall_at_25"] and cell["top1"] <= cell["mrr"] <= cell["recall_at_25"]
                        and cell["ci_low"] <= cell["ci_high"], "Impossible metric summary")
        require(set(result["candidate_coverage"][cohort]) == {"legacy", "expanded"}, "Missing budget coverage")
        for policy, cells in result["candidate_coverage"][cohort].items():
            require(set(cells) == classes, "Missing coverage class")
            for cls, cell in cells.items():
                require(cell["n_molecules"] == count and type(cell["covered"]) is int and 0 <= cell["covered"] <= count
                    and close(cell["oracle_mrr_ceiling"], cell["covered"] / count), "Coverage denominator/ceiling differs")
                require(type(cell["maximum_candidates"]) is int and cell["maximum_candidates"] >= 0
                    and type(cell["mean_candidates"]) in (int, float) and math.isfinite(cell["mean_candidates"])
                    and 0 <= cell["mean_candidates"] <= cell["maximum_candidates"]
                    and cell["mean_candidates"] * count + 1e-9 >= cell["covered"], "Impossible candidate-list accounting")
                require(summaries[policy]["none"][cls]["recall_at_25"] <= cell["covered"] / count + 1e-12, "Rank recall exceeds candidate coverage")
        require(all(summaries[p][m]["C3"]["recall_at_25"] == 0 for p in ("mass_only", "d1", "d1_d2") for m in modes), "Truth remained in original C3 pool")
    primary = result["primary"]
    coverage = result["candidate_coverage"]["internal"]
    difference = coverage["expanded"]["C3"]["covered"] - coverage["legacy"]["C3"]["covered"]
    require(primary["n_molecules"] == 150 and primary["lost"] == 0 and primary["gained"] == difference
        and close(primary["delta_coverage"], difference / 150) and primary["n_bootstrap"] == 4000
        and -1 <= primary["ci_low"] <= primary["ci_high"] <= 1, "Primary paired coverage comparison differs")
    require(primary["legacy_covered"] == coverage["legacy"]["C3"]["covered"]
        and primary["expanded_covered"] == coverage["expanded"]["C3"]["covered"] and primary["seed"] == protocol["seed"], "Primary numerator or bootstrap seed differs")
    effects = {}
    for cohort in ("internal", "external"):
        for cls in result["candidate_coverage"][cohort]["legacy"]:
            model = result["models"][cohort]
            effects[f"{cohort}_mrr_expanded_minus_legacy:{cls}"] = model["expanded"]["none"][cls]["mrr"] - model["legacy"]["none"][cls]["mrr"]
            if cohort == "external":
                effects[f"external_coverage_expanded_minus_legacy:{cls}"] = (result["candidate_coverage"][cohort]["expanded"][cls]["covered"] - result["candidate_coverage"][cohort]["legacy"][cls]["covered"]) / n
                for policy in ("legacy", "expanded"):
                    effects[f"external_{policy}_blank_minus_none:{cls}"] = model[policy]["blank"][cls]["mrr"] - model[policy]["none"][cls]["mrr"]
    joint = result["secondary_joint_intervals"]
    require(set(effects) == set(execution["secondary_family"]) == set(joint["outcomes"]) and joint["family_size"] == 11
        and joint["cohort_counts"] == {"internal": 150, "external": n} and joint["n_bootstrap"] == 4000
        and joint["seed"] == protocol["seed"], "Joint secondary family differs")
    radius = joint["joint_radius"]
    require(type(radius) in (int, float) and math.isfinite(radius) and radius >= 0, "Invalid joint uncertainty radius")
    for name, effect in effects.items():
        cell = joint["outcomes"][name]
        require(cell["cohort"] == name.split("_")[0] and close(cell["delta"], effect)
            and close(cell["ci_low"], effect - radius) and close(cell["ci_high"], effect + radius), "Secondary effect contradicts absolute outcomes")
    manifest = result["source_manifest"]
    required = set(required_validation_inputs(pool_tag="_purged_external", analog=True)) | {"structures.parquet", "pool_source_quality_purged_external.json"}
    require(len(manifest["files"]) == len(required) and {e["path"] for e in manifest["files"]} == required
        and manifest["files_digest_sha256"] == _canonical_files_digest(manifest["files"]), "Static source inventory differs")
    require(all(type(e.get("bytes")) is int and e["bytes"] > 0 and valid_hash(e.get("sha256")) for e in manifest["files"] + list(raw.values())), "Invalid input byte/hash record")
    require(result["fingerprint_provenance"]["status"] == "POOL_FINGERPRINTS_VERIFIED"
        and result["analog_index_provenance"]["cache_used"] is False, "Unverified candidate fingerprints or reused D2 cache")
    require(result["execution"] == {"platform": "Kaggle", "competition_source": "enveda-CASMI26-molecule-id-mass-spectra",
        "dataset_sources": DATASETS, "pretrained_models_attached": False, "gpu_enabled": False, "internet_enabled": False,
        "submitted_predictions": False, "protocol_registered_commit": REGISTERED_COMMIT, "new_model_test_results_before_registration": False,
        "source_amendment_sha256": text_hash("protocols/purged-external-v1-source-amendment.json")}, "Execution policy differs")
    require(result["environment"]["device"] == "cpu" and result["environment"]["rdkit"] == "2026.03.3", "Unexpected scientific runtime")
    engine = result["engine_config"]
    require(all(engine[k] is True for k in ("generate", "gen_consensus", "gen_expand", "np_score"))
        and engine["gen_max_analogs"] == 30 and engine["gen_cap"] == 900, "Candidate-budget implementation differs")
    def reject_private(value):
        if isinstance(value, dict):
            require(not set(value) & {"qid", "query_ids", "test_ids", "smiles", "normalized_smiles", "ranks", "spectrum_rows", "selected_rows", "ms2_mzs", "accessions"}, "Private row-level export detected")
            for nested in value.values(): reject_private(nested)
        elif isinstance(value, list):
            for nested in value: reject_private(nested)
    reject_private(result)


def main():
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi(); api.authenticate()
    local = json.loads((ROOT / "kaggle/audits/purged-external/notebook.ipynb").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="casmi-purged-readback-") as temporary:
        folder = Path(temporary)
        api.kernels_pull(KERNEL, path=temporary, metadata=True, quiet=True)
        metadata = json.loads((folder / "kernel-metadata.json").read_text())
        remote = json.loads((folder / metadata["code_file"]).read_text())
        status = api.kernels_status(KERNEL)
        if not isinstance(status, dict): status = json.loads(str(status))
        require(status.get("status") == "COMPLETE" and source_hash(local) == source_hash(remote) and not metadata["is_private"]
            and sorted(metadata["dataset_sources"]) == sorted(DATASETS) and not metadata["enable_gpu"] and not metadata["enable_internet"]
            and not metadata["kernel_sources"] and metadata["competition_sources"] == ["enveda-CASMI26-molecule-id-mass-spectra"],
            "Public completed source-matched CPU experiment and attachments not verified")
        api.kernels_output(KERNEL, path=temporary, file_pattern=r"^purged_external_study\.json$", quiet=True)
        artifact = folder / "purged_external_study.json"
        require(artifact.stat().st_size <= 2_000_000, "Unbounded aggregate export")
        raw = artifact.read_bytes(); validate_result(json.loads(raw))
    (ROOT / "reports/purged_external_study.json").write_bytes(raw)
    record = {"checked_at": datetime.now(UTC).isoformat(), "id": KERNEL, "url": f"https://www.kaggle.com/code/{KERNEL}",
        "latest_run_status": "COMPLETE", "public": True, "source_matches": True, "source_sha256": source_hash(local),
        "aggregate_artifact_sha256": hashlib.sha256(raw).hexdigest(), "status": "pass",
        "scope": "Fresh models with selected external families purged from all training/reference data; original-snapshot novelty feasibility remains unsuccessful; no independent laboratory replication"}
    (ROOT / "reports/purged_external_publication.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__": main()
