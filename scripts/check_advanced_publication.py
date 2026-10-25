"""Authenticate the completed fresh-weight internal study and bounded aggregate export."""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts.build_advanced_notebook import DATASETS, KERNEL
from scripts.check_kaggle_publications import source_hash
from scripts.research_snapshot import snapshot_hashes
from casmi.validate.source_manifest import required_validation_inputs, _canonical_files_digest

COMPETITION = "enveda-CASMI26-molecule-id-mass-spectra"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def validate_result(result, root=ROOT):
    expected = {"schema_version", "status", "checked_at", "elapsed_seconds", "protocol", "protocol_sha256",
                "partitions", "training_seal", "training_seal_sha256", "source_manifest", "candidate_source_quality",
                "fingerprint_provenance", "analog_index_provenance", "search_memo", "engine_config", "models",
                "candidate_coverage", "primary", "secondary_joint_intervals", "error_slices", "private_details_sha256",
                "environment", "limitations", "raw_inputs", "source_code_sha256", "execution", "installed_packages"}
    if set(result) != expected or result["schema_version"] != 1 or result["status"] != "COMPLETE_INTERNAL_FAMILY_HOLDOUT":
        raise ValueError("Unexpected or incomplete fresh-study schema")
    text = (root / "protocols/advanced-internal-v1.json").read_text(encoding="utf-8")
    protocol = json.loads(text)
    if result["protocol"] != protocol or result["protocol_sha256"] != hashlib.sha256(text.encode()).hexdigest():
        raise ValueError("Executed protocol differs from the predeclared protocol")
    if result["source_code_sha256"] != snapshot_hashes(root / "kaggle/audits/advanced-internal/notebook.ipynb"):
        raise ValueError("Executed source snapshot differs")
    partitions = result["partitions"]
    for name, count in protocol["split_counts"].items():
        block = partitions[name]
        if block["n_molecules"] != count or sum(block["libraries"].values()) != count:
            raise ValueError("Partition denominators differ from the declared population")
        for field in ("identity_sha256", "family_sha256", "metric_sha256"):
            if not re.fullmatch(r"[0-9a-f]{64}", block[field]):
                raise ValueError("Invalid private partition digest")
    if partitions["cross_split_identity_overlaps"] != 0:
        raise ValueError("Supervised split leakage detected")
    seal = result["training_seal"]
    if digest(seal) != result["training_seal_sha256"] or seal["status"] != "ALL_MODELS_SEALED_BEFORE_TEST_FEATURES":
        raise ValueError("Model seal identity differs")
    if (seal["test_identities_sha256"] != partitions["test"]["identity_sha256"]
            or seal["protocol_sha256"] != result["protocol_sha256"]
            or datetime.fromisoformat(seal["sealed_at"]) > datetime.fromisoformat(result["checked_at"])):
        raise ValueError("Test population or temporal model seal differs")
    variants = set(protocol["ablations"])
    if set(seal["models"]) != variants or len(seal["fingerprint_training"]) != 2:
        raise ValueError("Incomplete fresh training inventory")
    files = seal["model_files"]
    expected_files = {f"advanced_fp_fold{fold}.pt" for fold in (0, 1)}
    expected_files.update(f"advanced_{name}_{seed}.txt" for name in variants for seed in range(3))
    if set(files) != expected_files or any(not re.fullmatch(r"[0-9a-f]{64}", value) for value in files.values()):
        raise ValueError("Unexpected learned artifact inventory")
    for fold, lineage in enumerate(seal["fingerprint_training"]):
        if (lineage["fold"] != fold or lineage["dev_test_identities_used"] != 0
                or lineage["epochs_fixed_before_execution"] != protocol["fp_model"]["epochs"]
                or lineage["model_sha256"] != files[f"advanced_fp_fold{fold}.pt"]
                or not 0 < lineage["n_structures"] <= protocol["split_counts"]["train"]
                or lineage["n_spectra"] < 100
                or len(lineage["training_bce_by_epoch"]) != protocol["fp_model"]["epochs"]):
            raise ValueError("Fresh fingerprint fitting lineage differs")
        if any(not math.isfinite(value) or value < 0 for value in lineage["training_bce_by_epoch"]):
            raise ValueError("Invalid fingerprint training loss history")
    for name, training in seal["models"].items():
        if training["model_sha256"] != [files[f"advanced_{name}_{i}.txt"] for i in range(3)]:
            raise ValueError("Ranker weights differ from their sealed identities")
        if name != "membership_only" and any(column.startswith("prior_") for column in training["columns"]):
            raise ValueError("Undeclared source-membership features entered the scientific ranker")
        if len(training["best_iterations"]) != 3 or any(not 1 <= value <= protocol["ranker"]["rounds"] for value in training["best_iterations"]):
            raise ValueError("Invalid development stopping record")
    count = protocol["split_counts"]["test"]
    if set(result["models"]) != variants | {"mass_only", "d1", "d1_d2"}:
        raise ValueError("Incomplete model/ablation comparison")
    for conditions in result["models"].values():
        if set(conditions) != set(protocol["controls"]):
            raise ValueError("Missing spectral control")
        for classes in conditions.values():
            if set(classes) != {"C1", "C2", "C3"}:
                raise ValueError("Missing novelty condition")
            for values in classes.values():
                if set(values) != {"mrr", "ci_low", "ci_high", "top1", "top10", "recall_at_25", "n"} or values["n"] != count:
                    raise ValueError("Test denominators or aggregate fields differ")
                for field in ("mrr", "ci_low", "ci_high", "top1", "top10", "recall_at_25"):
                    if type(values[field]) not in (int, float) or not math.isfinite(values[field]) or not 0 <= values[field] <= 1:
                        raise ValueError("Invalid test metric")
                if not values["top1"] <= values["top10"] <= values["recall_at_25"] or values["ci_low"] > values["ci_high"]:
                    raise ValueError("Impossible rank summary")
    delta = result["models"]["full"]["none"]["C2"]["mrr"] - result["models"]["d1_d2"]["none"]["C2"]["mrr"]
    if not math.isclose(delta, result["primary"]["delta_mrr"], abs_tol=1e-12):
        raise ValueError("Primary contrast contradicts absolute model summaries")
    joint = result["secondary_joint_intervals"]
    if joint["family_size"] != 27 or len(joint["outcomes"]) != 27 or joint["n_molecules"] != count or joint["n_bootstrap"] != protocol["bootstrap_draws"]:
        raise ValueError("Secondary joint comparison family changed")
    for cls, coverage in result["candidate_coverage"].items():
        if not 0 <= coverage["pool_truth_present"] <= coverage["with_generation_truth_present"] <= count:
            raise ValueError("Impossible candidate coverage accounting")
    if result["candidate_coverage"]["C3"]["pool_truth_present"] != 0:
        raise ValueError("Truth remained in original C3 candidate pool")
    for blocks in result["error_slices"].values():
        if sum(block["n"] for block in blocks.values()) != count:
            raise ValueError("Error slice denominator differs")
        for block in blocks.values():
            if block["n"] < protocol["slice_min_n"] and (block["status"] != "SUPPRESSED_SMALL_CELL" or "classes" in block):
                raise ValueError("Small error slice was not suppressed")
    manifest = result["source_manifest"]
    required = set(required_validation_inputs(pool_tag="_advanced_internal", analog=True))
    required.update({"structures.parquet", "pool_source_quality_advanced_internal.json"})
    entries = manifest["files"]
    if (len(entries) != len(required) or {entry["path"] for entry in entries} != required
            or manifest["files_digest_sha256"] != _canonical_files_digest(entries)):
        raise ValueError("Static input inventory is incomplete or altered")
    if result["fingerprint_provenance"]["status"] != "POOL_FINGERPRINTS_VERIFIED" or result["analog_index_provenance"]["cache_used"] is not False:
        raise ValueError("Unverified fingerprints or stale reference index")
    execution = result["execution"]
    if (execution["platform"] != "Kaggle" or execution["competition_source"] != COMPETITION
            or execution["dataset_sources"] != DATASETS or execution["pretrained_models_attached"] is not False
            or execution["internet_enabled"] is not False or execution["submitted_predictions"] is not False
            or execution["protocol_preregistered_commit"] != "4f8ed74"):
        raise ValueError("Unexpected data, model attachment or submission execution")
    def reject_private(value):
        if isinstance(value, dict):
            if set(value) & {"qid", "query_ids", "test_ids", "smiles", "normalized_smiles", "ranks", "spectrum_rows", "selected_rows"}:
                raise ValueError("Private row-level evidence escaped the aggregate export")
            for nested in value.values():
                reject_private(nested)
        elif isinstance(value, list):
            for nested in value:
                reject_private(nested)
    reject_private(result)


def main():
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi()
    api.authenticate()
    local = json.loads((ROOT / "kaggle/audits/advanced-internal/notebook.ipynb").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="casmi-advanced-") as temporary:
        folder = Path(temporary)
        api.kernels_pull(KERNEL, path=temporary, metadata=True, quiet=True)
        metadata = json.loads((folder / "kernel-metadata.json").read_text(encoding="utf-8"))
        remote = json.loads((folder / metadata["code_file"]).read_text(encoding="utf-8"))
        status = api.kernels_status(KERNEL)
        if not isinstance(status, dict):
            status = json.loads(str(status))
        if (status.get("status") != "COMPLETE" or source_hash(local) != source_hash(remote) or metadata["is_private"]
                or metadata["competition_sources"] != [COMPETITION] or sorted(metadata["dataset_sources"]) != sorted(DATASETS)
                or metadata["kernel_sources"] or metadata["enable_internet"] or not metadata["enable_gpu"]):
            raise ValueError("Completed public execution, exact source and fresh-weight inputs are not verified")
        api.kernels_output(KERNEL, path=temporary, file_pattern=r"^advanced_study\.json$", quiet=True)
        artifact = folder / "advanced_study.json"
        if artifact.stat().st_size > 2_000_000:
            raise ValueError("Aggregate study artifact exceeds bounded size")
        raw = artifact.read_bytes()
        validate_result(json.loads(raw))
    (ROOT / "reports/advanced_study.json").write_bytes(raw)
    record = {"checked_at": datetime.now(UTC).isoformat(), "id": KERNEL, "url": f"https://www.kaggle.com/code/{KERNEL}",
              "latest_run_status": "COMPLETE", "public": True, "source_matches": True, "source_sha256": source_hash(local),
              "aggregate_artifact_sha256": hashlib.sha256(raw).hexdigest(), "status": "pass",
              "scope": "Fresh supervised weights and sealed internal chemical-family holdout; not external validation or competition scoring"}
    (ROOT / "reports/advanced_publication.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
