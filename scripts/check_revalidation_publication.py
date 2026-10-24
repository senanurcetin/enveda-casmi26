"""Verify the executed owned study before accepting aggregate research evidence."""
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
from scripts.build_revalidation_notebook import DATASETS, KERNEL, SCRIPTS
from scripts.check_kaggle_publications import source_hash
from casmi.validate.source_manifest import required_validation_inputs, _canonical_files_digest

COMPETITION = "enveda-CASMI26-molecule-id-mass-spectra"


def validate_result(result, root=ROOT):
    expected = {"schema_version", "status", "checked_at", "elapsed_seconds", "cohort",
                "source_artifact_manifest", "fingerprint_provenance", "analog_index_provenance",
                "knockout_audit", "model_configurations", "models", "paired_comparisons",
                "private_detail_sha256", "environment", "limitations", "raw_inputs",
                "candidate_source_quality", "queries_with_truth_absent_from_candidate_pool",
                "source_code_sha256", "execution"}
    if (set(result) != expected or result["schema_version"] != 1
            or result["status"] != "COMPLETE_RETROSPECTIVE_REVALIDATION"):
        raise ValueError("Unexpected or incomplete study report")
    sources = {path.relative_to(root).as_posix() for path in (root / "src").rglob("*.py")}
    sources.update("scripts/" + name for name in SCRIPTS)
    if set(result["source_code_sha256"]) != sources:
        raise ValueError("Executed source inventory differs")
    for name, digest in result["source_code_sha256"].items():
        if hashlib.sha256((root / name).read_text(encoding="utf-8").encode()).hexdigest() != digest:
            raise ValueError("Executed source hash differs: " + name)
    count = result["cohort"]["n_molecules"]
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("A nonempty molecule cohort is required")
    if result["cohort"]["library"] != "enveda-np-examples":
        raise ValueError("Unexpected evaluation population")
    manifest = result["source_artifact_manifest"]
    entries = manifest.get("files", [])
    required = set(required_validation_inputs(pool_tag="_corrected_study", analog=True))
    required.update({"structures.parquet", "pool_source_quality_corrected_study.json",
                     "models/fp_fold0.pt", "models/fp_fold1.pt"})
    required.update(f"models/ranker-v16-gen2_seed{seed}.txt" for seed in range(5))
    if (manifest.get("schema") != "casmi-local-artifact-manifest-v1"
            or len(entries) != len(required) or {entry.get("path") for entry in entries} != required
            or entries != sorted(entries, key=lambda entry: entry["path"])
            or manifest.get("files_digest_sha256") != _canonical_files_digest(entries)):
        raise ValueError("Incomplete or inconsistent static source/model inventory")
    for entry in entries:
        if type(entry.get("bytes")) is not int or entry["bytes"] <= 0 or not re.fullmatch(r"[0-9a-f]{64}", entry.get("sha256", "")):
            raise ValueError("Invalid static source/model hash")
    if result["fingerprint_provenance"].get("status") != "POOL_FINGERPRINTS_VERIFIED":
        raise ValueError("Unverified candidate fingerprints")
    analog = result["analog_index_provenance"]
    if analog.get("status") != "REBUILT_WITHOUT_CACHE" or analog.get("cache_used") is not False:
        raise ValueError("D2 index was not rebuilt without cache")
    quality = result["candidate_source_quality"]
    if (quality.get("mass_conflict_policy") != "quarantine_entire_identity"
            or quality.get("mass_conflict_groups") != 0
            or quality.get("quarantined_mass_conflict_groups", -1) < 0):
        raise ValueError("Candidate-source conflict policy is not verified")
    missing = result["queries_with_truth_absent_from_candidate_pool"]
    if type(missing) is not int or not 0 <= missing <= count:
        raise ValueError("Invalid missing-truth denominator accounting")
    if set(result["models"]) != {"d1", "d1_d2", "frozen_v16"}:
        raise ValueError("Incomplete model comparison")
    for block in result["models"].values():
        if set(block) != {"conditions", "controls"} or set(block["conditions"]) != {"none", "blank", "intensity_shuffle"}:
            raise ValueError("Incomplete model/control conditions")
        if set(block["controls"]) != {"blank", "intensity_shuffle"}:
            raise ValueError("Missing paired controls")
        for mode, control in block["controls"].items():
            if (control.get("status") != "DIAGNOSTIC_ONLY" or control.get("mode") != mode
                    or control.get("n_molecules") != count
                    or control.get("verified_original_source_pairing") is not True
                    or not 0 < control.get("changed_query_signatures", 0) <= count):
                raise ValueError("Unverified control pairing or spectral perturbation")
            for cls in ("C1", "C2", "C3"):
                values = control.get("classes", {}).get(cls, {})
                reference = block["conditions"]["none"][cls]["mrr"]
                treatment = block["conditions"][mode][cls]["mrr"]
                if (values.get("reference_mrr") != reference or values.get("control_mrr") != treatment
                        or not math.isclose(values.get("control_minus_reference", math.inf), treatment - reference, abs_tol=1e-12)):
                    raise ValueError("Paired control summaries contradict absolute results")
        for classes in block["conditions"].values():
            if set(classes) != {"C1", "C2", "C3"}:
                raise ValueError("Missing novelty class")
            for values in classes.values():
                if set(values) != {"mrr", "ci_low", "ci_high", "top1", "top10", "recall_at_25", "n"} or values["n"] != count:
                    raise ValueError("Inconsistent denominators or unexpected details")
                for field in ("mrr", "ci_low", "ci_high", "top1", "top10", "recall_at_25"):
                    if not isinstance(values[field], (int, float)) or not math.isfinite(values[field]) or not 0 <= values[field] <= 1:
                        raise ValueError("Invalid summary metric")
                if not values["ci_low"] <= values["mrr"] <= values["ci_high"]:
                    raise ValueError("Invalid confidence interval")
                if not values["top1"] <= values["top10"] <= values["recall_at_25"]:
                    raise ValueError("Invalid rank-hit accounting")
    if set(result["raw_inputs"]) != {"train.parquet", "coconut.csv", "pubchemlite.csv"}:
        raise ValueError("Incomplete raw input hashes")
    for record in result["raw_inputs"].values():
        if set(record) != {"bytes", "sha256"} or record["bytes"] <= 0 or not re.fullmatch(r"[0-9a-f]{64}", record["sha256"]):
            raise ValueError("Invalid input identity")
    execution = result["execution"]
    if (execution != {"platform": "Kaggle", "competition_source": COMPETITION,
                      "dataset_sources": DATASETS, "internet_enabled": False, "submitted_predictions": False}):
        raise ValueError("Unexpected execution inputs or submission state")
    # These row-level fields must never appear anywhere in an exported report.
    def reject_details(value):
        if isinstance(value, dict):
            if set(value) & {"query_ids", "query_signatures", "source_query_signatures", "ranks", "smiles", "spectrum_rows"}:
                raise ValueError("Unexpected private row-level details")
            for nested in value.values():
                reject_details(nested)
        elif isinstance(value, list):
            for nested in value:
                reject_details(nested)
    reject_details(result)


def main():
    from kaggle.api.kaggle_api_extended import KaggleApi

    api = KaggleApi()
    api.authenticate()
    local = json.loads((ROOT / "kaggle/audits/corrected-revalidation/notebook.ipynb").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="casmi-revalidation-") as temporary:
        folder = Path(temporary)
        api.kernels_pull(KERNEL, path=temporary, metadata=True, quiet=True)
        metadata = json.loads((folder / "kernel-metadata.json").read_text(encoding="utf-8"))
        remote = json.loads((folder / metadata["code_file"]).read_text(encoding="utf-8"))
        status = api.kernels_status(KERNEL)
        if not isinstance(status, dict):
            status = json.loads(str(status))
        if (status.get("status") != "COMPLETE" or source_hash(local) != source_hash(remote)
                or metadata["is_private"] or metadata["competition_sources"] != [COMPETITION]
                or sorted(metadata["dataset_sources"]) != sorted(DATASETS) or metadata["kernel_sources"]
                or metadata["enable_internet"] or not metadata["enable_gpu"]):
            raise ValueError("The public offline study and executed source are not verified")
        api.kernels_output(KERNEL, path=temporary, file_pattern=r"^revalidation_study\.json$", quiet=True)
        artifact = folder / "revalidation_study.json"
        if artifact.stat().st_size > 500_000:
            raise ValueError("Aggregate artifact exceeds bounded size")
        raw = artifact.read_bytes()
        validate_result(json.loads(raw))
    (ROOT / "reports/revalidation_study.json").write_bytes(raw)
    record = {"checked_at": datetime.now(UTC).isoformat(), "id": KERNEL,
              "url": f"https://www.kaggle.com/code/{KERNEL}", "latest_run_status": "COMPLETE",
              "public": True, "source_matches": True, "source_sha256": source_hash(local),
              "aggregate_artifact_sha256": hashlib.sha256(raw).hexdigest(), "status": "pass",
              "scope": "Retrospective frozen-model diagnostics; not a new untouched generalization score"}
    (ROOT / "reports/revalidation_publication.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
