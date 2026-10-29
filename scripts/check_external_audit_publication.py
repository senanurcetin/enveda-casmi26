"""Authenticate external-source feasibility separately from model performance."""
from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts.build_external_audit_notebook import KERNEL, DATASETS
from scripts.check_kaggle_publications import source_hash
from scripts.research_snapshot import snapshot_hashes


def validate_result(result, *, protocol_path="protocols/external-c3-v1.json",
                    prepared_path="reports/external_source_preparation.json",
                    notebook_path="kaggle/audits/external-independence/notebook.ipynb",
                    external_filename="massbank_cc0.parquet", registered_commit="6164c5e"):
    protocol_text = (ROOT / protocol_path).read_text(encoding="utf-8")
    protocol = json.loads(protocol_text)
    prepared = json.loads((ROOT / prepared_path).read_text(encoding="utf-8"))
    previous = json.loads((ROOT / "reports/public_overlap_audit.json").read_text(encoding="utf-8"))
    expected = {"schema_version", "status", "checked_at", "protocol_sha256", "inputs", "external_source", "identity_audit",
                "stored_content_audit", "model_performance_inspected", "limitations", "source_code_sha256", "execution"}
    if set(result) != expected or result["schema_version"] != 1 or result["model_performance_inspected"] is not False:
        raise ValueError("Unexpected external feasibility export or hidden model scoring")
    if (result["protocol_sha256"] != hashlib.sha256(protocol_text.encode()).hexdigest()
            or result["source_code_sha256"] != snapshot_hashes(ROOT / notebook_path)
            or result["external_source"] != prepared or prepared["source"] != protocol["external_source"]
            or result["inputs"]["train.parquet"] != previous["inputs"]["train.parquet"]
            or result["inputs"][external_filename]["sha256"] != prepared["parquet_sha256"]
            or result["inputs"][external_filename]["bytes"] != prepared["parquet_bytes"]):
        raise ValueError("Source, protocol or exact input provenance differs")
    audit, content = result["identity_audit"], result["stored_content_audit"]
    for key, value in audit.items():
        if not key.endswith("sha256") and (type(value) is not int or value < 0):
            raise ValueError("Invalid external audit counter")
    if (audit["training_rows_scanned"] != previous["train_rows_scanned"]
            or content["train_rows_scanned"] != audit["training_rows_scanned"]
            or content["public_spectra"] != prepared["counts"]["eligible_spectra"]
            or content["public_molecules"] != prepared["counts"]["eligible_metric_molecules"]):
        raise ValueError("Full-source denominator coverage differs")
    overlaps = [audit[key] for key in ("external_standard_overlap_molecules", "external_metric_overlap_molecules", "external_exact_content_overlap_molecules")]
    if (not max(overlaps) <= audit["excluded_external_molecules_union"] <= sum(overlaps)
            or audit["excluded_external_molecules_union"] + audit["novel_external_molecules"] != content["public_molecules"]
            or not audit["novel_external_molecules"] <= audit["novel_external_spectra"] <= content["public_spectra"]):
        raise ValueError("Impossible external overlap/novelty accounting")
    feasible = audit["unparseable_training_labels"] == 0 and audit["novel_external_molecules"] >= 30
    expected_status = "EXTERNAL_INDEPENDENCE_FEASIBLE" if feasible else "EXTERNAL_INDEPENDENCE_NOT_ESTABLISHED"
    if result["status"] != expected_status:
        raise ValueError("Feasibility conclusion contradicts observed counters")
    if result["execution"] != {"platform": "Kaggle", "internet_enabled": False, "gpu_enabled": False,
                               "model_training": False, "submitted_predictions": False, "protocol_registered_commit": registered_commit}:
        raise ValueError("Unexpected execution, fitting or submission policy")
    def private(value):
        if isinstance(value, dict):
            if set(value) & {"smiles", "normalized_smiles", "query_ids", "ranks", "accessions", "ms2_mzs", "ms2_normalized_intensities"}:
                raise ValueError("Row-level data escaped aggregate export")
            for nested in value.values(): private(nested)
        elif isinstance(value, list):
            for nested in value: private(nested)
    private(result)


def authenticate(*, kernel=KERNEL, datasets=DATASETS, output_prefix="external_audit", **validation_options):
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi(); api.authenticate()
    local = json.loads((ROOT / validation_options.get("notebook_path", "kaggle/audits/external-independence/notebook.ipynb")).read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="casmi-external-audit-") as temporary:
        folder = Path(temporary)
        api.kernels_pull(kernel, path=temporary, metadata=True, quiet=True)
        metadata = json.loads((folder / "kernel-metadata.json").read_text(encoding="utf-8"))
        remote = json.loads((folder / metadata["code_file"]).read_text(encoding="utf-8"))
        status = api.kernels_status(kernel)
        if not isinstance(status, dict): status = json.loads(str(status))
        if (status.get("status") != "COMPLETE" or source_hash(local) != source_hash(remote) or metadata["is_private"]
                or sorted(metadata["dataset_sources"]) != sorted(datasets) or metadata["enable_internet"]
                or metadata["enable_gpu"] or metadata["kernel_sources"]
                or metadata["competition_sources"] != ["enveda-CASMI26-molecule-id-mass-spectra"]):
            raise ValueError("Public source-matched completed audit and input attachments not verified")
        api.kernels_output(kernel, path=temporary, file_pattern=r"^external_audit\.json$", quiet=True)
        path = folder / "external_audit.json"
        if path.stat().st_size > 250_000: raise ValueError("Unbounded aggregate export")
        raw = path.read_bytes(); result = json.loads(raw); validate_result(result, **validation_options)
    (ROOT / f"reports/{output_prefix}.json").write_bytes(raw)
    record = {"checked_at": datetime.now(UTC).isoformat(), "id": kernel, "public": True,
              "url": f"https://www.kaggle.com/code/{kernel}", "latest_run_status": "COMPLETE", "source_matches": True,
              "source_sha256": source_hash(local), "aggregate_artifact_sha256": hashlib.sha256(raw).hexdigest(),
              "status": "pass", "scientific_status": result["status"], "scope": "External source independence feasibility only; no model outcome evidence"}
    (ROOT / f"reports/{output_prefix}_publication.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__": authenticate()
