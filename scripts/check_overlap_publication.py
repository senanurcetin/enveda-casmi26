"""Opt-in readback of the owned Kaggle content audit; download only its aggregate JSON."""
from __future__ import annotations

import hashlib
import json
import re
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.check_kaggle_publications import source_hash  # noqa: E402
KERNEL = "senanuretin/casmi26-exact-public-example-overlap-audit"
COMPETITION = "enveda-CASMI26-molecule-id-mass-spectra"
SOURCE_FILES = {"src/casmi/validate/spectrum_overlap.py", "src/casmi/validate/overlap_audit.py", "scripts/audit_public_overlap.py"}
COUNTS = {"train_rows_scanned", "public_spectra", "public_molecules", "exact_precursor_candidate_rows",
          "coarse_candidate_pairs_after_precursor_filter", "coarse_false_match_pairs", "exact_matching_pairs",
          "content_matched_spectra", "metadata_consistent_matched_spectra", "metadata_disagreement_pairs",
          "fully_content_matched_molecules", "fully_metadata_consistent_molecules"}


def validate_result(result, root=ROOT):
    permitted = COUNTS | {"schema_version", "status", "matching_method", "metadata_fields", "content_matches_by_library",
                          "limitations", "checked_at", "inputs", "audit_sources", "environment", "batch_size", "execution"}
    if set(result) != permitted or result["schema_version"] != 1 or result["status"] != "COMPLETE_CONTENT_AUDIT":
        raise ValueError("Unexpected audit result schema or incomplete scan")
    for name in COUNTS | {"batch_size"}:
        value = result[name]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError("Audit counts must be nonnegative integers")
    if not result["train_rows_scanned"] or not result["public_spectra"] or not result["batch_size"]:
        raise ValueError("Audit input populations and batch size must be positive")
    if not (0 <= result["metadata_consistent_matched_spectra"] <= result["content_matched_spectra"] <= result["public_spectra"]):
        raise ValueError("Impossible matched-spectrum counts")
    if not (result["fully_metadata_consistent_molecules"] <= result["fully_content_matched_molecules"] <= result["public_molecules"] <= result["public_spectra"]):
        raise ValueError("Impossible matched-molecule counts")
    if (result["coarse_false_match_pairs"] + result["exact_matching_pairs"] != result["coarse_candidate_pairs_after_precursor_filter"]
            or result["exact_matching_pairs"] < result["content_matched_spectra"]
            or result["metadata_disagreement_pairs"] > result["exact_matching_pairs"]
            or result["exact_precursor_candidate_rows"] > result["train_rows_scanned"]):
        raise ValueError("Inconsistent candidate/match accounting")
    for name, count in result["content_matches_by_library"].items():
        if (not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_. -]{1,80}", name)
                or isinstance(count, bool) or not isinstance(count, int) or not 0 <= count <= result["public_spectra"]):
            raise ValueError("Invalid aggregate source-library count")
    if set(result["inputs"]) != {"train.parquet", "test.parquet"}:
        raise ValueError("Unexpected input inventory")
    for entry in result["inputs"].values():
        if (set(entry) != {"bytes", "sha256"} or isinstance(entry["bytes"], bool)
                or not isinstance(entry["bytes"], int) or entry["bytes"] <= 0
                or not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"])):
            raise ValueError("Invalid input hash inventory")
    if set(result["audit_sources"]) != SOURCE_FILES:
        raise ValueError("Incomplete audit source inventory")
    for name, expected in result["audit_sources"].items():
        actual = hashlib.sha256((root / name).read_text(encoding="utf-8").encode()).hexdigest()
        if actual != expected:
            raise ValueError("Audit code hash differs from the executed source")
    if result["metadata_fields"] != ["adduct", "ionization_mode"]:
        raise ValueError("Unexpected acquisition-metadata comparison")
    execution = result["execution"]
    if execution.get("platform") != "Kaggle" or execution.get("competition_source") != COMPETITION:
        raise ValueError("Unexpected execution source")


def main():
    # Local invocation works as both python -m and python scripts/... .
    from kaggle.api.kaggle_api_extended import KaggleApi

    api = KaggleApi()
    api.authenticate()
    local_dir = ROOT / "kaggle/audits/public-overlap"
    local = json.loads((local_dir / "notebook.ipynb").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="casmi-public-overlap-") as temporary:
        folder = Path(temporary)
        api.kernels_pull(KERNEL, path=temporary, metadata=True, quiet=True)
        metadata = json.loads((folder / "kernel-metadata.json").read_text(encoding="utf-8"))
        remote = json.loads((folder / metadata["code_file"]).read_text(encoding="utf-8"))
        status = api.kernels_status(KERNEL)
        if not isinstance(status, dict):
            status = json.loads(str(status))
        if (status.get("status") != "COMPLETE" or source_hash(local) != source_hash(remote)
                or metadata["is_private"] or metadata["competition_sources"] != [COMPETITION]
                or metadata["dataset_sources"] or metadata["kernel_sources"]
                or metadata["enable_internet"] or metadata["enable_gpu"]):
            raise ValueError("Public CPU-only audit run, source and inputs are not yet verified")
        api.kernels_output(KERNEL, path=temporary, file_pattern=r"^public_overlap_audit\.json$", quiet=True)
        artifact = folder / "public_overlap_audit.json"
        if artifact.stat().st_size > 100_000:
            raise ValueError("Aggregate audit artifact exceeds expected bounded size")
        raw = artifact.read_bytes()
        result = json.loads(raw)
        validate_result(result)
    (ROOT / "reports/public_overlap_audit.json").write_bytes(raw)
    record = {"checked_at": datetime.now(UTC).isoformat(), "id": KERNEL,
              "url": f"https://www.kaggle.com/code/{KERNEL}", "latest_run_status": "COMPLETE",
              "public": True, "source_matches": True, "source_sha256": source_hash(local),
              "aggregate_artifact_sha256": hashlib.sha256(raw).hexdigest(), "status": "pass",
              "scope": "Executed public-example content audit; no hidden-test or model-performance claim"}
    (ROOT / "reports/public_overlap_publication.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
