"""Reject incomplete, stale-source and row-level research exports without datasets."""
import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path

from scripts.build_revalidation_notebook import DATASETS, SCRIPTS
from scripts.check_revalidation_publication import COMPETITION, validate_result
from casmi.validate.source_manifest import required_validation_inputs, _canonical_files_digest

ROOT = Path(__file__).resolve().parents[1]


def synthetic_report():
    paths = sorted((ROOT / "src").rglob("*.py")) + [ROOT / "scripts" / name for name in SCRIPTS]
    summary = {"mrr": .2, "ci_low": .1, "ci_high": .3,
               "top1": .1, "top10": .3, "recall_at_25": .4, "n": 10}
    result = dict.fromkeys(["checked_at", "elapsed_seconds", "source_artifact_manifest",
                           "fingerprint_provenance", "analog_index_provenance", "knockout_audit",
                           "model_configurations", "paired_comparisons", "private_detail_sha256",
                           "environment", "limitations", "candidate_source_quality",
                           "queries_with_truth_absent_from_candidate_pool"])
    result.update({"schema_version": 1, "status": "COMPLETE_RETROSPECTIVE_REVALIDATION",
                   "cohort": {"n_molecules": 10, "library": "enveda-np-examples"},
                   "models": {model: {"conditions": {mode: {cls: deepcopy(summary) for cls in ("C1", "C2", "C3")}
                                                       for mode in ("none", "blank", "intensity_shuffle")},
                                      "controls": {mode: {"status": "DIAGNOSTIC_ONLY", "mode": mode,
                                                           "n_molecules": 10, "verified_original_source_pairing": True,
                                                           "changed_query_signatures": 10,
                                                           "classes": {cls: {"reference_mrr": .2, "control_mrr": .2,
                                                                             "control_minus_reference": 0.0} for cls in ("C1", "C2", "C3")}}
                                                   for mode in ("blank", "intensity_shuffle")}}
                              for model in ("d1", "d1_d2", "frozen_v16")},
                   "source_code_sha256": {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_text(encoding="utf-8").encode()).hexdigest() for path in paths},
                   "raw_inputs": {name: {"bytes": 1, "sha256": "a" * 64} for name in ("train.parquet", "coconut.csv", "pubchemlite.csv")},
                   "execution": {"platform": "Kaggle", "competition_source": COMPETITION,
                                 "dataset_sources": DATASETS, "internet_enabled": False, "submitted_predictions": False}})
    result.update({"fingerprint_provenance": {"status": "POOL_FINGERPRINTS_VERIFIED"},
                   "analog_index_provenance": {"status": "REBUILT_WITHOUT_CACHE", "cache_used": False},
                   "candidate_source_quality": {"mass_conflict_policy": "quarantine_entire_identity",
                                                "mass_conflict_groups": 0, "quarantined_mass_conflict_groups": 2},
                   "queries_with_truth_absent_from_candidate_pool": 1})
    names = list(required_validation_inputs(pool_tag="_corrected_study", analog=True))
    names += ["structures.parquet", "pool_source_quality_corrected_study.json", "models/fp_fold0.pt", "models/fp_fold1.pt"]
    names += [f"models/ranker-v16-gen2_seed{seed}.txt" for seed in range(5)]
    entries = [{"path": name, "bytes": 1, "sha256": "a" * 64} for name in sorted(names)]
    result["source_artifact_manifest"] = {"schema": "casmi-local-artifact-manifest-v1", "files": entries,
                                         "files_digest_sha256": _canonical_files_digest(entries)}
    return result


class StudyReadbackTests(unittest.TestCase):
    def test_rejects_incomplete_stale_and_private_details(self):
        result = synthetic_report()
        validate_result(result)
        for field, value in (("status", "RUNNING"), ("query_ids", ["private"])):
            changed = deepcopy(result)
            changed[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_result(changed)
        changed = deepcopy(result)
        changed["source_code_sha256"]["scripts/run_revalidation_study.py"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "source hash"):
            validate_result(changed)
        changed = deepcopy(result)
        changed["models"]["d1"]["controls"]["blank"]["ranks"] = [1]
        with self.assertRaisesRegex(ValueError, "private row-level"):
            validate_result(changed)

    def test_rejects_changed_denominators_and_invalid_rank_metrics(self):
        for field, value in (("n", 9), ("mrr", float("nan")), ("top1", .8)):
            result = synthetic_report()
            result["models"]["d1"]["conditions"]["none"]["C1"][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_result(result)


if __name__ == "__main__":
    unittest.main()
