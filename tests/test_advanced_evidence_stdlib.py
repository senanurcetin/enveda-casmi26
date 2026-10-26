"""Authenticate fresh-study evidence and reject scientifically inconsistent exports."""
import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path

from scripts.check_advanced_publication import digest, validate_result
from scripts.check_kaggle_publications import source_hash

ROOT = Path(__file__).resolve().parents[1]


class AdvancedEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads((ROOT / "reports/advanced_study.json").read_text(encoding="utf-8"))

    def test_completed_artifact_source_and_raw_input_binding(self):
        validate_result(self.result)
        publication = json.loads((ROOT / "reports/advanced_publication.json").read_text(encoding="utf-8"))
        raw = (ROOT / "reports/advanced_study.json").read_text(encoding="utf-8").encode("utf-8")
        self.assertEqual(hashlib.sha256(raw).hexdigest(), publication["aggregate_artifact_sha256"])
        notebook = json.loads((ROOT / "kaggle/audits/advanced-internal/notebook.ipynb").read_text(encoding="utf-8"))
        self.assertEqual(source_hash(notebook), publication["source_sha256"])
        self.assertEqual(publication["latest_run_status"], "COMPLETE")
        self.assertTrue(publication["public"] and publication["source_matches"])
        historical = json.loads((ROOT / "reports/revalidation_study.json").read_text(encoding="utf-8"))
        self.assertEqual(self.result["raw_inputs"], historical["raw_inputs"])

    def test_rejects_source_partition_model_and_private_row_tampering(self):
        mutations = [
            lambda r: r.update(status="RUNNING"),
            lambda r: r["source_code_sha256"].update({"scripts/run_advanced_study.py": "0" * 64}),
            lambda r: r["partitions"].update(cross_split_identity_overlaps=1),
            lambda r: r["partitions"]["test"].update(n_molecules=249),
            lambda r: r["training_seal"]["model_files"].update({"advanced_fp_fold0.pt": "0" * 64}),
            lambda r: r["models"]["full"]["none"]["C2"].update(n=249),
            lambda r: r["execution"].update(pretrained_models_attached=True),
            lambda r: r["execution"].update(source_amendment_sha256="0" * 64),
            lambda r: r["environment"].update(query_ids=["private"]),
            lambda r: r["raw_inputs"]["train.parquet"].update(sha256="invalid"),
        ]
        for index, mutate in enumerate(mutations):
            changed = deepcopy(self.result)
            mutate(changed)
            with self.subTest(case=index), self.assertRaises(ValueError):
                validate_result(changed)

    def test_rejects_resealed_untrained_or_shortcut_scientific_models(self):
        for field in ("dev_test_identities_used", "columns"):
            changed = deepcopy(self.result)
            if field == "columns":
                changed["training_seal"]["models"]["full"]["columns"].append("prior_coconut")
            else:
                changed["training_seal"]["fingerprint_training"][0][field] = 1
            changed["training_seal_sha256"] = digest(changed["training_seal"])
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_result(changed)

    def test_rejects_contradictory_effects_coverage_and_uncertainty(self):
        mutations = [
            lambda r: r["primary"].update(delta_mrr=0),
            lambda r: r["primary"].update(ci_low=2),
            lambda r: r["secondary_joint_intervals"]["outcomes"]["no_d1_minus_full:C2"].update(delta_mrr=0),
            lambda r: r["secondary_joint_intervals"].update(joint_radius=0),
            lambda r: r["candidate_coverage"]["C3"].update(pool_truth_present=1),
            lambda r: r["candidate_coverage"]["C2"].update(pool_truth_present=0, with_generation_truth_present=0),
            lambda r: r["models"]["full"]["none"]["C2"].update(top1=1),
        ]
        for index, mutate in enumerate(mutations):
            changed = deepcopy(self.result)
            mutate(changed)
            with self.subTest(case=index), self.assertRaises(ValueError):
                validate_result(changed)


if __name__ == "__main__":
    unittest.main()
