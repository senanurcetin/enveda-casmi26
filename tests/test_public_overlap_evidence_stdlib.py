"""Bind committed real-data evidence to the executed source and readback artifact."""
import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path

from scripts.check_overlap_publication import validate_result
from scripts.check_kaggle_publications import source_hash

ROOT = Path(__file__).resolve().parents[1]


class PublishedOverlapEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.result = json.loads((ROOT / "reports/public_overlap_audit.json").read_text(encoding="utf-8"))

    def test_committed_evidence_matches_executed_source_and_artifact(self):
        validate_result(self.result)
        record = json.loads((ROOT / "reports/public_overlap_publication.json").read_text(encoding="utf-8"))
        canonical = (ROOT / "reports/public_overlap_audit.json").read_text(encoding="utf-8").encode()
        self.assertEqual(hashlib.sha256(canonical).hexdigest(), record["aggregate_artifact_sha256"])
        notebook = json.loads((ROOT / "kaggle/audits/public-overlap/notebook.ipynb").read_text(encoding="utf-8"))
        self.assertEqual(source_hash(notebook), record["source_sha256"])
        self.assertTrue(record["source_matches"] and record["public"])
        self.assertEqual(record["latest_run_status"], "COMPLETE")

    def test_altered_code_hash_impossible_counts_and_extra_fields_rejected(self):
        for change in ({"content_matched_spectra": self.result["public_spectra"] + 1},
                       {"public_spectra": True}, {"exact_matching_pairs": 0},
                       {"raw_molecule_ids": ["private"]}):
            result = deepcopy(self.result)
            result.update(change)
            with self.subTest(change=tuple(change)), self.assertRaises(ValueError):
                validate_result(result)
        result = deepcopy(self.result)
        result["audit_sources"]["scripts/audit_public_overlap.py"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "code hash"):
            validate_result(result)


if __name__ == "__main__":
    unittest.main()
