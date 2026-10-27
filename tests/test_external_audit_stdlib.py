"""Protect negative external evidence against accidental promotion or tampering."""
import copy
import json
import unittest
from pathlib import Path

from scripts.check_external_audit_publication import validate_result

ROOT = Path(__file__).resolve().parents[1]


class ExternalAuditTests(unittest.TestCase):
    def setUp(self):
        self.report = json.loads((ROOT / "reports/external_audit.json").read_text(encoding="utf-8"))

    def test_authenticated_negative_evidence_is_valid(self):
        validate_result(self.report)
        self.assertEqual(self.report["status"], "EXTERNAL_INDEPENDENCE_NOT_ESTABLISHED")
        self.assertEqual(self.report["identity_audit"]["novel_external_molecules"], 0)

    def test_tampered_evidence_is_rejected(self):
        alterations = [
            ("status", "EXTERNAL_INDEPENDENCE_FEASIBLE"),
            ("model_performance_inspected", True),
            ("protocol_sha256", "0" * 64),
        ]
        for field, value in alterations:
            with self.subTest(field=field):
                report = copy.deepcopy(self.report)
                report[field] = value
                with self.assertRaises(ValueError):
                    validate_result(report)
        for field, value in [("novel_external_molecules", 30), ("training_rows_scanned", 1),
                             ("excluded_external_molecules_union", 490)]:
            with self.subTest(counter=field):
                report = copy.deepcopy(self.report)
                report["identity_audit"][field] = value
                with self.assertRaises(ValueError):
                    validate_result(report)


if __name__ == "__main__":
    unittest.main()
