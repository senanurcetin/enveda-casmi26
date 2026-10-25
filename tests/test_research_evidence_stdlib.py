"""Protect evidence attribution and reject unsupported or unsafe ledger entries."""
import copy
import json
import unittest
from pathlib import Path

from scripts.research_evidence import report, validate

ROOT = Path(__file__).resolve().parents[1]


class ResearchEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.ledger = json.loads((ROOT / "reports/research_evidence.json").read_text())

    def test_export_distinguishes_model_attribution_and_pending_work(self):
        text = report(self.ledger, ROOT)
        self.assertIn("0.309 | own", text)
        self.assertIn("0.422 | fusion", text)
        self.assertIn("revalidation — measured", text)
        self.assertIn("Untouched external holdout — pending", text)
        self.assertIn("not external verification", text)

    def test_unsafe_missing_sources_and_unknown_states_are_rejected(self):
        for change in ({"source": "../outside.md"}, {"source": "missing.md"}, {"state": "verified-live"}):
            ledger = copy.deepcopy(self.ledger)
            ledger["claims"][0].update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate(ledger, ROOT)

    def test_invalid_scores_and_duplicate_ids_are_rejected(self):
        for value in (float("nan"), float("inf"), -0.1, 1.01, True, "0.422"):
            ledger = copy.deepcopy(self.ledger)
            ledger["experiments"][0]["score"] = value
            with self.subTest(score=value), self.assertRaises(ValueError):
                validate(ledger, ROOT)
        self.ledger["claims"][0]["id"] = self.ledger["experiments"][0]["id"]
        with self.assertRaises(ValueError):
            validate(self.ledger, ROOT)


if __name__ == "__main__":
    unittest.main()
