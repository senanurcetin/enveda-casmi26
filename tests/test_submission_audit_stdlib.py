"""Fast, dependency-free checks for the CASMI26 pre-submission audit."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts.audit_submission import audit_rows

ROOT = Path(__file__).resolve().parents[1]


class SubmissionAuditTests(unittest.TestCase):
    def test_complete_valid_submission(self):
        rows = [
            {"molecule_id": "m1", "smiles": "CCO;O"},
            {"molecule_id": "m2", "smiles": "c1ccccc1"},
        ]
        report = audit_rows(rows, ["m1", "m2"])
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["errors"], 0)
        self.assertEqual(report["candidate_count"]["max"], 2)
        self.assertEqual(report["coverage"]["coverage_fraction"], 1.0)

    def test_coverage_duplicates_and_invalid_slots(self):
        rows = [
            {"molecule_id": "m1", "smiles": "CCO;;CCO"},
            {"molecule_id": "m1", "smiles": "O"},
            {"molecule_id": "extra", "smiles": "C"},
        ]
        report = audit_rows(rows, ["m1", "m2"])
        codes = {item["code"] for item in report["issues"]}
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue({"blank_candidate", "duplicate_smiles",
                         "duplicate_molecule_id", "missing_molecule",
                         "unexpected_molecule"} <= codes)
        self.assertEqual(report["coverage"]["matched_molecules"], 1)

    def test_enforces_top_25(self):
        report = audit_rows([{"molecule_id": "m1", "smiles": ";".join(["C"] * 26)}])
        self.assertIn("over_25", {i["code"] for i in report["issues"]})
        self.assertEqual(report["candidate_count"]["max"], 26)

    def test_chemistry_check_can_be_injected_without_rdkit(self):
        def key(smiles):
            return {"A": "same", "B": "same", "X": None}.get(smiles)
        report = audit_rows([{"molecule_id": "m1", "smiles": "A;B;X"}],
                            chemistry_key=key)
        self.assertEqual(report["invalid_smiles"], 1)
        self.assertEqual(report["duplicate_candidates"], 1)
        self.assertTrue(report["chemistry_checked"])

    def test_cli_roundtrip_and_json_export(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            submission = d / "submission.csv"
            with submission.open("w", encoding="utf-8", newline="") as file:
                writer = csv.writer(file)
                writer.writerows([["molecule_id", "smiles"], ["m1", "CCO"]])
            (d / "expected.txt").write_text("m1\n", encoding="utf-8")
            out = d / "audit.json"
            proc = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "audit_submission.py"),
                 str(submission), "--expected-ids", str(d / "expected.txt"),
                 "--json-out", str(out)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(json.loads(out.read_text())["status"], "PASS")

    def test_cli_rejects_wrong_header(self):
        with tempfile.TemporaryDirectory() as d:
            submission = Path(d) / "bad.csv"
            submission.write_text("smiles,molecule_id\nCCO,m1\n")
            proc = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "audit_submission.py"),
                 str(submission)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(proc.returncode, 2)
            self.assertIn("expected columns", proc.stderr)


if __name__ == "__main__":
    unittest.main()
