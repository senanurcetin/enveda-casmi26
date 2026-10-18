"""Independent warehouse-vs-file reconciliation regression coverage."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

from scripts.reconcile_bigquery_profile import parse_aggregate, reconcile

ROOT = Path(__file__).resolve().parents[1]


class ReconcileTests(unittest.TestCase):
    def setUp(self):
        self.profile = {
            "train": {
                "spectra": 6,
                "structures": 4,
                "by_library": {
                    "lib-a": {"spectra": 4, "structures": 3},
                    "lib-b": {"spectra": 2, "structures": 2},
                },
            }
        }
        self.observed = {
            "lib-a": {"spectra": 4, "structures": 3},
            "lib-b": {"spectra": 2, "structures": 2},
        }

    def test_match_and_nonadditive_structure_count(self):
        result = reconcile(self.profile, self.observed)
        self.assertEqual(result["status"], "MATCH")
        self.assertEqual(result["aggregated_spectra"], 6)
        self.assertTrue(result["distinct_structure_counts_are_nonadditive"])

    def test_count_drift_is_explicit(self):
        changed = deepcopy(self.observed)
        changed["lib-b"]["spectra"] = 3
        report = reconcile(self.profile, changed)
        self.assertEqual(report["status"], "DRIFT")
        codes = [i["issue"] for i in report["issues"]]
        self.assertIn("spectra_mismatch", codes)
        self.assertIn("total_spectra_mismatch", codes)

    def test_missing_or_extra_source_is_explicit(self):
        observed = {"lib-a": {"spectra": 4, "structures": 3}}
        report = reconcile(self.profile, observed)
        self.assertIn("missing_source_library", [i["issue"] for i in report["issues"]])
        observed["new-source"] = {"spectra": 2, "structures": 2}
        report = reconcile(self.profile, observed)
        self.assertIn("unseen_source_library", [i["issue"] for i in report["issues"]])

    def test_csv_rejects_duplicate_libraries(self):
        with tempfile.TemporaryDirectory() as path:
            p = Path(path) / "in.csv"
            p.write_text(
                "source_library,spectrum_count,distinct_2d_structures\n"
                "A,2,1\nA,3,2\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "duplicate"):
                parse_aggregate(p)

    def test_csv_rejects_impossible_counts(self):
        with tempfile.TemporaryDirectory() as path:
            p = Path(path) / "in.csv"
            p.write_text(
                "source_library,spectrum_count,distinct_2d_structures\n"
                "A,2,3\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "impossible"):
                parse_aggregate(p)

    def test_cli_historical_snapshot_roundtrip(self):
        profile = json.loads((ROOT / "reports/data_profile.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            csv_path, report_path = root / "aggregate.csv", root / "result.json"
            with csv_path.open("w", newline="", encoding="utf-8") as file:
                writer = csv.writer(file)
                writer.writerow(("source_library", "spectrum_count", "distinct_2d_structures"))
                writer.writerows(
                    (name, row["spectra"], row["structures"])
                    for name, row in profile["train"]["by_library"].items()
                )
            proc = subprocess.run(
                [sys.executable, str(ROOT / "scripts/reconcile_bigquery_profile.py"),
                 "--profile", str(ROOT / "reports/data_profile.json"),
                 "--aggregate", str(csv_path), "--json-out", str(report_path)],
                capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(json.loads(report_path.read_text())["status"], "MATCH")


if __name__ == "__main__":
    unittest.main()
