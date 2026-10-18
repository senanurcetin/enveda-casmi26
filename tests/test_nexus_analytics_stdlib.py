"""NEXUS analytics asset and SQL renderer integrity smoke tests."""
import json
import unittest
from pathlib import Path

from scripts.render_bigquery_sql import render

ROOT = Path(__file__).resolve().parents[1]


class AnalyticsAssetTests(unittest.TestCase):
    def test_bigquery_query_rendering(self):
        sql = (ROOT / "sql/bigquery/library_quality.sql").read_text()
        rendered = render(sql, "my-casmi-project", "casmi26", "train_spectra")
        self.assertIn("my-casmi-project.casmi26.train_spectra", rendered)
        self.assertIn("COUNT(DISTINCT inchikey14)", rendered)
        self.assertNotIn("{{", rendered)

    def test_reject_identifier_injection(self):
        sql = "SELECT * FROM {{PROJECT_ID}}.{{DATASET_ID}}.{{TABLE_NAME}}"
        for parts in [
            ("bad;DROP", "good", "good"),
            ("my-project", "test;DROP", "good"),
            ("my-project", "good", "name with spaces"),
            ("my-project", "good", "DROP TABLE"),
        ]:
            with self.subTest(parts=parts):
                with self.assertRaises(ValueError):
                    render(sql, *parts)

    def test_profile_train_library_partition_is_consistent(self):
        profile = json.loads((ROOT / "reports/data_profile.json").read_text())
        train = profile["train"]
        total = sum(row["spectra"] for row in train["by_library"].values())
        self.assertEqual(total, train["spectra"])
        self.assertGreater(train["structures"], 250000)
        self.assertEqual(
            profile["test_public_located_in_train"]["matched_spectra"],
            profile["test_public"]["spectra"],
        )

    def test_site_contains_approved_logo_and_source_links(self):
        page = (ROOT / "site/index.html").read_text()
        logo = (ROOT / "site/assets/nexus-official.svg").read_text()
        self.assertIn('src="assets/nexus-official.svg"', page)
        self.assertIn("../reports/data_profile.json", page)
        self.assertIn("NOT LIVE DATA", page)
        self.assertTrue(logo.strip().startswith("<svg"))
        self.assertIn("#00E0FF", logo)


if __name__ == "__main__":
    unittest.main()
