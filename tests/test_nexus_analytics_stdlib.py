"""NEXUS analytics asset and SQL renderer integrity smoke tests."""
import json
import unittest
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

from scripts.render_bigquery_sql import render

ROOT = Path(__file__).resolve().parents[1]


class SiteLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if attributes.get("id"):
            self.ids.add(attributes["id"])
        if tag in {"a", "link", "img", "script"}:
            address = attributes.get("href") or attributes.get("src")
            if address:
                self.links.append(address)


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


    def test_static_site_links_resolve_without_network(self):
        site = ROOT / "site/index.html"
        parser = SiteLinks()
        parser.feed(site.read_text(encoding="utf-8"))
        self.assertIn("#overview", parser.links)
        for address in parser.links:
            parts = urlsplit(address)
            if parts.scheme or parts.netloc:
                continue
            with self.subTest(address=address):
                if parts.path:
                    target = (site.parent / unquote(parts.path)).resolve()
                    self.assertTrue(target.is_relative_to(ROOT))
                    self.assertTrue(target.exists())
                elif parts.fragment:
                    self.assertIn(unquote(parts.fragment), parser.ids)


if __name__ == "__main__":
    unittest.main()
