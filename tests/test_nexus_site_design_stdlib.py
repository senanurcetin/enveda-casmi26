"""Offline tests for the static NEXUS CASMI26 research evidence interface.

No browser, network, JS runtime, Kaggle account, private dataset or CI needed.
"""
import json
import unittest
from html.parser import HTMLParser
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


class Outline(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.fragments = []
        self.local_assets = []
    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if "id" in attributes:
            self.ids.add(attributes["id"])
        if tag == "a" and attributes.get("href", "").startswith("#"):
            self.fragments.append(attributes["href"][1:])
        if tag in ("img", "link") and "src" in attributes or tag == "link" and "href" in attributes:
            address = attributes.get("src") or attributes.get("href") or ""
            if address and not address.startswith(("http:", "https:", "//")):
                self.local_assets.append(address)


class ResearchSiteTests(unittest.TestCase):
    def setUp(self):
        self.index = (ROOT / "site" / "index.html").read_text(encoding="utf-8")
        self.style = (ROOT / "site" / "assets" / "research-finish.css").read_text(encoding="utf-8")
        self.evidence = json.loads((ROOT / "site" / "evidence-status.json").read_text(encoding="utf-8"))
        self.parser = Outline()
        self.parser.feed(self.index)

    def test_navigation_and_accessibility(self):
        for expected in ("overview", "models", "data", "quality", "execution", "sources", "main-content"):
            self.assertIn(expected, self.parser.ids)
        for target in self.parser.fragments:
            self.assertIn(target, self.parser.ids, f"Broken navigation anchor: {target}")
        self.assertIn('href="#main-content"', self.index)
        self.assertIn("aria-live=\"polite\"", self.index)
        self.assertIn(":focus-visible", self.style)
        self.assertIn("@media(max-width:650px)", self.style)
        self.assertIn("@media print", self.style)
        self.assertIn("@media(prefers-contrast:more)", self.style)

    def test_logo_is_valid_official_svg_and_asset_references_exist(self):
        root = ET.parse(ROOT / "site" / "assets" / "nexus-official.svg").getroot()
        self.assertTrue(root.tag.endswith("svg"))
        self.assertIn('src="assets/nexus-official.svg"', self.index)
        self.assertIn('href="assets/research-finish.css"', self.index)
        for reference in self.parser.local_assets:
            self.assertTrue((ROOT / "site" / reference).exists(),
                            f"Missing page asset: {reference}")

    def test_research_status_uses_no_undocumented_live_claims(self):
        self.assertEqual(self.evidence["schema"], "nexus-casmi-research-evidence-v1")
        self.assertEqual(self.evidence["reviewed_at"], "2026-10-19")
        self.assertEqual(self.evidence["historical_score"]["value"], 0.422)
        self.assertEqual(self.evidence["historical_score"]["as_of"], "2026-10-08")
        states = {record["id"]: record["state"] for record in self.evidence["gates"]}
        self.assertEqual(states["kaggle"], "blocked")
        self.assertEqual(states["dataset"], "pending")
        self.assertEqual(states["ci"], "pending")
        self.assertEqual(states["leaderboard"], "historical")
        self.assertIn('fetch("evidence-status.json"', self.index)
        self.assertIn("textContent=gate.detail", self.index)
        self.assertNotIn("KGAT_", self.index + self.style + json.dumps(self.evidence))
        self.assertNotIn("live leaderboard", self.index.lower())

    def test_source_documents_resolve_and_static_design_does_not_execute_kaggle(self):
        for path in [
            "docs/KAGGLE_ACCESS_AND_FP_PROVENANCE.md",
            "docs/D2_FRESH_INDEX_EVIDENCE.md",
            "docs/FINGERPRINT_ZERO_VECTOR_AUDIT.md",
            "reports/data_profile.json",
        ]:
            self.assertTrue((ROOT / path).is_file(), path)
        self.assertNotIn("kaggle.api.", self.index)
        self.assertNotIn("fetch(\"https://www.kaggle.com", self.index)


if __name__ == "__main__":
    unittest.main()
