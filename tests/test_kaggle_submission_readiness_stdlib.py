"""Stdlib-only regression for local competition input readiness."""
import csv
import tempfile
import unittest
from pathlib import Path

from scripts.kaggle_submission_readiness import evaluate


class ReadinessTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.repo = self.root / "repo"
        for name in ("kaggle/submission_driver.py", "kaggle/ours_channel.py",
                     "scripts/audit_submission.py"):
            target = self.repo / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("synthetic", encoding="utf-8")

    def test_missing_competition_files_are_blocked_without_sensitive_paths(self):
        result = evaluate(self.root, repo_root=self.repo)
        self.assertEqual(result["status"], "BLOCKED")
        self.assertEqual(len([x for x in result["blockers"] if "official competition" in x]), 3)
        self.assertFalse(result["submitted"])
        self.assertFalse(result["scored"])
        self.assertNotIn(str(self.root), str(result))

    def test_wrong_submission_header_fails_closed(self):
        (self.root / "sample_submission.csv").write_text(
            "molecule_id,wrong\nm_001,CCO\n", encoding="utf-8"
        )
        result = evaluate(self.root, repo_root=self.repo)
        self.assertIn("Unexpected sample submission columns", result["blockers"])

    def test_missing_notebook_components_are_reported(self):
        (self.repo / "kaggle/ours_channel.py").unlink()
        result = evaluate(self.root, repo_root=self.repo)
        self.assertIn(
            "Missing repository notebook component: kaggle/ours_channel.py",
            result["blockers"],
        )


if __name__ == "__main__":
    unittest.main()
