"""Protect registered live evidence before any model-result artifact exists."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import build_purged_external_notebook as builder
from scripts.check_kaggle_publications import source_hash
from scripts.check_purged_external_publication import validate_result


class PurgedPublicationTests(unittest.TestCase):
    def test_live_bundle_is_verified_without_regeneration(self):
        notebook = {"cells": [{"cell_type": "code", "source": "print('synthetic evidence fixture')"}]}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); folder = root / "bundle"; folder.mkdir()
            (root / "reports").mkdir()
            path = folder / "notebook.ipynb"
            path.write_text(json.dumps(notebook), encoding="utf-8")
            before = path.read_bytes()
            (root / "reports/purged_external_execution_history.json").write_text(
                json.dumps({"attempts": [{"source_sha256": source_hash(notebook)}]}), encoding="utf-8")
            with patch.object(builder, "ROOT", root), patch.object(builder, "OUTPUT", folder):
                builder.build()
                self.assertEqual(path.read_bytes(), before)
                notebook["cells"][0]["source"] = "print('altered fixture')"
                path.write_text(json.dumps(notebook), encoding="utf-8")
                with self.assertRaises(ValueError): builder.build()

    def test_completion_label_without_required_evidence_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_result({"status": "COMPLETE_C3_AND_PURGED_REFERENCE_EXTERNAL_HOLDOUT"})


if __name__ == "__main__": unittest.main()
