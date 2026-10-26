import json
import tempfile
import unittest
from pathlib import Path

from scripts.research_snapshot import embedded_sources


class SnapshotSafetyTests(unittest.TestCase):
    def test_never_executes_source_or_accepts_traversal(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "notebook.ipynb"
            for source in ("sources = json.loads(open('secret').read())", "sources = malicious()",
                           "sources = json.loads('{\"../escape.py\":\"pass\"}')"):
                path.write_text(json.dumps({"cells": [{"cell_type": "code", "source": source}]}))
                with self.subTest(source=source), self.assertRaises((ValueError, TypeError)):
                    embedded_sources(path)
