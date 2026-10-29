import json
import tempfile
import unittest
from pathlib import Path

from scripts.research_snapshot import embedded_sources


class SnapshotSafetyTests(unittest.TestCase):
    def test_python_and_json_inventory_is_read_as_data(self):
        sources = {'src/module.py': 'pass', 'protocols/design.json': '{"seed": 1}'}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'notebook.ipynb'
            path.write_text(json.dumps({'cells': [{'cell_type': 'code',
                'source': 'sources = json.loads(' + repr(json.dumps(sources)) + ')'}]}))
            self.assertEqual(embedded_sources(path), sources)

    def test_never_executes_source_or_accepts_traversal(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "notebook.ipynb"
            for source in ("sources = json.loads(open('secret').read())", "sources = malicious()",
                           "sources = json.loads('{\"../escape.py\":\"pass\"}')",
                           "sources = json.loads('{\"../escape.json\":\"{}\"}')",
                           "sources = json.loads('{\"script.exe\":\"binary\"}')"):
                path.write_text(json.dumps({"cells": [{"cell_type": "code", "source": source}]}))
                with self.subTest(source=source), self.assertRaises((ValueError, TypeError)):
                    embedded_sources(path)
