"""Protect registered live evidence before any model-result artifact exists."""
import json
import copy
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import build_purged_external_notebook as builder
from scripts.check_kaggle_publications import source_hash
from scripts.check_purged_external_publication import validate_result


class PurgedPublicationTests(unittest.TestCase):
    def test_actual_completed_artifact_and_publication_digest(self):
        root = Path(__file__).resolve().parents[1]
        raw = (root / 'reports/purged_external_study.json').read_bytes()
        publication = json.loads((root / 'reports/purged_external_publication.json').read_text())
        self.assertEqual(hashlib.sha256(raw).hexdigest(), publication['aggregate_artifact_sha256'])
        validate_result(json.loads(raw))

    def test_actual_result_rejects_tampering(self):
        root = Path(__file__).resolve().parents[1]
        original = json.loads((root / 'reports/purged_external_study.json').read_text())
        mutations = [
            ('source', lambda r: r['source_code_sha256'].clear()),
            ('purge', lambda r: r['purged_training'].update(retained_external_family_overlap=1)),
            ('population', lambda r: r['external_population'].update(n_molecules=149)),
            ('seal', lambda r: r['training_seal']['model_files'].pop('external_fp_fold0.pt')),
            ('input', lambda r: r['raw_inputs']['original_train.parquet'].update(sha256='0'*64)),
            ('primary', lambda r: r['primary'].update(delta_coverage=0.5)),
            ('secondary', lambda r: r['secondary_joint_intervals']['outcomes'].pop(next(iter(r['secondary_joint_intervals']['outcomes'])))),
            ('privacy', lambda r: r['limitations'].append({'smiles': 'private fixture'})),
            ('denominator', lambda r: r['models']['external']['expanded']['none']['C2'].update(n=149)),
            ('runtime', lambda r: r['environment'].update(rdkit='2025.3.3')),
        ]
        for name, mutate in mutations:
            with self.subTest(name=name):
                result = copy.deepcopy(original)
                mutate(result)
                with self.assertRaises(ValueError): validate_result(result)

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
