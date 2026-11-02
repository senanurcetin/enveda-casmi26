import copy
import hashlib
import json
from pathlib import Path
import unittest

from scripts.check_unopened_coverage_publication import validate

ROOT = Path(__file__).resolve().parents[1]


class UnopenedPublicationTests(unittest.TestCase):
    def test_actual_result_and_authenticated_bytes(self):
        raw = (ROOT / 'reports/unopened_coverage_study.json').read_bytes()
        record = json.loads((ROOT / 'reports/unopened_coverage_publication.json').read_text())
        self.assertEqual(hashlib.sha256(raw).hexdigest(), record['aggregate_sha256'])
        validate(json.loads(raw))

    def test_real_result_rejects_evidence_tampering(self):
        result = json.loads((ROOT / 'reports/unopened_coverage_study.json').read_text())
        changes = [
            lambda r: r['source_code_sha256'].clear(),
            lambda r: r['protection'].update(selected_protected_metric_overlap=1),
            lambda r: r['population']['test'].update(n_molecules=299),
            lambda r: r['seal'].update(supervised_models_fitted=1),
            lambda r: r['primary'].update(delta_coverage=.5),
            lambda r: r['secondary']['outcomes'].pop('C1_MRR_optimal_cosine_minus_entropy'),
            lambda r: r['C1']['optimal_cosine'].update(mrr=.99),
            lambda r: r['coverage']['extended'].update(covered=301),
            lambda r: r['raw_inputs']['train.parquet'].update(sha256='0' * 64),
            lambda r: r['limitations'].append({'smiles': 'private fixture'}),
        ]
        for index, change in enumerate(changes):
            with self.subTest(case=index):
                altered = copy.deepcopy(result)
                change(altered)
                with self.assertRaises((ValueError, KeyError)): validate(altered)
