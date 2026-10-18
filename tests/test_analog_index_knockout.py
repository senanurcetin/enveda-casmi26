"""End-to-end D2 AnalogIndex.search under a simulated knockout.

The reference rankings are synthetic. No Kaggle spectra, labels, or GPU are
needed. The code imports scipy/numba to exercise the real search method.
"""
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.channels.analog import AnalogConfig, AnalogIndex  # noqa: E402


class FakeMatrix:
    def __init__(self, scores):
        self.scores = np.asarray(scores, np.float32)

    def __matmul__(self, other):
        return self.scores.copy()


def index_for(scores):
    instance = AnalogIndex.__new__(AnalogIndex)
    n = len(scores)
    instance.cfg = AnalogConfig(prefilter_top=1, keep_structures=5)
    instance.matrix = FakeMatrix(scores)
    instance.rows = np.arange(20, 20 + n, dtype=np.int64)
    instance.nmass = np.full(n, 300.0, dtype=float)
    instance.adduct = np.array(["[M+H]+"] * n)
    instance.prec = np.full(n, 301.0, dtype=float)
    instance.key = np.array(["self", "real-analog", "fallback"][:n], dtype=object)
    instance._clean = lambda row: (
        np.array([80., 120.], np.float32),
        np.array([0.5, 1.], np.float32),
    )
    return instance


QUERY = [{
    "mz_c": np.array([80., 120.], dtype=np.float32),
    "it_c": np.array([0.5, 1.], dtype=np.float32),
    "precursor_mz": 301.0,
    "adduct": "[M+H]+",
}]


class AnalogIndexKnockoutTests(unittest.TestCase):
    def test_excluded_best_match_does_not_consume_top_slot(self):
        idx = index_for([0.99, 0.8, 0.7])
        with patch("casmi.channels.analog.entropy_similarity", return_value=0.8):
            result = idx.search(QUERY, query_mass=300.0, excluded_rows=np.array([20]))
        self.assertEqual(list(result), ["real-analog"])

    def test_excluded_all_return_empty_without_exception(self):
        idx = index_for([0.99, 0.8, 0.7])
        with patch("casmi.channels.analog.entropy_similarity") as similarity:
            result = idx.search(QUERY, 300., excluded_rows=[20, 21, 22])
            self.assertEqual(result, {})
            similarity.assert_not_called()

    def test_no_reference_rows_return_empty(self):
        idx = index_for([])
        with patch("casmi.channels.analog.entropy_similarity") as similarity:
            self.assertEqual(idx.search(QUERY, 300.0), {})
            similarity.assert_not_called()

    def test_adduct_filter_precedes_topk(self):
        idx = index_for([0.99, 0.8, 0.7])
        idx.adduct[0] = "[M-H]-"
        with patch("casmi.channels.analog.entropy_similarity", return_value=0.75):
            result = idx.search(QUERY, query_mass=300.0)
        self.assertEqual(list(result), ["real-analog"])

    def test_mass_filter_precedes_topk(self):
        idx = index_for([0.99, 0.8, 0.7])
        idx.nmass[0] = 900.
        with patch("casmi.channels.analog.entropy_similarity", return_value=0.75):
            result = idx.search(QUERY, query_mass=300.0)
        self.assertEqual(list(result), ["real-analog"])


if __name__ == "__main__":
    unittest.main()
