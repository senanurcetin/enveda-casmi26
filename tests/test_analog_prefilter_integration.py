"""Synthetic regressions for analog top-k knockout ordering.

These are algorithm fixtures, not model performance or competition results.
Requires NumPy only.
"""
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.channels.prefilter import eligible_top_k  # noqa: E402


class AnalogPrefilterTests(unittest.TestCase):
    def test_knockout_applied_before_top_k_selection(self):
        # Historically indices 0/1 won the top-2 shortlist and were
        # discarded afterwards, hiding true eligible analogs 2 and 3.
        scores = np.array([0.99, 0.95, 0.85, 0.75])
        allowed = np.array([False, False, True, True])
        np.testing.assert_array_equal(eligible_top_k(scores, allowed, 2), [2, 3])

    def test_same_structure_reference_can_not_displace_unseen_analog(self):
        rows = np.array([10, 11, 12, 13])
        excluded_rows = np.array([10, 11])
        allowed = ~np.isin(rows, excluded_rows)
        top = eligible_top_k(np.array([5, 4, 3, 2.5]), allowed, 1)
        self.assertEqual(rows[top].tolist(), [12])

    def test_tied_scores_break_by_source_index(self):
        scores = np.array([0.7, 0.7, 0.7, 0.9, 0.2])
        mask = np.array([True, True, True, False, True])
        self.assertEqual(eligible_top_k(scores, mask, 2).tolist(), [0, 1])
        self.assertEqual(eligible_top_k(scores, mask, 3).tolist(), [0, 1, 2])

    def test_positive_finite_only(self):
        scores = np.array([np.nan, np.inf, -np.inf, -1, 0, 0.3, 0.1])
        allowed = np.ones(len(scores), bool)
        self.assertEqual(eligible_top_k(scores, allowed, 20).tolist(), [5, 6])

    def test_empty_and_all_excluded_return_empty(self):
        self.assertEqual(eligible_top_k([], np.array([], bool), 10).size, 0)
        self.assertEqual(eligible_top_k(np.array([0.9]), np.array([False]), 1).size, 0)
        self.assertEqual(eligible_top_k(np.array([0.9]), np.array([True]), 0).size, 0)

    def test_input_arrays_remain_intact(self):
        scores = np.array([0.4, 0.9, 0.3], dtype=float)
        allowed = np.array([True, False, True])
        original_scores, original_mask = scores.copy(), allowed.copy()
        self.assertEqual(eligible_top_k(scores, allowed, 2).tolist(), [0, 2])
        np.testing.assert_array_equal(scores, original_scores)
        np.testing.assert_array_equal(allowed, original_mask)

    def test_invalid_parameters_fail_closed(self):
        for scores, allowed, k in [
            (np.array([[1.]]), np.array([True]), 1),
            (np.array([1.]), np.array([True, False]), 1),
            (np.array(["high"]), np.array([True]), 1),
            (np.array([0.3]), np.array([1]), 1),
            (np.array([0.3]), np.array([True]), -1),
            (np.array([0.3]), np.array([True]), 1.5),
        ]:
            with self.subTest(scores=scores, allowed=allowed, k=k):
                with self.assertRaises(ValueError):
                    eligible_top_k(scores, allowed, k)

    def test_large_score_array_matches_reference_full_sort(self):
        rng = np.random.default_rng(23)
        scores = rng.choice(np.array([0., 0.1, 0.3, 0.5, 0.7, 0.9]), size=10000)
        mask = rng.random(len(scores)) > 0.3
        expected = [i for i in sorted(
            np.flatnonzero(mask & (scores > 0)), key=lambda i: (-scores[i], i)
        )[:125]]
        self.assertEqual(eligible_top_k(scores, mask, 125).tolist(), expected)


if __name__ == "__main__":
    unittest.main()
