"""Synthetic numerical checks for paired simultaneous novelty-mixture inference.

These are numerical fixtures, never reported as Kaggle or scientific results.
"""
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.simultaneous_ci import (  # noqa: E402
    compare_mixtures, paired_novelty_deltas, DEFAULT_GRID,
)


def ranks(n=60, baseline=0, candidate=1):
    a = {cls: [baseline] * n for cls in ("C1", "C2", "C3")}
    b = {cls: [candidate] * n for cls in ("C1", "C2", "C3")}
    return a, b


class SimultaneousInferenceTests(unittest.TestCase):
    def test_positive_and_negative_signal(self):
        a, b = ranks()
        result = compare_mixtures(a, b, n_boot=200, seed=4)
        self.assertEqual(result["verdict"], "ACCEPT")
        self.assertEqual(len(result["scenarios"]), len(DEFAULT_GRID))
        self.assertEqual(result["n_molecules"], 60)
        self.assertEqual(compare_mixtures(b, a, n_boot=200)["verdict"], "REJECT")

    def test_exact_tie_cannot_claim_improvement(self):
        a, b = ranks(candidate=0)
        result = compare_mixtures(a, b, n_boot=200)
        self.assertEqual(result["verdict"], "NOT PROVEN")
        self.assertEqual(result["half_width"], 0.0)

    def test_reproducible_and_simultaneous_wider_than_pointwise(self):
        a, b = ranks(n=75)
        # C1 effects are heterogeneous; mixture draws remain paired.
        b["C1"] = [1 if i % 3 else 0 for i in range(75)]
        b["C2"] = [1 if i % 2 else 0 for i in range(75)]
        b["C3"] = [1 if i % 5 else 0 for i in range(75)]
        r = compare_mixtures(a, b, n_boot=300, seed=1234)
        self.assertEqual(r, compare_mixtures(a, b, n_boot=300, seed=1234))
        effects = paired_novelty_deltas(a, b)
        weights = np.array([[f1, 1 - f1 - f3, f3] for f1, f3 in DEFAULT_GRID])
        means = effects.mean(axis=0) @ weights.T
        idx = np.random.default_rng(1234).integers(0, len(effects), (300, len(effects)))
        sampled = effects[idx].mean(axis=1) @ weights.T
        largest_pointwise = max(
            np.quantile(np.abs(sampled[:, i] - means[i]), 0.95)
            for i in range(len(DEFAULT_GRID))
        )
        self.assertGreaterEqual(r["half_width"] + 1e-12, largest_pointwise)

    def test_mixed_signals_do_not_overclaim(self):
        a, b = ranks(n=100)
        b["C1"] = [1] * 100
        b["C2"] = [0] * 100
        a["C3"] = [1] * 100
        b["C3"] = [0] * 100
        self.assertEqual(
            compare_mixtures(a, b, n_boot=200)["verdict"], "NOT PROVEN"
        )

    def test_invalid_grid_and_parameters_rejected(self):
        a, b = ranks(n=5)
        cases = [
            {"grid": []},
            {"grid": [(-0.1, 0.4)]},
            {"grid": [(0.8, 0.3)]},
            {"grid": [(float("nan"), 0.4)]},
            {"n_boot": 99},
            {"seed": -1},
            {"confidence": 1},
        ]
        for kwargs in cases:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    compare_mixtures(a, b, **kwargs)

    def test_mismatched_arrays_and_rank_type_rejected(self):
        a, b = ranks(n=10)
        b["C2"] = [1] * 9
        with self.assertRaisesRegex(ValueError, "Unpaired"):
            compare_mixtures(a, b, n_boot=200)
        b["C2"] = [1.5] * 10
        with self.assertRaisesRegex(ValueError, "integers"):
            compare_mixtures(a, b, n_boot=200)


if __name__ == "__main__":
    unittest.main()
