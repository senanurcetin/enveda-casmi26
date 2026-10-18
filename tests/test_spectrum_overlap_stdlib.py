"""Regression tests for full-spectrum leakage checks."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.spectrum_overlap import (  # noqa: E402
    coarse_spectrum_bucket,
    full_spectrum_match,
)


class FullSpectrumOverlapTests(unittest.TestCase):
    def test_true_full_peak_and_intensity_match(self):
        self.assertTrue(full_spectrum_match(
            350.15, [40.2, 100.3], [0.25, 1.0],
            350.15, [40.2, 100.3], [0.25, 1.0],
        ))

    def test_bucket_collision_is_not_exact_overlap(self):
        a_mz, b_mz = [40.2, 100.3], [40.2, 101.3]
        self.assertEqual(
            coarse_spectrum_bucket(350.15, a_mz),
            coarse_spectrum_bucket(350.15, b_mz),
        )
        self.assertFalse(full_spectrum_match(
            350.15, a_mz, [0.25, 1.0],
            350.15, b_mz, [0.25, 1.0],
        ))

    def test_intensity_difference_invalidates_match(self):
        self.assertFalse(full_spectrum_match(
            350.15, [40.2, 100.3], [0.25, 1.0],
            350.15, [40.2, 100.3], [0.5, 1.0],
        ))

    def test_precursor_and_array_length_mismatch(self):
        self.assertFalse(full_spectrum_match(350.15, [40.2], [1.0],
                                             350.151, [40.2], [1.0]))
        self.assertFalse(full_spectrum_match(350.15, [40.2], [1.0],
                                             350.15, [40.2, 100.3], [1.0, 0.2]))
        self.assertFalse(full_spectrum_match(350.15, [40.2], [1.0, 2.0],
                                             350.15, [40.2], [1.0]))

    def test_empty_and_non_finite_not_evidence(self):
        self.assertIsNone(coarse_spectrum_bucket(350.1, []))
        self.assertIsNone(coarse_spectrum_bucket(float("nan"), [40.0]))
        self.assertFalse(full_spectrum_match(350.15, [], [],
                                             350.15, [], []))
        self.assertFalse(full_spectrum_match(350.15, [40.2], [float("nan")],
                                             350.15, [40.2], [float("nan")]))


if __name__ == "__main__":
    unittest.main()
