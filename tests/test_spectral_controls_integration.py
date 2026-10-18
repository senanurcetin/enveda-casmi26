"""Deterministic molecule-level MS/MS negative-control regressions.

Synthetic inputs only; diagnostics are not Kaggle validation results.
"""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.spectral_controls import perturb_spectra  # noqa: E402


def input_spectra():
    return [{
        "mz": np.array([70., 100., 130., 170.], np.float32),
        "it": np.array([0.1, 1.0, 0.6, 0.2], np.float32),
        "precursor_mz": 250.0,
        "adduct": "[M+H]+",
        "ce_ev": 25.0,
        "mode": "positive",
    }, {
        "mz": np.array([50., 60., 90.], np.float32),
        "it": np.array([0.8, 0.25, 0.1], np.float32),
        "precursor_mz": 250.0,
        "adduct": "[M+H]+",
        "ce_ev": 40.0,
        "mode": "positive",
    }]


class SpectralControlTests(unittest.TestCase):
    def test_blank_preserves_mass_adduct_and_source(self):
        original = input_spectra()
        before = deepcopy(original)
        controlled = perturb_spectra(original, mode="blank", query_id="lib:m1")
        for a, b in zip(original, controlled):
            self.assertEqual(b["mz"].size, 0)
            self.assertEqual(b["it"].size, 0)
            for key in ("precursor_mz", "adduct", "ce_ev", "mode"):
                self.assertEqual(a[key], b[key])
        for a, b in zip(original, before):
            np.testing.assert_array_equal(a["mz"], b["mz"])
            np.testing.assert_array_equal(a["it"], b["it"])

    def test_shuffle_is_deterministic_per_query_not_order(self):
        original = input_spectra()
        a = perturb_spectra(original, mode="intensity_shuffle",
                            seed=5, query_id="lib:m1")
        b = perturb_spectra(deepcopy(original), mode="intensity_shuffle",
                            seed=5, query_id="lib:m1")
        for x, y, raw in zip(a, b, original):
            np.testing.assert_array_equal(x["it"], y["it"])
            np.testing.assert_array_equal(x["mz"], raw["mz"])
            np.testing.assert_array_equal(np.sort(x["it"]), np.sort(raw["it"]))
            self.assertFalse(np.array_equal(x["it"], raw["it"]))
            self.assertEqual(x["precursor_mz"], raw["precursor_mz"])
        for x, raw in zip(a, original):
            self.assertIsNot(x["it"], raw["it"])

    def test_none_still_deep_copies_peak_arrays(self):
        original = input_spectra()
        clone = perturb_spectra(original, mode="none", query_id="lib:m1")
        for x, y in zip(clone, original):
            np.testing.assert_array_equal(x["it"], y["it"])
            self.assertIsNot(x["mz"], y["mz"])

    def test_modes_and_seed_validated(self):
        raw = input_spectra()
        for kwargs in [
            {"mode": "unknown", "query_id": "lib:m1"},
            {"mode": "blank", "query_id": ""},
            {"mode": "blank", "query_id": "lib:m1", "seed": -1},
            {"mode": "blank", "query_id": "lib:m1", "seed": 1.5},
        ]:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    perturb_spectra(raw, **kwargs)
        malformed = input_spectra()
        malformed[0]["it"] = malformed[0]["it"][:-1]
        with self.assertRaisesRegex(ValueError, "matching 1D"):
            perturb_spectra(malformed, mode="blank", query_id="lib:m1")


if __name__ == "__main__":
    unittest.main()
