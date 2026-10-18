"""Offline tests for D2 cross-library index-retention metadata audit."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_analog_retention import audit_retention  # noqa: E402


def metadata():
    rows = [
        ("A", "querylib", "[M+H]+", 10, True),
        ("A", "querylib", "[M+H]+", 9, True),
        ("A", "other", "[M+H]+", 8, False),
        ("B", "querylib", "[M+H]+", 8, True),
        ("B", "other", "[M+H]+", 7, False),
        ("C", "other", "[M+H]+", 9, False),
        ("D", "querylib", "[M+Na]+", 9, True),
        ("D", "other", "[M+H]+", 7, False),
    ]
    return pd.DataFrame({
        "train_ik14": [r[0] for r in rows],
        "ingest_lib": [r[1] for r in rows],
        "adduct": [r[2] for r in rows],
        "n_peaks": [r[3] for r in rows],
        "is_timstof": [r[4] for r in rows],
        "test_adduct": [True] * len(rows),
        "ppm": [0.0] * len(rows),
    })


class IndexRetentionAuditTests(unittest.TestCase):
    def test_detects_cross_library_displacement(self):
        result = audit_retention(metadata(), "querylib", per_structure_adduct=2)
        self.assertEqual(result["status"], "GAP_FOUND")
        self.assertEqual(result["cross_library_candidate_pairs"], 2)
        self.assertEqual(result["cross_library_pairs_retained"], 1)
        self.assertEqual(result["cross_library_pairs_lost_by_index"], 1)
        self.assertEqual(result["structures_with_lost_cross_library_pairs"], 1)
        self.assertEqual(result["cross_library_pair_loss_fraction"], 0.5)

    def test_higher_reference_cap_recovers_cross_library_pair(self):
        result = audit_retention(metadata(), "querylib", per_structure_adduct=3)
        self.assertEqual(result["status"], "NO_GAP_FOUND")
        self.assertEqual(result["cross_library_pairs_lost_by_index"], 0)

    def test_filtering_quality_and_adduct_applied(self):
        meta = metadata()
        meta.loc[4, "ppm"] = 50.0  # other-library B reference fails quality
        result = audit_retention(meta, "querylib", per_structure_adduct=2)
        self.assertEqual(result["cross_library_candidate_pairs"], 1)
        self.assertEqual(result["cross_library_pairs_lost_by_index"], 1)

    def test_no_possible_cross_source_pairs(self):
        meta = metadata()
        meta["ingest_lib"] = "querylib"
        result = audit_retention(meta, "querylib", per_structure_adduct=2)
        self.assertIsNone(result["cross_library_pair_loss_fraction"])
        self.assertEqual(result["cross_library_pairs_lost_by_index"], 0)

    def test_invalid_inputs(self):
        meta = metadata()
        for kwargs in [
            {"query_library": ""},
            {"query_library": "querylib", "per_structure_adduct": 0},
        ]:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    audit_retention(meta, **kwargs)
        with self.assertRaises(ValueError):
            audit_retention(meta.drop(columns=["ppm"]), "querylib")

    def test_nan_and_missing_structure_ignored(self):
        meta = metadata()
        meta.loc[0, "train_ik14"] = None
        meta.loc[1, "ppm"] = np.nan
        result = audit_retention(meta, "querylib", per_structure_adduct=2)
        self.assertEqual(result["status"], "NO_GAP_FOUND")


if __name__ == "__main__":
    unittest.main()
