"""Exact content, denominator and metadata boundaries for public-example overlap."""
import json
import sys
import unittest
from copy import deepcopy
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.overlap_audit import PublicOverlapAudit  # noqa: E402


def spectrum(sid="private-spectrum", molecule="private-molecule", **changes):
    row = {"spectrum_id": sid, "molecule_id": molecule, "precursor_mz": 350.15,
           "ms2_mzs": [40.2, 100.3], "ms2_normalized_intensities": [0.25, 1.0],
           "adduct": "[M+H]+", "ionization_mode": "positive", "ingest_lib": "library-A"}
    row.update(changes)
    return row


class PublicOverlapTests(unittest.TestCase):
    def test_full_arrays_not_coarse_keys_and_complete_molecule_denominator(self):
        first = spectrum()
        second = spectrum("other-private-spectrum", ms2_mzs=[40.2, 200.3])
        audit = PublicOverlapAudit([first, second])
        # Same coarse bucket and precursor, different later peak or intensity.
        collision = spectrum(ms2_normalized_intensities=[0.5, 1.0])
        audit.scan_batch([first, collision], scanned_count=10)
        result = audit.finish(expected_train_rows=10)
        self.assertEqual(result["public_spectra"], 2)
        self.assertEqual(result["content_matched_spectra"], 1)
        self.assertEqual(result["fully_content_matched_molecules"], 0)
        self.assertEqual(result["coarse_false_match_pairs"], 3)
        self.assertEqual(result["train_rows_scanned"], 10)
        text = json.dumps(result)
        self.assertNotIn("private-spectrum", text)
        self.assertNotIn("private-molecule", text)
        self.assertNotIn("350.15", text)

    def test_duplicate_references_do_not_inflate_spectrum_counts(self):
        row = spectrum()
        audit = PublicOverlapAudit([row])
        audit.scan_batch([row, row, dict(row, ingest_lib="library-B")])
        result = audit.finish(expected_train_rows=3)
        self.assertEqual(result["exact_matching_pairs"], 3)
        self.assertEqual(result["content_matched_spectra"], 1)
        self.assertEqual(result["content_matches_by_library"], {"library-A": 1, "library-B": 1})

    def test_numeric_equality_is_distinct_from_acquisition_metadata(self):
        row = spectrum()
        audit = PublicOverlapAudit([row])
        audit.scan_batch([dict(row, ionization_mode="negative")])
        result = audit.finish(expected_train_rows=1)
        self.assertEqual(result["content_matched_spectra"], 1)
        self.assertEqual(result["metadata_consistent_matched_spectra"], 0)
        self.assertEqual(result["metadata_disagreement_pairs"], 1)

    def test_invalid_targets_and_incomplete_scans_fail_closed(self):
        for rows in ([], [spectrum(), spectrum()], [spectrum(sid=None)],
                     [spectrum(ms2_normalized_intensities=[float("nan"), 1])],
                     [spectrum(ms2_normalized_intensities=[1])]):
            with self.subTest(rows=len(rows)), self.assertRaises(ValueError):
                PublicOverlapAudit(rows)
        audit = PublicOverlapAudit([spectrum()])
        with self.assertRaises(ValueError):
            audit.scan_batch([spectrum()], scanned_count=0)
        audit.scan_batch([spectrum()])
        for count in (0, 2, True):
            with self.subTest(expected=count), self.assertRaises(ValueError):
                audit.finish(expected_train_rows=count)

    def test_matching_training_rows_require_library_provenance(self):
        audit = PublicOverlapAudit([spectrum()])
        with self.assertRaises(ValueError):
            audit.scan_batch([spectrum(ingest_lib=None)])


if __name__ == "__main__":
    unittest.main()
