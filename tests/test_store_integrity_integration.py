"""Offset and peak-file integrity regressions on synthetic CASMI store files.

NumPy-only, no Parquet reader, labeled spectra, or competition files needed.
"""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.spectra.store_integrity import inspect_store_layout, scan_store_peak_content  # noqa: E402


def make_store(root, *, offsets=(0, 2, 2, 5), peak_count=5):
    np.asarray(offsets, dtype=np.int64).tofile(root / "offsets.i64")
    np.arange(peak_count, dtype=np.float32).tofile(root / "mz.f32")
    np.arange(peak_count, dtype=np.float32).tofile(root / "it.f32")


class StoreLayoutTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        make_store(self.root)

    def test_consistent_layout_and_metadata(self):
        valid = inspect_store_layout(
            self.root, expected_spectra=3, metadata_n_peaks=[2, 0, 3]
        )
        self.assertEqual(valid.spectrum_count, 3)
        self.assertEqual(valid.peak_count, 5)
        self.assertEqual(valid.empty_spectra, 1)
        self.assertEqual(valid.max_peaks_per_spectrum, 3)
        self.assertEqual(valid.summary()["status"], "STRUCTURE_VERIFIED")

    def test_metadata_count_mismatch_does_not_silently_reindex(self):
        with self.assertRaisesRegex(ValueError, "row count"):
            inspect_store_layout(self.root, expected_spectra=4)
        with self.assertRaisesRegex(ValueError, "n_peaks disagrees"):
            inspect_store_layout(self.root, metadata_n_peaks=[2, 1, 2])
        with self.assertRaisesRegex(ValueError, "must align"):
            inspect_store_layout(self.root, metadata_n_peaks=[2, 3])

    def test_invalid_offsets_are_rejected_before_mmap(self):
        for offsets, error in [
            ((1, 2, 2, 5), "First spectrum"),
            ((0, 3, 2, 5), "nondecreasing"),
            ((0, -1, 2, 5), "nonnegative"),
            ((0, 2, 2, 4), "Final spectrum"),
            ((0, 2, 2, 6), "Final spectrum"),
        ]:
            with self.subTest(offsets=offsets):
                make_store(self.root, offsets=offsets)
                with self.assertRaisesRegex(ValueError, error):
                    inspect_store_layout(self.root)
        make_store(self.root)

    def test_raw_file_byte_misalignment(self):
        (self.root / "offsets.i64").write_bytes(b"\x00" * 15)
        with self.assertRaisesRegex(ValueError, "offsets.i64 byte length"):
            inspect_store_layout(self.root)
        make_store(self.root)
        with (self.root / "it.f32").open("ab") as file:
            file.write(b"x")
        with self.assertRaisesRegex(ValueError, "float32-aligned"):
            inspect_store_layout(self.root)
        make_store(self.root)
        np.arange(4, dtype=np.float32).tofile(self.root / "it.f32")
        with self.assertRaisesRegex(ValueError, "different peak counts"):
            inspect_store_layout(self.root)

    def test_invalid_metadata_values(self):
        for bad, error in [
            ([2, None, 3], "missing or invalid"),
            ([2, -1, 4], "missing or invalid"),
            ([2, 0.5, 2.5], "missing or invalid"),
            ([2, 0, "unknown"], "Invalid metadata"),
        ]:
            with self.subTest(counts=bad):
                with self.assertRaisesRegex(ValueError, error):
                    inspect_store_layout(self.root, metadata_n_peaks=bad)

    def test_empty_peak_arrays_are_structurally_valid(self):
        make_store(self.root, offsets=(0, 0, 0), peak_count=0)
        result = inspect_store_layout(
            self.root, expected_spectra=2, metadata_n_peaks=[0, 0]
        )
        self.assertEqual(result.peak_count, 0)
        self.assertEqual(result.empty_spectra, 2)

    def test_optional_scan_flags_numeric_corruption_across_small_chunks(self):
        original = scan_store_peak_content(self.root, chunk_peaks=2)
        self.assertEqual(original["status"], "INVALID_CLEANED_PEAKS")
        # make_store() puts 0 in the first fragment m/z and intensity.
        self.assertEqual(original["nonpositive_mz"], 1)
        self.assertEqual(original["nonpositive_intensity"], 1)
        self.assertEqual(original["invalid_peak_pairs"], 1)

        mz = np.array([10., -5., np.nan, 40., 50.], np.float32)
        it = np.array([1., 2., 3., 0., np.inf], np.float32)
        mz.tofile(self.root / "mz.f32")
        it.tofile(self.root / "it.f32")
        report = scan_store_peak_content(self.root, chunk_peaks=2)
        self.assertEqual(report["status"], "INVALID_CLEANED_PEAKS")
        self.assertEqual(report["total_peaks_scanned"], 5)
        self.assertEqual(report["nonfinite_mz"], 1)
        self.assertEqual(report["nonpositive_mz"], 1)
        self.assertEqual(report["nonfinite_intensity"], 1)
        self.assertEqual(report["nonpositive_intensity"], 1)
        self.assertEqual(report["invalid_peak_pairs"], 4)

    def test_content_scan_clean_empty_and_chunk_validation(self):
        mz = np.array([40., 70., 140., 180., 210.], np.float32)
        it = np.array([0.1, 0.2, 1., 0.3, 0.5], np.float32)
        mz.tofile(self.root / "mz.f32")
        it.tofile(self.root / "it.f32")
        ok = scan_store_peak_content(self.root, chunk_peaks=1)
        self.assertEqual(ok["status"], "VALID_CLEANED_PEAKS")
        self.assertEqual(ok["invalid_peak_pairs"], 0)
        make_store(self.root, offsets=(0, 0, 0), peak_count=0)
        empty = scan_store_peak_content(self.root)
        self.assertEqual(empty["status"], "VALID_CLEANED_PEAKS")
        self.assertEqual(empty["total_peaks_scanned"], 0)
        for invalid in (0, -1, 2.5, True, 1048577):
            with self.subTest(invalid=invalid):
                with self.assertRaisesRegex(ValueError, "chunk_peaks"):
                    scan_store_peak_content(self.root, chunk_peaks=invalid)

    def test_truncated_or_missing_offsets_fail_closed(self):
        (self.root / "offsets.i64").write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "offsets.i64"):
            inspect_store_layout(self.root)
        (self.root / "offsets.i64").unlink()
        with self.assertRaises(FileNotFoundError):
            inspect_store_layout(self.root)


if __name__ == "__main__":
    unittest.main()
