"""Only standard library: no Parquet, Kaggle data, or CI required."""
import os
from pathlib import Path
import tempfile
import unittest

from casmi.spectra.build_safety import fresh_store_target


class FrozenStoreBuilderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.target = self.root / "store"

    def test_publishes_only_after_successful_verified_build(self):
        with fresh_store_target(self.target) as staging:
            self.assertFalse(self.target.exists())
            (staging / "meta.parquet").write_bytes(b"synthetic-meta")
            (staging / "offsets.i64").write_bytes(b"synthetic-offsets")
            self.assertTrue(staging.exists())
        self.assertEqual(
            (self.target / "meta.parquet").read_bytes(), b"synthetic-meta"
        )
        self.assertEqual(
            (self.target / "offsets.i64").read_bytes(), b"synthetic-offsets"
        )
        self.assertFalse(staging.exists())

    def test_existing_snapshot_cannot_be_overwritten(self):
        self.target.mkdir()
        marker = self.target / "historical.dat"
        marker.write_bytes(b"frozen")
        with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
            with fresh_store_target(self.target):
                self.fail("Should not enter builder")
        self.assertEqual(marker.read_bytes(), b"frozen")

    def test_failed_validation_removes_incomplete_stage(self):
        with self.assertRaisesRegex(ValueError, "invalid artifact"):
            with fresh_store_target(self.target) as stage:
                (stage / "mz.f32").write_bytes(b"partial")
                raise ValueError("invalid artifact")
        self.assertFalse(self.target.exists())
        self.assertFalse(list(self.root.glob(".store-building-*")))

    def test_output_created_during_build_is_not_replaced(self):
        with self.assertRaisesRegex(FileExistsError, "appeared during build"):
            with fresh_store_target(self.target) as stage:
                (stage / "meta.parquet").write_bytes(b"unpublished")
                self.target.mkdir()
                (self.target / "another-process").write_bytes(b"retain")
        self.assertEqual(
            (self.target / "another-process").read_bytes(), b"retain"
        )
        self.assertFalse(list(self.root.glob(".store-building-*")))

    def test_dangling_symlink_cannot_be_used_as_output(self):
        link = self.root / "store"
        try:
            link.symlink_to(self.root / "missing-target")
        except (OSError, NotImplementedError):
            self.skipTest("Symlink privilege unavailable")
        self.assertFalse(link.exists())
        with self.assertRaises(FileExistsError):
            with fresh_store_target(link):
                pass

    def test_nested_parent_created_but_target_remains_fresh(self):
        target = self.root / "new" / "experiment-02" / "store"
        with fresh_store_target(target) as stage:
            (stage / "meta.parquet").write_bytes(b"v2")
        self.assertEqual((target / "meta.parquet").read_bytes(), b"v2")


if __name__ == "__main__":
    unittest.main()
