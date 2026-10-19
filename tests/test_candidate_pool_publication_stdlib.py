"""Only standard library: never overwrite published CASMI research candidate pools."""
from pathlib import Path
import tempfile
import unittest

from casmi.pool.build_safety import fresh_pool_file


class CandidatePoolExclusivePublicationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = self.root / "pool_research_v2.parquet"

    def test_new_file_published_after_successful_exit(self):
        with fresh_pool_file(self.output) as staging:
            self.assertFalse(self.output.exists())
            staging.write_bytes(b"fake-parquet-data")
        self.assertEqual(self.output.read_bytes(), b"fake-parquet-data")
        self.assertFalse(staging.exists())

    def test_existing_pool_never_overwritten(self):
        self.output.write_bytes(b"original-immutable")
        with self.assertRaises(FileExistsError):
            with fresh_pool_file(self.output):
                self.fail("Existing files must prevent entering build")
        self.assertEqual(self.output.read_bytes(), b"original-immutable")

    def test_failed_build_cleans_partial_output(self):
        with self.assertRaisesRegex(ValueError, "build-failed"):
            with fresh_pool_file(self.output) as staging:
                staging.write_bytes(b"partial-write")
                raise ValueError("build-failed")
        self.assertFalse(self.output.exists())
        self.assertFalse(list(self.root.glob(".pool_research_v2.parquet-building-*")))

    def test_second_writer_collision_cannot_replace_existing_snapshot(self):
        with self.assertRaises(FileExistsError):
            with fresh_pool_file(self.output) as staging:
                staging.write_bytes(b"candidate-new-data")
                self.output.write_bytes(b"first-writer-wins")
        self.assertEqual(self.output.read_bytes(), b"first-writer-wins")
        self.assertFalse(list(self.root.glob(".pool_research_v2.parquet-building-*")))

    def test_empty_output_is_rejected_and_unpublished(self):
        with self.assertRaisesRegex(ValueError, "empty output"):
            with fresh_pool_file(self.output):
                pass
        self.assertFalse(self.output.exists())

    def test_nested_path_is_allowed_only_when_fresh(self):
        path = self.root / "exp-v3" / "pool.parquet"
        with fresh_pool_file(path) as staging:
            staging.write_bytes(b"fake-parquet")
        self.assertEqual(path.read_bytes(), b"fake-parquet")


if __name__ == "__main__":
    unittest.main()
