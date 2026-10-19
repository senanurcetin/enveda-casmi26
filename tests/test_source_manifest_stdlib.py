"""Offline provenance tests; no spectrum store, Parquet or credentials required."""
import json
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.source_manifest import (  # noqa: E402
    create_manifest, verify_manifest, required_validation_inputs,
    require_validation_coverage,
)


class SourceManifestTests(unittest.TestCase):
    def setUp(self):
        self._temp = tempfile.TemporaryDirectory()
        self.addCleanup(self._temp.cleanup)
        self.root = Path(self._temp.name)
        (self.root / "store").mkdir()
        (self.root / "store" / "meta.parquet").write_bytes(b"fake-parquet-a")
        (self.root / "store" / "mz.f32").write_bytes(b"fake-spectrum-b")
        self.paths = ["store/meta.parquet", "store/mz.f32"]

    def test_inventory_is_content_sensitive_and_path_independent(self):
        manifest = create_manifest(self.root, self.paths[::-1])
        self.assertEqual(
            [f["path"] for f in manifest["files"]],
            sorted(self.paths),
        )
        self.assertEqual(len(manifest["files_digest_sha256"]), 64)
        verified = verify_manifest(self.root, manifest)
        self.assertEqual(verified["status"], "VERIFIED_LOCAL_ARTIFACTS")
        self.assertEqual(verified["files_checked"], 2)
        self.assertNotIn(str(self.root), json.dumps(manifest))
        self.assertNotIn("fake-spectrum", json.dumps(manifest))

    def test_changed_bytes_even_if_same_size_fail_closed(self):
        manifest = create_manifest(self.root, self.paths)
        file = self.root / "store" / "meta.parquet"
        file.write_bytes(b"fake-parquet-b")
        with self.assertRaisesRegex(ValueError, "differ"):
            verify_manifest(self.root, manifest)

    def test_missing_file_and_fake_digest_fail_closed(self):
        manifest = create_manifest(self.root, self.paths)
        modified = deepcopy(manifest)
        modified["files"][0]["bytes"] = 0
        with self.assertRaisesRegex(ValueError, "list/digest"):
            verify_manifest(self.root, modified)
        modified = deepcopy(manifest)
        modified["files"][0]["sha256"] = "z" * 64
        with self.assertRaisesRegex(ValueError, "hash"):
            verify_manifest(self.root, modified)
        (self.root / "store" / "mz.f32").unlink()
        with self.assertRaises(FileNotFoundError):
            verify_manifest(self.root, manifest)

    def test_duplicate_traversal_and_absolute_paths_fail(self):
        bad = [
            ["store/meta.parquet", "store/meta.parquet"],
            ["../outside.bin"],
            ["/tmp/private.bin"],
            ["store/../store/meta.parquet"],
            ["store\\meta.parquet"],
        ]
        for paths in bad:
            with self.subTest(paths=paths):
                with self.assertRaises(ValueError):
                    create_manifest(self.root, paths)

    def test_symlink_outside_root_rejected(self):
        with tempfile.TemporaryDirectory() as outside:
            (Path(outside) / "secret.dat").write_bytes(b"PRIVATE")
            (self.root / "store" / "linked").symlink_to(
                Path(outside) / "secret.dat"
            )
            with self.assertRaisesRegex(ValueError, "outside"):
                create_manifest(self.root, ["store/linked"])

    def test_required_validation_coverage_is_complete_and_tag_aware(self):
        basic = required_validation_inputs()
        self.assertEqual(len(basic), 5)
        self.assertIn("pool.parquet", basic)
        self.assertEqual(len(required_validation_inputs(analog=True)), 7)
        tagged = required_validation_inputs(pool_tag="_pcl", analog=True)
        self.assertIn("pool_pcl.parquet", tagged)
        self.assertIn("pool_fp2048_pcl.npy", tagged)
        with self.assertRaisesRegex(ValueError, "Invalid candidate-pool tag"):
            required_validation_inputs(pool_tag="../escape")
        manifest = {"files": [{"path": name} for name in basic]}
        self.assertEqual(require_validation_coverage(manifest), basic)
        with self.assertRaisesRegex(ValueError, "omits required"):
            require_validation_coverage(manifest, analog=True)
        with self.assertRaisesRegex(ValueError, "omits required"):
            require_validation_coverage(manifest, pool_tag="_pcl")

    def test_empty_manifest_rejected(self):
        with self.assertRaisesRegex(ValueError, "at least one"):
            create_manifest(self.root, [])
        with self.assertRaisesRegex(ValueError, "schema"):
            verify_manifest(self.root, {})
        with self.assertRaisesRegex(ValueError, "artifact entries"):
            verify_manifest(self.root, {"schema": "casmi-local-artifact-manifest-v1", "files": []})


if __name__ == "__main__":
    unittest.main()
