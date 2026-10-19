"""Synthetic, local-only pool-to-Morgan-fingerprint provenance regression tests."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from casmi.pool.fingerprint_provenance import (
    SCHEMA,
    make_fingerprint_provenance,
    pool_row_digest,
    verify_fingerprint_provenance,
)


def fixtures():
    pool = pd.DataFrame({
        "pool_key": ["AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB"],
        "smiles": ["CCO", "O=C=O"],
        "mass": [46.0418648, 43.9898292],
    })
    keys = pool.pool_key.to_numpy()
    fps = np.zeros((2, 256), dtype=np.uint8)
    fps[0, :3] = [1, 5, 8]
    fps[1, :3] = [2, 7, 11]
    return pool, keys, fps


class FingerprintProvenanceTests(unittest.TestCase):
    def test_manifest_verifies_exact_pool_and_fingerprints(self):
        pool, keys, fps = fixtures()
        signed = make_fingerprint_provenance(pool, keys, fps)
        self.assertEqual(signed["schema"], SCHEMA)
        self.assertEqual(signed["n_structures"], 2)
        verified = verify_fingerprint_provenance(pool, keys, fps, signed)
        self.assertEqual(verified["status"], "POOL_FINGERPRINTS_VERIFIED")
        self.assertNotIn("CCO", json.dumps(signed))
        self.assertNotIn("O=C=O", json.dumps(signed))

    def test_same_keys_with_changed_smiles_are_not_compatible(self):
        pool, keys, fps = fixtures()
        signed = make_fingerprint_provenance(pool, keys, fps)
        changed = pool.copy()
        changed.loc[0, "smiles"] = "CO"
        self.assertFalse(pool_row_digest(pool) == pool_row_digest(changed))
        with self.assertRaisesRegex(ValueError, "pool_rows_sha256"):
            verify_fingerprint_provenance(changed, keys, fps, signed)

    def test_mass_or_order_changes_also_invalidate_provenance(self):
        pool, keys, fps = fixtures()
        signed = make_fingerprint_provenance(pool, keys, fps)
        changed_mass = pool.copy()
        changed_mass.loc[0, "mass"] = 50.
        with self.assertRaisesRegex(ValueError, "pool_rows_sha256"):
            verify_fingerprint_provenance(changed_mass, keys, fps, signed)
        with self.assertRaisesRegex(ValueError, "keys do not match"):
            make_fingerprint_provenance(pool, keys[::-1], fps)
        with self.assertRaisesRegex(ValueError, "row count"):
            make_fingerprint_provenance(pool, keys[:1], fps)

    def test_changed_fingerprint_bytes_fail_even_when_keys_and_smiles_match(self):
        pool, keys, fps = fixtures()
        signed = make_fingerprint_provenance(pool, keys, fps)
        corrupted = fps.copy()
        corrupted[0, 0] ^= 1
        with self.assertRaisesRegex(ValueError, "fingerprint_bytes_sha256"):
            verify_fingerprint_provenance(pool, keys, corrupted, signed)

    def test_invalid_schema_and_vector_shapes_rejected(self):
        pool, keys, fps = fixtures()
        signed = make_fingerprint_provenance(pool, keys, fps)
        with self.assertRaisesRegex(ValueError, "schema"):
            verify_fingerprint_provenance(pool, keys, fps, {})
        altered = dict(signed, n_structures=True)
        with self.assertRaisesRegex(ValueError, "n_structures"):
            verify_fingerprint_provenance(pool, keys, fps, altered)
        with self.assertRaisesRegex(ValueError, "256 columns"):
            make_fingerprint_provenance(pool, keys, fps[:, :255])
        with self.assertRaisesRegex(ValueError, "uint8"):
            make_fingerprint_provenance(pool, keys, fps.astype(np.float32))

    def test_mmap_read_only_array_supported(self):
        pool, keys, fps = fixtures()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fps.npy"
            np.save(path, fps)
            mmap_fps = np.load(path, allow_pickle=False, mmap_mode="r")
            try:
                signed = make_fingerprint_provenance(pool, keys, mmap_fps)
                self.assertEqual(
                    verify_fingerprint_provenance(pool, keys, mmap_fps, signed)["n_structures"],
                    2,
                )
            finally:
                # Windows forbids removing a file with an open NumPy mmap.
                mmap_fps._mmap.close()


if __name__ == "__main__":
    unittest.main()
