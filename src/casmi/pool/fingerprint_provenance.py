"""Content-bound provenance for CASMI candidate-pool Morgan fingerprints.

Fingerprint row order alone does not prove that the bit vectors were computed
from the CURRENT candidate SMILES: a corrected de-salted COCONUT parent may
share a key with an old molecule but have different input structure bytes.
An immutable research snapshot binds pool_key + smiles + mass AND fingerprint
bytes to the same recorded generation, without disclosing source rows.
"""
from __future__ import annotations

import hashlib
import json
import math

import numpy as np

SCHEMA = "casmi-morgan-fingerprint-provenance-v1"
FP_SIZE_BYTES = 256  # radius 2, Morgan 2048 bits


def pool_row_digest(pool) -> str:
    """Hash ordered (key, SMILES, mass) without materializing a JSON dataset."""
    required = ("pool_key", "smiles", "mass")
    if any(col not in pool for col in required):
        raise ValueError("Candidate pool lacks key, SMILES or mass")
    sha = hashlib.sha256()
    sha.update(b"casmi-pool-rows-v1\n")
    for key, smiles, mass in pool.loc[:, list(required)].itertuples(index=False, name=None):
        if not isinstance(key, str) or not key:
            raise ValueError("Invalid candidate pool key")
        if not isinstance(smiles, str) or not smiles:
            raise ValueError("Missing candidate pool SMILES")
        try:
            mass_value = float(mass)
        except (ValueError, TypeError) as exc:
            raise ValueError("Invalid candidate mass") from exc
        if not math.isfinite(mass_value) or mass_value <= 0:
            raise ValueError("Invalid candidate mass")
        line = json.dumps([key, smiles, mass_value.hex()],
                          separators=(",", ":"), ensure_ascii=True)
        sha.update(line.encode("utf-8"))
        sha.update(b"\n")
    return sha.hexdigest()


def _array_digest(fps) -> str:
    if not isinstance(fps, np.ndarray) or fps.dtype != np.uint8:
        raise ValueError("Packed fingerprints must be a NumPy uint8 array")
    if fps.ndim != 2 or fps.shape[1] != FP_SIZE_BYTES:
        raise ValueError("Expected 2048-bit packed fingerprints with 256 columns")
    if not fps.flags.c_contiguous:
        raise ValueError("Packed fingerprint array must be C-contiguous")
    payload = memoryview(fps).cast("B")
    digest = hashlib.sha256()
    for start in range(0, len(payload), 1024 * 1024):
        digest.update(payload[start:start + 1024 * 1024])
    return digest.hexdigest()


def make_fingerprint_provenance(pool, keys, fps) -> dict:
    """Reject mismatched row ordering and produce a content-bound descriptor."""
    names = np.asarray(keys)
    expected = pool["pool_key"].to_numpy()
    if names.ndim != 1 or len(names) != len(expected):
        raise ValueError("Fingerprint-key row count differs from candidate pool")
    if not np.array_equal(names, expected):
        raise ValueError("Fingerprint keys do not match ordered candidate pool keys")
    if len(fps) != len(expected):
        raise ValueError("Fingerprint-vector row count differs from candidate pool")
    return {
        "schema": SCHEMA,
        "n_structures": len(expected),
        "fingerprint": "Morgan radius=2, 2048 bits, packed uint8",
        "pool_rows_sha256": pool_row_digest(pool),
        "fingerprint_bytes_sha256": _array_digest(fps),
        "limitations": (
            "Local byte-alignment attestation only. Does not authenticate "
            "Kaggle provenance, reconstruct chemistry or prove model accuracy."
        ),
    }


def verify_fingerprint_provenance(pool, keys, fps, manifest) -> dict:
    """Fail closed for stale SMILES, reordered keys, mutated bits or old schema."""
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise ValueError("Missing or unsupported fingerprint provenance schema")
    actual = make_fingerprint_provenance(pool, keys, fps)
    for field in ("n_structures", "fingerprint", "pool_rows_sha256",
                  "fingerprint_bytes_sha256"):
        if type(manifest.get(field)) is not type(actual[field]) or manifest[field] != actual[field]:
            raise ValueError(f"Fingerprint provenance mismatch: {field}")
    return {
        "status": "POOL_FINGERPRINTS_VERIFIED",
        "schema": SCHEMA,
        "n_structures": actual["n_structures"],
        "pool_rows_sha256": actual["pool_rows_sha256"],
        "fingerprint_bytes_sha256": actual["fingerprint_bytes_sha256"],
    }
