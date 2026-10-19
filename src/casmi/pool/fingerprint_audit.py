"""Read-only quality summary for packed CASMI Morgan fingerprints.

A SHA-256 inventory can consistently sign invalid zero-bit vectors. This
audit adds a numeric sanity check without changing or disclosing source rows.
It does not recompute RDKit chemistry from SMILES and cannot certify model
accuracy or upstream dataset provenance.
"""
from __future__ import annotations

import numpy as np

BITS_BYTES = 256


def fingerprint_quality_summary(pool_keys, fp_keys, fps, *, batch_rows: int = 65536) -> dict:
    """Verify row identity, vector format, and count zero vectors in chunks."""
    if type(batch_rows) is not int or not 1 <= batch_rows <= 1000000:
        raise ValueError("batch_rows must be an integer from 1 to 1000000")
    if (not isinstance(fps, np.ndarray) or fps.dtype != np.uint8
            or fps.ndim != 2 or fps.shape[1] != BITS_BYTES):
        raise ValueError("Fingerprint vectors must be packed uint8 Morgan 2048 bits")
    source = np.asarray(pool_keys)
    recorded = np.asarray(fp_keys)
    if source.ndim != 1 or recorded.ndim != 1:
        raise ValueError("Fingerprint keys must be one-dimensional")
    if len(source) != len(recorded) or len(source) != len(fps):
        raise ValueError("Pool keys and fingerprint vectors have different row counts")
    if not np.array_equal(source, recorded):
        raise ValueError("Fingerprint key order differs from candidate pool")

    zero_rows = 0
    for start in range(0, len(fps), batch_rows):
        fragment = fps[start:start + batch_rows]
        zero_rows += int(np.count_nonzero(~fragment.any(axis=1)))
    return {
        "status": "PASS_NONZERO_FINGERPRINTS" if zero_rows == 0
                  else "INVALID_ZERO_FINGERPRINTS",
        "candidate_count": len(fps),
        "zero_fingerprint_rows": zero_rows,
        "zero_fraction": float(zero_rows / len(fps)) if len(fps) else 0.0,
        "method": "Read-only, row-aligned uint8 bits scanned in bounded batches",
        "warning": (
            "Does not authenticate Kaggle source, recalculate chemistry or "
            "prove candidate SMILES are structurally correct"
        ),
    }
