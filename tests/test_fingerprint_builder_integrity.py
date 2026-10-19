"""Synthetic candidate fingerprint integrity checks, no competition data."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts import build_fingerprints


def test_fp_chunk_retains_row_order_and_positive_morgan_bits():
    packed = build_fingerprints._fp_chunk((12, ["CCO", "O=C=O", "c1ccccc1"]))
    assert packed.shape == (3, 256)
    assert packed.dtype == np.uint8
    assert np.all(packed.any(axis=1))
    # Unique molecules should not all have the same packed bits.
    assert not np.array_equal(packed[0], packed[1])


@pytest.mark.parametrize("bad", ["", "not_a_smiles", "*", None, "C1("])
def test_fp_chunk_rejects_bad_smiles_without_leaking_source(bad):
    with pytest.raises(ValueError, match="pool row 41") as error:
        build_fingerprints._fp_chunk((40, ["C", bad, "CC"]))
    if isinstance(bad, str) and bad:
        assert bad not in str(error.value)


def test_build_checks_candidate_mass_before_generating_any_vectors(tmp_path, monkeypatch):
    # New builder uses the same CandidatePool invariants as inference.
    df = pd.DataFrame({
        "pool_key": ["AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB"],
        "smiles": ["CC", "CCO"],
        "mass": [30.0, 0.0],
    })
    monkeypatch.setattr(build_fingerprints, "ARTIFACT_DIR", tmp_path)
    monkeypatch.setattr(build_fingerprints, "POOL_FILE", "pool.parquet")
    monkeypatch.setattr(build_fingerprints, "POOL_FP_FILE", "pool_fp2048.npy")
    monkeypatch.setattr(build_fingerprints, "POOL_FP_KEYS_FILE", "pool_fp_keys.parquet")
    monkeypatch.setattr(build_fingerprints, "POOL_TAG", "")
    monkeypatch.setattr(build_fingerprints.pd, "read_parquet", lambda *a, **k: df)
    with pytest.raises(ValueError, match="finite and strictly positive"):
        build_fingerprints.main()
    assert not list(tmp_path.iterdir())


def test_failed_chunk_cannot_publish_zero_fingerprint_artifacts(tmp_path, monkeypatch):
    df = pd.DataFrame({
        "pool_key": ["AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB"],
        "smiles": ["CC", "invalid_smiles!"],
        "mass": [30.0, 100.0],
    })
    monkeypatch.setattr(build_fingerprints, "ARTIFACT_DIR", tmp_path)
    monkeypatch.setattr(build_fingerprints, "POOL_FILE", "pool.parquet")
    monkeypatch.setattr(build_fingerprints, "POOL_FP_FILE", "pool_fp2048.npy")
    monkeypatch.setattr(build_fingerprints, "POOL_FP_KEYS_FILE", "pool_fp_keys.parquet")
    monkeypatch.setattr(build_fingerprints, "POOL_TAG", "")
    monkeypatch.setattr(build_fingerprints.pd, "read_parquet", lambda *a, **k: df)

    class InProcessExecutor:
        def __init__(self, *args, **kwargs):
            pass
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def map(self, fn, jobs):
            return map(fn, jobs)

    monkeypatch.setattr(build_fingerprints, "ProcessPoolExecutor", InProcessExecutor)
    with pytest.raises(ValueError, match="pool row 1"):
        build_fingerprints.main()
    assert not (tmp_path / "pool_fp2048.npy").exists()
    assert not (tmp_path / "pool_fp_keys.parquet").exists()
    assert not (tmp_path / "pool_fp_provenance.json").exists()
