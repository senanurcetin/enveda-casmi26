"""Molecule-generation overlay integrity tests (synthetic RDKit structures)."""
import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from casmi.pool.pool import CandidatePool


def make_pool():
    return CandidatePool(pd.DataFrame({
        "pool_key": ["AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB"],
        "smiles": ["CC", "CO"],
        "mass": [30.04695, 32.02621],
        "sources": [1, 1],
    }))


def generated_key(smiles):
    return "G" + Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))[:14]


def test_valid_generated_morgan_vectors_and_existing_pool_fps():
    pool = make_pool()
    molecules = ["CCO", "c1ccccc1"]
    keys = [generated_key(smi) for smi in molecules]
    pool.register_generated(keys, molecules)
    assert len(pool.extra_smiles) == 2
    assert all(fp.shape == (256,) and fp.dtype == np.uint8 and fp.any()
               for fp in pool.extra_fp.values())
    original = np.ones((len(pool), 256), dtype=np.uint8)
    result = pool.packed_fps(
        ["AAAAAAAAAAAAAA", keys[0], "BBBBBBBBBBBBBB", keys[1]],
        original,
    )
    assert result.shape == (4, 256)
    np.testing.assert_array_equal(result[0], original[0])
    np.testing.assert_array_equal(result[1], pool.extra_fp[keys[0]])
    np.testing.assert_array_equal(result[2], original[1])


def test_invalid_generated_row_does_not_partially_publish_valid_rows():
    pool = make_pool()
    key = generated_key("CCO")
    with pytest.raises(ValueError, match="row 1") as error:
        pool.register_generated([key, "G" + "C" * 14], ["CCO", "not_a_smiles"])
    assert "not_a_smiles" not in str(error.value)
    assert not pool.extra_smiles and not pool.extra_fp


def test_fake_valid_smiles_key_rejected_and_not_registered():
    pool = make_pool()
    with pytest.raises(ValueError, match="row 0"):
        pool.register_generated(["G" + "A" * 14], ["CCO"])
    assert not pool.extra_smiles


@pytest.mark.parametrize("smi", ["", "[He]", "not a smiles", None])
def test_generated_requires_actual_heavy_atom_structure(smi):
    pool = make_pool()
    with pytest.raises(ValueError, match="row 0"):
        pool.register_generated(["G" + "Z" * 14], [smi])
    assert not pool.extra_fp


def test_unbalanced_generated_identity_lists_rejected_without_writes():
    pool = make_pool()
    with pytest.raises(ValueError, match="counts differ"):
        pool.register_generated([generated_key("CCO")], [])
    assert not pool.extra_fp


def test_repeated_generated_candidate_keeps_first_valid_record():
    pool = make_pool()
    key = generated_key("CCO")
    pool.register_generated([key], ["CCO"])
    original = pool.extra_fp[key].copy()
    pool.register_generated([key], ["CCO"])
    assert len(pool.extra_fp) == 1
    np.testing.assert_array_equal(pool.extra_fp[key], original)


@pytest.mark.parametrize("bad", [
    np.zeros((2, 256), dtype=np.float32),
    np.zeros((1, 256), dtype=np.uint8),
    np.zeros((2, 255), dtype=np.uint8),
])
def test_packed_fingerprint_format_mismatch_never_silently_broadcasts(bad):
    pool = make_pool()
    with pytest.raises(ValueError, match="row-aligned packed"):
        pool.packed_fps(["AAAAAAAAAAAAAA"], bad)


def test_unregistered_key_cannot_use_empty_or_unrelated_fingerprint():
    pool = make_pool()
    with pytest.raises(KeyError, match="absent"):
        pool.packed_fps(["unknown-key"], np.zeros((len(pool), 256), np.uint8))
