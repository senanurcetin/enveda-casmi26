"""Regression tests for cross-source candidate identity reconciliation.

Synthetic molecules and annotations only; no CASMI dataset access.
"""
import numpy as np
import pandas as pd
import pytest

from casmi.pool.source_reconciliation import reconcile_source_rows


def records():
    return pd.DataFrame({
        "pool_key": ["AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB",
                     "AAAAAAAAAAAAAA", "CCCCCCCCCCCCCC"],
        "smiles": ["CC", "CO", "CC", "CCC"],
        "mass": [30.04695, 32.02622, 30.046951, 44.0626],
        "sources": [1, 1, 4, 4],
        "np_likeness": [np.nan, np.nan, 0.5, 0.3],
    })


def test_consistent_source_keys_and_small_roundoff_preserve_source_order():
    out, report = reconcile_source_rows(records())
    assert report["status"] == "SOURCE_ROWS_RECONCILED"
    assert report["input_rows"] == 4
    assert report["distinct_pool_keys"] == 3
    assert report["shared_pool_keys"] == 1
    assert report["mass_conflict_groups"] == 0
    assert len(out) == 4
    assert out.iloc[0].sources == 1  # training precedence remains unchanged


def test_material_same_key_mass_collision_rejected_without_leaking_smiles():
    frame = records()
    frame.loc[2, "mass"] = 44.0626
    with pytest.raises(ValueError, match="Conflicting neutral mass for 1"):
        reconcile_source_rows(frame)


def test_explicit_quarantine_removes_all_sources_of_conflicting_identity():
    frame = records()
    frame.loc[2, "mass"] = 44.0626
    accepted, report = reconcile_source_rows(frame, quarantine_conflicts=True)
    assert accepted.pool_key.tolist() == ["BBBBBBBBBBBBBB", "CCCCCCCCCCCCCC"]
    assert report["quarantined_mass_conflict_groups"] == 1
    assert report["quarantined_source_rows"] == 2
    assert report["complete_rows"] == 2
    assert report["distinct_pool_keys"] == 2
    assert report["mass_conflict_policy"] == "quarantine_entire_identity"


def test_quarantine_cannot_publish_empty_pool():
    frame = records().iloc[[0, 2]].copy()
    frame.loc[2, "mass"] += 10
    with pytest.raises(ValueError, match="No complete candidate rows remain"):
        reconcile_source_rows(frame, quarantine_conflicts=True)


def test_incomplete_row_cannot_be_spliced_into_a_different_source_row():
    frame = records()
    frame.loc[0, "mass"] = np.nan
    accepted, report = reconcile_source_rows(frame)
    assert report["incomplete_rows_excluded"] == 1
    group = accepted.loc[accepted.pool_key == "AAAAAAAAAAAAAA"]
    assert len(group) == 1
    assert group.iloc[0].sources == 4
    assert group.iloc[0].smiles == "CC"
    assert group.iloc[0].mass == pytest.approx(30.046951)


@pytest.mark.parametrize("bad_mass", [0, -10, np.inf, -np.inf, "not-a-mass"])
def test_present_invalid_mass_rejected(bad_mass):
    frame = records()
    frame.loc[0, "mass"] = bad_mass
    with pytest.raises(ValueError, match="mass"):
        reconcile_source_rows(frame)


@pytest.mark.parametrize("bad_key", ["bad", "aaaaaaaaaaaaaa", "A" * 15, "12345678901234"])
def test_malformed_source_key_rejected(bad_key):
    frame = records()
    frame.loc[0, "pool_key"] = bad_key
    with pytest.raises(ValueError, match="InChIKey14"):
        reconcile_source_rows(frame)


def test_invalid_and_blank_smiles_rejected():
    frame = records()
    frame.loc[0, "smiles"] = ""
    with pytest.raises(ValueError, match="SMILES"):
        reconcile_source_rows(frame)
    frame.loc[0, "smiles"] = None
    accepted, report = reconcile_source_rows(frame)
    assert len(accepted) == 3
    assert report["incomplete_rows_excluded"] == 1


def test_all_incomplete_sources_cannot_build_a_pool():
    frame = records()
    frame.loc[:, "mass"] = np.nan
    with pytest.raises(ValueError, match="No complete"):
        reconcile_source_rows(frame)


def test_ppm_tolerance_does_not_hide_large_absolute_mass_conflicts():
    frame = records()
    frame.loc[2, "mass"] += 0.02
    with pytest.raises(ValueError, match="Conflicting neutral mass"):
        reconcile_source_rows(frame)
    reconciled, summary = reconcile_source_rows(frame, max_spread_da=0.03)
    assert len(reconciled) == 4
    assert summary["mass_conflict_groups"] == 0
