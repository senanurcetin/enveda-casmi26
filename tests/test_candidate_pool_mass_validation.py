"""Candidate-pool physical mass and ppm window invariants. Synthetic data only."""
import numpy as np
import pandas as pd
import pytest

from casmi.pool.pool import CandidatePool


def table(masses=(100.0, 200.0, 300.0)):
    return pd.DataFrame({
        "pool_key": ["AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB", "CCCCCCCCCCCCCC"],
        "smiles": ["C", "CC", "CCC"],
        "mass": masses,
        "sources": [1, 1, 1],
    })


@pytest.mark.parametrize("mass", [0, -50, float("inf"), float("-inf"), "invalid"])
def test_invalid_present_candidate_mass_fails_before_ppm(mass):
    t = table()
    t.loc[1, "mass"] = mass
    with pytest.raises(ValueError, match="masses|non-numeric"):
        CandidatePool(t)


def test_valid_mass_window_excludes_array_and_set_ids():
    pool = CandidatePool(table())
    result = pool.window(200.0, tol_ppm=10)
    np.testing.assert_array_equal(result.idx, np.array([1]))
    np.testing.assert_allclose(result.ppm, [0.0])
    assert pool.window(200.0, 10, exclude=np.array(["BBBBBBBBBBBBBB"])).idx.size == 0
    assert pool.window(200.0, 10, exclude={"BBBBBBBBBBBBBB"}).idx.size == 0
    np.testing.assert_array_equal(pool.window(200.0, 10, exclude=set()).idx, [1])


@pytest.mark.parametrize("value", [
    float("nan"), float("inf"), float("-inf"), 0.0, -1.0,
    [], True, "not-a-mass",
])
def test_invalid_query_mass_is_rejected(value):
    with pytest.raises(ValueError, match="Query mass"):
        CandidatePool(table()).window(value)


@pytest.mark.parametrize("value", [
    float("nan"), float("inf"), float("-inf"), -1.0,
    [], True, "not-ppm",
])
def test_invalid_mass_tolerance_is_rejected(value):
    with pytest.raises(ValueError, match="Mass tolerance ppm"):
        CandidatePool(table()).window(200, tol_ppm=value)


def test_missing_mass_records_still_drop_but_not_impossible_mass():
    t = table()
    t.loc[1, "mass"] = np.nan
    pool = CandidatePool(t)
    assert len(pool) == 2


def test_required_columns_fail_informatively():
    with pytest.raises(ValueError, match="missing columns"):
        CandidatePool(table().drop(columns=["smiles"]))
