import numpy as np
import pytest

from casmi.validate.cohort_bootstrap import joint_cohort_intervals


def test_different_cohort_sizes_preserve_pairing_and_determinism():
    effects = {"internal": {"a": [0, 1, 0, -1]}, "external": {"b": [1, 0, 1]}}
    result = joint_cohort_intervals(effects, n_boot=200, seed=13)
    assert result == joint_cohort_intervals(effects, n_boot=200, seed=13)
    assert result["cohort_counts"] == {"internal": 4, "external": 3}
    assert result["family_size"] == 2
    assert result["outcomes"]["a"]["delta"] == 0
    assert result["outcomes"]["b"]["delta"] == pytest.approx(2 / 3)


def test_constant_effects_have_zero_bootstrap_radius():
    result = joint_cohort_intervals({"source": {"effect": [.5] * 3}}, n_boot=100)
    assert result["joint_radius"] == 0
    assert result["outcomes"]["effect"]["ci_low"] == .5


@pytest.mark.parametrize("effects", [{}, {"a": {}}, {"a": {"x": [0]}},
    {"a": {"x": [0, np.nan]}}, {"a": {"x": [0, 2]}},
    {"a": {"x": [0, 1], "y": [0, 1, 0]}},
    {"a": {"x": [0, 1]}, "b": {"x": [0, 1]}}])
def test_invalid_cohorts_are_rejected(effects):
    with pytest.raises(ValueError):
        joint_cohort_intervals(effects)
