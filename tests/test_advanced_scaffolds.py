import numpy as np
import pytest

from casmi.validate.scaffold_study import family_key, partition, joint_intervals


def test_generic_families_block_substitution_and_heteroatom_changes():
    assert family_key("c1ccccc1O") == family_key("c1ccncc1C")
    assert family_key("CCO") == family_key("CCC")
    assert family_key("CC(C)C") != family_key("CCCC")
    assert family_key("not_smiles") is None


def test_partition_unions_metric_equivalents_before_assigning_splits():
    rows = [{"key": str(i), "family": f"family{i}", "metric": f"metric{i}"} for i in range(10)]
    rows[1]["metric"] = rows[0]["metric"]
    result = partition(rows, {"train": 3, "dev": 2, "test": 2}, "frozen")
    assert result == partition(list(reversed(rows)), {"train": 3, "dev": 2, "test": 2}, "frozen")
    ids = [row["metric"] for values in result.values() for row in values]
    assert len(ids) == len(set(ids)) == 7
    with pytest.raises(ValueError, match="Insufficient"):
        partition(rows, {"train": 8, "dev": 1, "test": 1}, "frozen")


def test_joint_bootstrap_retains_pairing_and_covers_zero_constant_effect():
    out = joint_intervals({"zero": np.zeros(20), "gain": np.full(20, .2)}, n_boot=200)
    assert out["family_size"] == 2
    assert abs(out["outcomes"]["zero"]["delta_mrr"]) < 1e-10
    assert abs(out["outcomes"]["gain"]["delta_mrr"] - .2) < 1e-10
    assert out["joint_radius"] < 1e-10
