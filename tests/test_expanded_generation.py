import numpy as np
import pandas as pd
import pytest
from types import SimpleNamespace
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from casmi.engine import BaselineEngine, EngineConfig
from casmi.generate.transforms import matching, propose_agg
from casmi.generate.expanded import expanded_proposals
from casmi.validate.budget_features import budget_features


def test_reverse_routes_are_added_without_reordering_historical_routes():
    legacy = matching(30.01056, .001)
    expanded = matching(30.01056, .001, ordered_pairs=True)
    assert expanded[:len(legacy)] == legacy
    paths = {tuple(t.name for t in path) for path in expanded}
    assert ("methyl", "hydroxyl") in paths and ("hydroxyl", "methyl") in paths
    assert len(expanded) > len(legacy)


def test_real_chemistry_expansion_keeps_exact_legacy_features():
    seed = "Oc1ccccc1"
    mass = ExactMolWt(Chem.MolFromSmiles(seed))
    analogs = [(seed, .7, mass)]
    target_mass = mass + 30.01056
    legacy = propose_agg(analogs, target_mass)
    expanded = expanded_proposals(analogs, target_mass)
    assert legacy
    assert set(legacy) <= set(expanded)
    for smiles, record in legacy.items():
        assert {k: v for k, v in expanded[smiles].items() if k != "audit_legacy"} == record
        assert expanded[smiles]["audit_legacy"] == 1


def test_filtered_reference_seeds_do_not_change_the_legacy_budget(monkeypatch):
    captured = {}
    def propose(*args, **kwargs):
        captured.update(kwargs)
        return {}
    monkeypatch.setattr("casmi.generate.expanded.expanded_proposals", propose)
    engine = BaselineEngine.__new__(BaselineEngine)
    engine.cfg = EngineConfig(generate=True, gen_consensus=True, gen_expand=True, gen_max_analogs=30, gen_cap=900)
    engine.pool = SimpleNamespace(row_of={str(i): i for i in range(30) if i not in (0, 2)}, smiles_of=lambda key: "CC" + "C" * int(key))
    analogs = {str(i): (1 - i / 100, 100.0) for i in range(30)}
    assert engine._generated(analogs, 150, set()).empty
    assert len(captured["legacy_analogs"]) == 13
    assert captured["legacy_analogs"][-1][0] == engine.pool.smiles_of("14")


def test_subset_context_is_recomputed_without_origin_shortcuts():
    frame = pd.DataFrame({"pool_key": ["a", "b", "new"], "audit_legacy": [1, 1, 0],
                          "d1_max": [.1, .2, 1], "d3_ll": [1, 2, 10], "d4_int": [.2, .4, 1]})
    legacy = budget_features(frame, "legacy")
    expanded = budget_features(frame, "expanded")
    assert not any(c.startswith("audit_") for c in legacy)
    assert list(legacy.d3_ll_rank) == [2, 1]
    assert list(expanded.d3_ll_rank) == [3, 2, 1]
    np.testing.assert_allclose(legacy.d4_int_z, [-1, 1])
    assert set(legacy.n_cand) == {2}
    assert set(expanded.n_cand) == {3}
    assert list(frame.audit_legacy) == [1, 1, 0]
    with pytest.raises(ValueError):
        budget_features(frame.drop(columns="audit_legacy"), "legacy")


@pytest.mark.parametrize("kwargs", [{"max_analogs": 14}, {"cap": 299}, {"cap": True}])
def test_reduced_expansion_budgets_are_rejected(kwargs):
    with pytest.raises(ValueError):
        expanded_proposals([], 150, **kwargs)


def test_nested_budget_preserves_prior_candidates_under_a_hard_cap():
    from casmi.generate.expanded import preserve_budget
    assert preserve_budget(['a', 'b'], ['b', 'c', 'd'], 3) == ['a', 'b', 'c']
    with pytest.raises(ValueError): preserve_budget(['a', 'b'], ['c'], 1)
