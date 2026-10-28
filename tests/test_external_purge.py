import pandas as pd
import pytest

from casmi.chem.canon import canonical_ik14
from casmi.validate.external_purge import select_external_families, blocked_training_identities
from casmi.validate.scaffold_study import family_key


def test_external_family_selection_is_order_invariant_and_keeps_metric_groups():
    source = pd.DataFrame([
        {"metric_key": metric, "family": family, "accession": accession, "standard_ik14": metric}
        for metric, family, accession in [("a", "f1", "a1"), ("a", "f1", "a2"), ("b", "f1", "b1"),
                                          ("c", "f2", "c1"), ("d", "f3", "d1")]])
    chosen, families = select_external_families(source, "fixed-synthetic-salt", maximum=2, minimum=2)
    shuffled, shuffled_families = select_external_families(source.sample(frac=1, random_state=9), "fixed-synthetic-salt", maximum=2, minimum=2)
    pd.testing.assert_frame_equal(chosen, shuffled)
    assert families == shuffled_families and chosen.metric_key.nunique() == 2
    assert chosen.groupby("family").metric_key.nunique().max() == 1
    with pytest.raises(ValueError, match="Insufficient"):
        select_external_families(source, "fixed-synthetic-salt")
    with pytest.raises(ValueError, match="identity contract"):
        select_external_families(pd.concat([source, source.iloc[:1]]), "fixed-synthetic-salt", maximum=2, minimum=2)


def test_isotopic_tautomer_match_blocks_entire_declared_training_identity():
    source = pd.DataFrame([{"metric_key": canonical_ik14("CC(=O)C"), "standard_ik14": "fixture-source-key", "smiles": "CC(=O)C"}])
    rows = [("ambiguous-key", "[13CH3]C(O)=C"), ("ambiguous-key", "CCCCCC"), ("clean-key", "CCN")]
    blocked = blocked_training_identities(rows, source, [])
    assert "ambiguous-key" in blocked and "clean-key" not in blocked
    assert [row for row in rows if row[0] not in blocked] == [("clean-key", "CCN")]


def test_entire_generic_external_family_is_removed_without_metric_equality():
    source = pd.DataFrame([{"metric_key": canonical_ik14("CCO"), "standard_ik14": "fixture-source-key", "smiles": "CCO"}])
    assert canonical_ik14("CCN") != canonical_ik14("CCO")
    blocked = blocked_training_identities([("same-family", "CCN"), ("other-family", "c1ccccc1")], source, [family_key("CCO")])
    assert "same-family" in blocked and "other-family" not in blocked


def test_undefined_family_quarantines_all_rows_of_its_declared_identity():
    source = pd.DataFrame([{"metric_key": canonical_ik14("CCO"), "standard_ik14": "fixture-source-key", "smiles": "CCO"}])
    rows = [("unassigned-key", "[H+]"), ("unassigned-key", "c1ccccc1"), ("clean-key", "CCN")]
    assert family_key("[H+]") is None
    blocked, audit = blocked_training_identities(rows, source, [], return_audit=True)
    assert "unassigned-key" in blocked and "clean-key" not in blocked
    assert audit == {"unassignable_family_labels_removed": 1, "unassignable_family_training_identities_removed": 1}
