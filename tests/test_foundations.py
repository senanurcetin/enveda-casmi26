"""Foundations every later number depends on: masses, adducts, keys and the metric."""
import pandas as pd
import pytest

from casmi.chem.canon import canonical_ik14
from casmi.chem.masses import formula_mass, ion_mz, neutral_mass, parse_adduct
from casmi.validate.metric import mrr, rank_of_truth, summarize


# ---------------------------------------------------------------- masses (reference values from NIST / PubChem)
@pytest.mark.parametrize("formula, mass", [
    ("C8H10N4O2", 194.080376),   # caffeine
    ("C6H12O6", 180.063388),     # glucose
    ("C15H10O4", 254.057909),    # chrysin
])
def test_formula_mass(formula, mass):
    assert formula_mass(formula) == pytest.approx(mass, abs=2e-6)


@pytest.mark.parametrize("adduct, expected", [
    ("[M+H]+", 195.087652),        # caffeine [M+H]+
    ("[M+Na]+", 217.069596),
    ("[M-H]-", 193.073100),
    ("[2M+H]+", 389.168028),
])
def test_ion_mz_caffeine(adduct, expected):
    assert ion_mz(194.080376, adduct) == pytest.approx(expected, abs=2e-5)


def test_every_test_adduct_parses_and_inverts():
    from casmi.chem.masses import TEST_ADDUCTS
    for a in TEST_ADDUCTS:
        assert parse_adduct(a) is not None, a
        assert neutral_mass(ion_mz(348.1, a), a) == pytest.approx(348.1, abs=1e-9)


# ---------------------------------------------------------------- metric, cases copied from the official doctests
SOL = pd.DataFrame({"molecule_id": ["m_1", "m_2"], "smiles": ["CCO", "c1ccccc1"]})


def test_metric_position_is_rank():
    sub = pd.DataFrame({"molecule_id": ["m_1", "m_2"], "smiles": ["CCO;CC(=O)O", "CC(=O)O;c1ccccc1"]})
    assert mrr(SOL, sub) == 0.75


def test_metric_first_hit_only_and_miss_is_zero():
    sub = pd.DataFrame({"molecule_id": ["m_1", "m_2"], "smiles": ["CC(=O)O", "c1ccccc1;c1ccccc1;c1ccccc1"]})
    assert mrr(SOL, sub) == 0.5


def test_metric_ignores_stereo():
    sol = pd.DataFrame({"molecule_id": ["m_1"], "smiles": ["C[C@H](N)C(=O)O"]})
    sub = pd.DataFrame({"molecule_id": ["m_1"], "smiles": ["CC(N)C(=O)O"]})
    assert mrr(sol, sub) == 1.0


def test_unparseable_guess_still_costs_its_slot():
    sub = pd.DataFrame({"molecule_id": ["m_1", "m_2"], "smiles": ["not_a_smiles;CCO", "c1ccccc1"]})
    assert mrr(SOL, sub) == pytest.approx((0.5 + 1.0) / 2)



def test_blank_submission_slot_does_not_improve_mrr():
    # The truth in the third position must score 1/3, not 1/2.
    # Blank entries are invalid for submission but must not be dropped by
    # the local scorer when auditing a malformed ranking.
    sub = pd.DataFrame({
        "molecule_id": ["m_1", "m_2"],
        "smiles": ["CC(=O)O;;CCO", "c1ccccc1"],
    })
    assert mrr(SOL, sub) == pytest.approx((1.0 / 3 + 1.0) / 2)


def test_empty_guess_slots_count_toward_top25_limit():
    sub = pd.DataFrame({
        "molecule_id": ["m_1", "m_2"],
        "smiles": [";" * 25 + "CCO", "c1ccccc1"],
    })
    with pytest.raises(ValueError, match="More than 25 guesses"):
        mrr(SOL, sub)


def test_tautomers_share_a_key():
    # 2-hydroxypyridine and 2-pyridone are the same InChIKey14 after tautomer canonicalisation
    assert canonical_ik14("Oc1ccccn1") == canonical_ik14("O=c1cccc[nH]1")


def test_key_level_rank_and_summary():
    assert rank_of_truth(["A", "B", "C"], "B") == 2
    assert rank_of_truth(["A"], "Z") == 0
    s = summarize([1, 2, 0, 4], n_boot=200)
    assert s.mrr == pytest.approx((1 + .5 + 0 + .25) / 4)


def test_skeleton_key_is_a_valid_prefilter():
    from casmi.chem.canon import skeleton_key
    # tautomers share the skeleton (necessary condition for sharing a metric key)
    assert skeleton_key("Oc1ccccn1") == skeleton_key("O=c1cccc[nH]1")
    assert skeleton_key("CC(=O)C") == skeleton_key("CC(O)=C")
    # constitutional isomers do not
    assert skeleton_key("CCO") != skeleton_key("COC")
