"""Isotope/tautomer metric equivalents must survive external-audit prefiltering."""
from casmi.chem.canon import canonical_ik14, heavy_atom_key, skeleton_key


def test_composition_prefilter_preserves_isotope_tautomer_equivalence():
    keto, isotope_enol = "CC(=O)C", "[13CH3]C(O)=C"
    assert canonical_ik14(keto) == canonical_ik14(isotope_enol)
    assert skeleton_key(keto) != skeleton_key(isotope_enol)
    assert heavy_atom_key(keto) == heavy_atom_key(isotope_enol)
