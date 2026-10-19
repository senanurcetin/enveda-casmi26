"""Chemically consistent neutral parent identities for multi-fragment inputs.

COCONUT supplies some multi-component SMILES strings. Taking the largest
fragment without recomputing its standard InChIKey and monoisotopic mass
creates an inconsistent candidate identity. Candidate pool keys are standard
(non-tautomer-canonical) InChIKey first blocks, as in the train.parquet pool.
"""
from __future__ import annotations


def neutral_parent_identity(smiles: str) -> tuple[str, str, float, str] | None:
    """Return (parent canonical SMILES, standard InChIKey14, exact mass, formula).

    None signifies malformed SMILES, missing InChI support, or a charged
    parent. No placeholder molecule is ever substituted for invalid input.
    Largest *heavy-atom* fragment is selected, preserving the original
    source's parent-fragment rule. Caller must decide which rows to discard.
    """
    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt
    from rdkit.Chem import rdMolDescriptors
    import math

    if not isinstance(smiles, str) or not smiles.strip():
        return None
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return None
        parts = Chem.GetMolFrags(molecule, asMols=True)
        if not parts:
            return None
        parent = max(parts, key=lambda mol: mol.GetNumHeavyAtoms())
        if parent.GetNumHeavyAtoms() == 0 or Chem.GetFormalCharge(parent) != 0:
            return None
        canonical = Chem.MolToSmiles(parent, canonical=True)
        standard_key = Chem.MolToInchiKey(parent)
        mass = float(ExactMolWt(parent))
        formula = rdMolDescriptors.CalcMolFormula(parent)
    except (ValueError, RuntimeError, TypeError):
        return None
    if (not isinstance(standard_key, str) or len(standard_key) < 14
            or not standard_key[:14].isalpha()
            or not math.isfinite(mass) or mass <= 0):
        return None
    return canonical, standard_key[:14], mass, formula
