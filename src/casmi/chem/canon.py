"""Structure identity exactly as the official metric computes it.

The metric canonicalises every SMILES with RDKit 2026.03.3's TautomerEnumerator and compares the
first InChIKey block. Our pool keys, deduplication and validation labels all go through the same
function, so "same key here" means "same key on the leaderboard".
"""
import warnings
from functools import lru_cache

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

EXPECTED_RDKIT = "2026.03.3"
if rdkit.__version__ != EXPECTED_RDKIT:
    warnings.warn(f"RDKit {rdkit.__version__} != metric's {EXPECTED_RDKIT}; InChIKey14 keys may not match the leaderboard.")

RDLogger.DisableLog("rdApp.*")
_ENUMERATOR = rdMolStandardize.TautomerEnumerator()


def ik14_from_mol(mol):
    try:
        return Chem.MolToInchiKey(_ENUMERATOR.Canonicalize(mol))[:14]
    except (ValueError, RuntimeError, TypeError):
        return None


@lru_cache(maxsize=500_000)
def canonical_ik14(smiles):
    """Tautomer-canonical InChIKey14 (metric-equivalent), or None if RDKit cannot parse it."""
    if not isinstance(smiles, str):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except (ValueError, RuntimeError, TypeError):
        return None
    return None if mol is None else ik14_from_mol(mol)


def heavy_atom_key(smiles):
    """Heavy-atom composition, the metric's cheap pre-filter before tautomer canonicalisation."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    counts = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "H":
            counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
    return tuple(sorted(counts.items()))


@lru_cache(maxsize=1_000_000)
def skeleton_key(smiles):
    """Heavy-atom connectivity only: bonds single, no aromaticity, charges, H or stereo.
    Tautomers (and protomers) share it, so two structures can only have the same metric key if
    they share this key. It costs ~0.2 ms against 10–300 ms for tautomer canonicalisation, so it
    serves as an exact pre-filter before the expensive comparison."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    # The official metric compares only InChIKey's first block, which excludes
    # isotopic layers. Clear isotopes before removing explicit isotope hydrogens.
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    rw = Chem.RWMol(Chem.RemoveHs(mol, sanitize=False))
    for atom in rw.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in rw.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    return Chem.MolToSmiles(rw, canonical=True)
