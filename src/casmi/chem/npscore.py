"""RDKit natural-product likeness (Ertl et al. 2008, rdkit.Contrib.NP_Score), cached per key.

Computed from the structure for every candidate alike: unlike COCONUT's own np_likeness column it
carries no information about which database a candidate came from.
"""
_MODEL = None
_CACHE = {}


def np_likeness(key, smiles):
    global _MODEL
    v = _CACHE.get(key)
    if v is None:
        from rdkit import Chem
        from rdkit.Contrib.NP_Score import npscorer
        if _MODEL is None:
            _MODEL = npscorer.readNPModel()
        mol = Chem.MolFromSmiles(smiles)
        v = float(npscorer.scoreMol(mol, _MODEL)) if mol is not None else 0.0
        _CACHE[key] = v
    return v
