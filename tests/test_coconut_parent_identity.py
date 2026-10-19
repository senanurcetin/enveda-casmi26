"""COCONUT parent normalization: SMILES/key/mass must describe one molecule.

All fixtures are synthetic and license-free. Run locally with pinned RDKit
to check chemistry; no Kaggle sources or CI workflow required.
"""
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from casmi.pool.structure_identity import neutral_parent_identity
from scripts import build_pool


ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
MIXTURE = ASPIRIN + ".[Na+].[Cl-]"
ANION = "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]"


def key(smiles):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))[:14]


def test_parent_identity_recomputes_all_three_fields_consistently():
    parent = neutral_parent_identity(MIXTURE)
    assert parent is not None
    smiles, parent_key, exact_mass = parent
    assert parent_key == key(ASPIRIN)
    assert parent_key != key(MIXTURE)
    assert exact_mass == pytest.approx(ExactMolWt(Chem.MolFromSmiles(ASPIRIN)))
    assert len(Chem.GetMolFrags(Chem.MolFromSmiles(smiles))) == 1
    assert key(smiles) == parent_key


@pytest.mark.parametrize("invalid", [
    "", "not_a_smiles.[Na+]", "[Na+].[Cl-]", ANION, None,
])
def test_invalid_or_charged_parent_never_becomes_placeholder_molecule(invalid):
    assert neutral_parent_identity(invalid) is None


def test_single_component_parent_roundtrip_is_stable():
    parent = neutral_parent_identity(ASPIRIN)
    assert parent is not None
    assert parent[1] == key(ASPIRIN)
    assert parent[2] == pytest.approx(180.0422587, rel=1e-6)


def test_coconut_loader_merges_salt_and_neutral_parent_by_corrected_key(
    tmp_path, monkeypatch,
):
    rows = [
        # Wrong for the de-salted parent: the source key describes the salt.
        (MIXTURE, key(MIXTURE), "C9H8ClNaO4", 0),
        (ASPIRIN, key(ASPIRIN), "C9H8O4", 0),
        # A charged molecular parent is not a neutral candidate.
        (ANION, key(ANION), "C9H7NaO4", 0),
        ("not_a_smiles.[Na+]", "AAAAAAAAAAAAAA", "C9H8O4", 0),
    ]
    csv = pd.DataFrame({
        "identifier": [f"synthetic-{i}" for i in range(len(rows))],
        "canonical_smiles": [r[0] for r in rows],
        "standard_inchi_key": [r[1] for r in rows],
        "molecular_formula": [r[2] for r in rows],
        "exact_molecular_weight": [180.0] * len(rows),
        "formal_charge": [r[3] for r in rows],
        "np_likeness": [0.5] * len(rows),
        "np_classifier_pathway": ["test"] * len(rows),
        "np_classifier_superclass": ["test"] * len(rows),
        "contains_sugar": [False] * len(rows),
    })
    path = tmp_path / "coconut_synthetic.csv"
    csv.to_csv(path, index=False)
    monkeypatch.setattr(build_pool, "COCONUT_ZIP", path)
    result = build_pool.load_coconut()
    assert len(result) == 1
    assert result.iloc[0].pool_key == key(ASPIRIN)
    assert key(result.iloc[0].smiles) == result.iloc[0].pool_key
    assert result.iloc[0].mass == pytest.approx(
        ExactMolWt(Chem.MolFromSmiles(ASPIRIN)), rel=1e-8
    )
