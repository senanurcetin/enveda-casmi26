from types import SimpleNamespace

import pandas as pd
import pytest

from scripts import make_submission


def spectrum(intensities, label="CCO"):
    return {"precursor_mz": 100.0, "ms2_mzs": [40.0, 60.0],
            "ms2_normalized_intensities": intensities, "adduct": "[M+H]+",
            "ionization_mode": "positive", "normalized_smiles": label}


def test_full_arrays_reject_coarse_collisions_and_preserve_duplicate_targets(tmp_path, monkeypatch):
    real = spectrum([1.0, .5])
    collision = spectrum([1.0, .7], "CCC")
    pd.DataFrame([real, collision]).to_parquet(tmp_path / "train.parquet")
    monkeypatch.setattr(make_submission, "DATA_DIR", tmp_path)
    targets = pd.DataFrame([{**real, "molecule_id": "a"}, {**real, "molecule_id": "b"}])
    answers = make_submission.sanity_answers(targets)
    assert answers.molecule_id.tolist() == ["a", "b"]
    assert answers.smiles.tolist() == ["CCO", "CCO"]


def test_conflicting_labels_and_unmatched_targets_fail_closed(tmp_path, monkeypatch):
    real = spectrum([1.0, .5])
    pd.DataFrame([real, spectrum([1.0, .5], "CCC")]).to_parquet(tmp_path / "train.parquet")
    monkeypatch.setattr(make_submission, "DATA_DIR", tmp_path)
    targets = pd.DataFrame([{**real, "molecule_id": "a"}])
    with pytest.raises(ValueError, match="conflicting canonical labels"):
        make_submission.sanity_answers(targets)
    pd.DataFrame([spectrum([1.0, .7])]).to_parquet(tmp_path / "train.parquet")
    with pytest.raises(ValueError, match="missing for 1"):
        make_submission.sanity_answers(targets)
