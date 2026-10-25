import pandas as pd
import pytest

from casmi.pool.structure_quality import quarantine_invalid_candidates


def test_quarantine_never_publishes_invalid_or_empty_fingerprints():
    table = pd.DataFrame({"pool_key": ["valid", "invalid", "empty", "hydrogen"],
                          "smiles": ["CCO", "invalid", "", "[H][H]"]})
    retained, report = quarantine_invalid_candidates(table)
    assert retained.pool_key.tolist() == ["valid"]
    assert report["quarantined_candidates"] == 3
    assert sum(report["reasons"].values()) == 3
    assert len(table) == 4
    assert "pool_key" not in report and "smiles" not in report


def test_quarantine_cannot_publish_empty_or_duplicate_identity_pool():
    with pytest.raises(ValueError, match="No candidates remain"):
        quarantine_invalid_candidates(pd.DataFrame({"pool_key": ["bad"], "smiles": ["bad"]}))
    with pytest.raises(ValueError, match="one candidate per identity"):
        quarantine_invalid_candidates(pd.DataFrame({"pool_key": ["same", "same"], "smiles": ["CC", "CO"]}))
