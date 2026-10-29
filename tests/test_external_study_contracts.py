import pytest

from scripts.run_external_c3_study import development_choice, primary_coverage_interval
from scripts.prepare_gnps_external import parse_record


def test_primary_coverage_preserves_paired_gain_loss_accounting():
    result = primary_coverage_interval([0, 1, 0, 0], [1, 1, 0, 0], seed=3, n_boot=200)
    assert result["delta_coverage"] == .25
    assert result["gained"] == 1 and result["lost"] == 0
    assert result == primary_coverage_interval([0, 1, 0, 0], [1, 1, 0, 0], seed=3, n_boot=200)


def test_development_gate_needs_a_real_gain_and_no_loss():
    assert development_choice([0, 1], [1, 1])["selected_policy"] == "expanded"
    assert development_choice([0, 1], [0, 1])["selected_policy"] == "legacy"
    assert development_choice([0, 1], [1, 0])["selected_policy"] == "legacy"
    with pytest.raises(ValueError): development_choice([0], [0, 1])
    with pytest.raises(ValueError): primary_coverage_interval([0, 2], [1, 1], seed=3)


def test_direct_gnps_record_preserves_source_metadata_and_rejects_inference():
    record = {"metadata": {"MSLEVEL": "2", "INSTRUMENT": "ESI QTOF", "IONMODE": "Positive", "PEPMASS": "151",
              "SMILES": "CC", "ADDUCT": "[M+H]+", "PI": "Synthetic fixture"},
              "peaks": [[10 + i, i + 1] for i in range(5)]}
    result = parse_record(record, "SYNTHETIC", "https://example.org/synthetic")
    assert result["authors"] == "Synthetic fixture"
    assert result["ms2_normalized_intensities"][-1] == 100
    assert result["accession"].startswith("GNPS2-SNAPSHOT:")
    record["metadata"]["COMMENT"] = "predicted spectrum"
    with pytest.raises(ValueError, match="unsupported_acquisition"):
        parse_record(record, "SYNTHETIC", "https://example.org/synthetic")
