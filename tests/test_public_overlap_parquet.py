"""Exercise Arrow filtering and batching against tiny synthetic Parquet inputs."""
import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from scripts.audit_public_overlap import run_audit


def row(**changes):
    result = {"spectrum_id": "synthetic-spectrum", "molecule_id": "synthetic-molecule",
              "precursor_mz": 300.0, "ms2_mzs": [40.0, 70.0],
              "ms2_normalized_intensities": [0.5, 1.0], "adduct": "[M+H]+",
              "ionization_mode": "positive", "ingest_lib": "synthetic-library"}
    result.update(changes)
    return result


def test_exact_precursor_prefilter_is_complete_and_batch_size_independent(tmp_path):
    train, public = tmp_path / "train.parquet", tmp_path / "test.parquet"
    pq.write_table(pa.Table.from_pylist([row(precursor_mz=300.00001), row(),
                                       row(ms2_mzs=[40.0, 71.0]), row(precursor_mz=None)]), train)
    pq.write_table(pa.Table.from_pylist([row()]), public)
    small, large = (run_audit(train, public, batch_size=size) for size in (1, 100))
    assert small["train_rows_scanned"] == 4
    assert small["exact_precursor_candidate_rows"] == 2
    assert small["content_matched_spectra"] == 1
    assert small["coarse_false_match_pairs"] == 1
    assert small["metadata_consistent_matched_spectra"] == 1
    for name in ("content_matched_spectra", "exact_matching_pairs", "inputs", "audit_sources"):
        assert small[name] == large[name]
    serialized = json.dumps(small)
    assert str(tmp_path) not in serialized
    assert "synthetic-spectrum" not in serialized
    assert small["inputs"]["train.parquet"]["bytes"] == train.stat().st_size


def test_same_input_and_missing_schema_are_rejected(tmp_path):
    source = tmp_path / "wrong.parquet"
    pq.write_table(pa.table({"wrong": [1]}), source)
    with pytest.raises(ValueError, match="different"):
        run_audit(source, source)
    other = tmp_path / "other.parquet"
    pq.write_table(pa.table({"wrong": [2]}), other)
    with pytest.raises(ValueError, match="columns"):
        run_audit(source, other)
