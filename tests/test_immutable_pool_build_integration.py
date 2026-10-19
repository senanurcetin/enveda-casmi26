"""Research snapshots cannot overwrite historical candidate/pool fingerprints.

No Parquet, COCONUT, process pool or Kaggle data are accessed by these tests.
"""
import pytest

from scripts import build_fingerprints, build_pool


def test_pool_rebuild_refuses_existing_candidate_snapshot(tmp_path, monkeypatch):
    existing = tmp_path / "pool.parquet"
    existing.write_bytes(b"frozen-candidate-pool")
    monkeypatch.setattr(build_pool, "ARTIFACT_DIR", tmp_path)
    monkeypatch.setattr(build_pool, "POOL_FILE", "pool.parquet")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        build_pool.main()
    assert existing.read_bytes() == b"frozen-candidate-pool"


@pytest.mark.parametrize("existing_file", [
    "pool_fp2048.npy", "pool_fp_keys.parquet", "pool_fp_provenance.json",
])
def test_fingerprint_rebuild_refuses_any_existing_output(
    tmp_path, monkeypatch, existing_file,
):
    target = tmp_path / existing_file
    target.write_bytes(b"historical-evidence")
    monkeypatch.setattr(build_fingerprints, "ARTIFACT_DIR", tmp_path)
    monkeypatch.setattr(build_fingerprints, "POOL_FP_FILE", "pool_fp2048.npy")
    monkeypatch.setattr(build_fingerprints, "POOL_FP_KEYS_FILE", "pool_fp_keys.parquet")
    monkeypatch.setattr(build_fingerprints, "POOL_TAG", "")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        build_fingerprints.main()
    assert target.read_bytes() == b"historical-evidence"
