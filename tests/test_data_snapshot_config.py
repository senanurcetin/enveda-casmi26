from pathlib import Path

import casmi.config as config


def test_explicit_snapshot_overrides_kaggle_mount_discovery(monkeypatch, tmp_path):
    monkeypatch.setenv("CASMI_DATA", str(tmp_path / "purged-snapshot"))
    monkeypatch.setattr(config, "ON_KAGGLE", True)
    def forbidden(*args, **kwargs):
        raise AssertionError("Mounted original competition data must not override an explicit snapshot")
    monkeypatch.setattr(Path, "rglob", forbidden)
    assert config._find_competition_data() == tmp_path / "purged-snapshot"


def test_implicit_kaggle_discovery_still_works(monkeypatch):
    monkeypatch.delenv("CASMI_DATA", raising=False)
    monkeypatch.setattr(config, "ON_KAGGLE", True)
    monkeypatch.setattr(Path, "rglob", lambda self, pattern: iter([Path("/kaggle/input/declared-competition/test.parquet")]))
    assert config._find_competition_data() == Path("/kaggle/input/declared-competition")
