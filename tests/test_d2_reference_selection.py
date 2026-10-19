"""Scientific regression checks for opt-in D2 source-balanced reference selection.

Synthetic indices only. They do not establish a leaderboard improvement.
"""
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import pytest

from casmi.channels.analog import AnalogConfig, AnalogIndex
from casmi.channels.reference_selection import select_reference_rows


def meta():
    # A: both best-ranked references are from the held-out timsTOF
    # source, so the legacy cap=2 omits the independent other library.
    rows = [
        ("A", "querylib", 10, True),
        ("A", "querylib", 9, True),
        ("A", "other", 8, False),
        ("A", "other2", 7, False),
        ("B", "querylib", 9, True),
        ("B", "other", 7, False),
        ("C", "other", 6, False),
    ]
    t = pd.DataFrame(rows, columns=["train_ik14", "ingest_lib", "n_peaks", "is_timstof"])
    t["adduct"] = "[M+H]+"
    t["test_adduct"] = True
    t["ppm"] = 1.0
    t["precursor_mz"] = 300.0
    t["neutral_mass"] = 299.0
    return t


class DummyStore:
    def __init__(self, frame):
        self.meta = frame

    def peaks(self, index):
        return np.array([50, 80, 120, 150, 190, 220], dtype=float), np.array(
            [0.3, 0.8, 1, 0.2, 0.3, 0.1], dtype=float
        )


def test_legacy_is_same_tims_priority_reference_set():
    t = meta()
    legacy = select_reference_rows(t, t.train_ik14.to_numpy(), per_structure_adduct=2)
    assert list(legacy) == [0, 1, 4, 5, 6]


def test_balanced_retains_independent_library_under_same_cap():
    t = meta()
    keys = t.train_ik14.to_numpy()
    balanced = select_reference_rows(
        t, keys, per_structure_adduct=2, policy="source_balanced"
    )
    assert list(balanced) == [0, 2, 4, 5, 6]
    assert len(balanced) == len(select_reference_rows(t, keys, per_structure_adduct=2))


def test_balanced_fills_unused_slots_by_source_first_then_quality():
    t = meta()
    result = select_reference_rows(
        t, t.train_ik14.to_numpy(), per_structure_adduct=3,
        policy="source_balanced",
    )
    assert list(result) == [0, 2, 3, 4, 5, 6]


def test_balanced_cache_separated_from_historical_index_and_revalidated():
    t = meta()
    store = DummyStore(t)
    keys = t.train_ik14.to_numpy()
    with TemporaryDirectory() as directory:
        cache = Path(directory) / "analog_index.npz"
        old = AnalogIndex(store, keys, AnalogConfig(per_structure_adduct=2), cache=cache)
        assert list(old.rows[:3]) == [0, 1, 4]
        balanced_cfg = AnalogConfig(
            per_structure_adduct=2, reference_policy="source_balanced"
        )
        new = AnalogIndex(store, keys, balanced_cfg, cache=cache)
        balanced_cache = Path(directory) / "analog_index_source_balanced_v1.npz"
        assert cache.exists() and balanced_cache.exists()
        assert list(new.rows[:3]) == [0, 2, 4]
        reloaded = AnalogIndex(store, keys, balanced_cfg, cache=cache)
        np.testing.assert_array_equal(reloaded.rows, new.rows)
        np.testing.assert_allclose(reloaded.matrix.toarray(), new.matrix.toarray())

        # Identical-length snapshot with changed eligibility must invalidate.
        t.loc[0, "ppm"] = 50.0
        altered = AnalogIndex(store, keys, balanced_cfg, cache=cache)
        assert 0 not in altered.rows
        assert 1 in altered.rows
        with np.load(balanced_cache, allow_pickle=False) as z:
            np.testing.assert_array_equal(z["rows"], altered.rows)
        # Legacy cached artifact remains unchanged by the balanced variant.
        with np.load(cache, allow_pickle=False) as z:
            assert 0 in z["rows"]


def test_corrupted_balanced_cache_is_rebuilt_without_touching_legacy():
    t = meta()
    store = DummyStore(t)
    keys = t.train_ik14.to_numpy()
    cfg = AnalogConfig(per_structure_adduct=2, reference_policy="source_balanced")
    with TemporaryDirectory() as directory:
        old_cache = Path(directory) / "analog_index.npz"
        AnalogIndex(store, keys, AnalogConfig(per_structure_adduct=2), cache=old_cache)
        old_bytes = old_cache.read_bytes()
        expected = AnalogIndex(store, keys, cfg, cache=old_cache)
        balanced_cache = Path(directory) / "analog_index_source_balanced_v1.npz"
        balanced_cache.write_bytes(b"corrupt-cache")
        rebuilt = AnalogIndex(store, keys, cfg, cache=old_cache)
        np.testing.assert_array_equal(rebuilt.rows, expected.rows)
        assert balanced_cache.stat().st_size > len(b"corrupt-cache")
        assert old_cache.read_bytes() == old_bytes


def test_balanced_cache_invalidates_changed_config_without_row_changes():
    t = meta()
    keys = t.train_ik14.to_numpy()
    with TemporaryDirectory() as directory:
        cache = Path(directory) / "analog_index.npz"
        first = AnalogIndex(
            DummyStore(t), keys,
            AnalogConfig(reference_policy="source_balanced", per_structure_adduct=2),
            cache=cache,
        )
        changed = AnalogIndex(
            DummyStore(t), keys,
            AnalogConfig(reference_policy="source_balanced", per_structure_adduct=2,
                         tol_da=0.03),
            cache=cache,
        )
        np.testing.assert_array_equal(first.rows, changed.rows)
        path = Path(directory) / "analog_index_source_balanced_v1.npz"
        with np.load(path, allow_pickle=False) as z:
            assert str(z["selection_digest"]) == changed._cache_stamp(
                changed.store, keys, changed.rows
            )


def test_balanced_policy_invalid_inputs():
    t = meta()
    keys = t.train_ik14.to_numpy()
    with pytest.raises(ValueError, match="Unknown selection"):
        select_reference_rows(t, keys, policy="impossible")
    with pytest.raises(ValueError, match="positive"):
        select_reference_rows(t, keys, per_structure_adduct=0)
    with pytest.raises(ValueError, match="lengths"):
        select_reference_rows(t, keys[:-1])


def test_non_tims_source_reservation_beats_second_tims_library():
    t = meta()
    t.loc[1, "ingest_lib"] = "querylib2"  # top tims references from 2 sources
    result = select_reference_rows(
        t, t.train_ik14.to_numpy(), policy="source_balanced",
        per_structure_adduct=2,
    )
    assert list(result[:2]) == [0, 2]


def test_empty_eligible_reference_library_builds_valid_index():
    t = meta()
    t["test_adduct"] = False
    with TemporaryDirectory() as directory:
        index = AnalogIndex(
            DummyStore(t),
            t.train_ik14.to_numpy(),
            AnalogConfig(reference_policy="source_balanced"),
            cache=Path(directory) / "analog_index.npz",
        )
        assert index.rows.size == 0
        assert index.matrix.shape[0] == 0
        assert index.search([], 300.0) == {}




def test_fresh_index_provenance_rebuilt_identically_without_cache():
    store = DummyStore(meta())
    keys = store.meta.train_ik14.to_numpy()
    cfg = AnalogConfig(per_structure_adduct=2, reference_policy="legacy")
    first = AnalogIndex(store, keys, cfg, cache=None)
    second = AnalogIndex(store, keys, cfg, cache=None)
    evidence = first.fresh_index_provenance()
    assert evidence["status"] == "REBUILT_WITHOUT_CACHE"
    assert evidence["schema"] == "casmi-d2-fresh-index-v1"
    assert evidence["cache_used"] is False
    assert evidence["n_references"] == len(first.rows)
    assert evidence["n_nonzero"] == first.matrix.nnz
    assert first.cache_origin == "recomputed_without_cache"
    assert second.fresh_index_provenance()["index_sha256"] == evidence["index_sha256"]
    # Different reference policies can deliberately produce different D2
    # features; they should have different content digests.
    balanced = AnalogIndex(
        store, keys,
        AnalogConfig(per_structure_adduct=2, reference_policy="source_balanced"),
        cache=None,
    )
    assert balanced.fresh_index_provenance()["index_sha256"] != evidence["index_sha256"]


def test_cached_legacy_analog_index_cannot_claim_fresh_provenance():
    store = DummyStore(meta())
    keys = store.meta.train_ik14.to_numpy()
    with TemporaryDirectory() as directory:
        path = Path(directory) / "analog_index.npz"
        written = AnalogIndex(store, keys, AnalogConfig(), cache=path)
        assert written.cache_origin == "built_with_cache"
        with pytest.raises(ValueError, match="not freshly rebuilt"):
            written.fresh_index_provenance()
        loaded = AnalogIndex(store, keys, AnalogConfig(), cache=path)
        assert loaded.cache_origin == "loaded_cache"
        with pytest.raises(ValueError, match="not freshly rebuilt"):
            loaded.fresh_index_provenance()


def test_config_default_is_legacy_for_frozen_baseline():
    assert AnalogConfig().reference_policy == "legacy"
