"""D1 scientific regressions: never compare spectra across ionization modes.

Synthetic metadata and peaks only. No CASMI labels or competition inputs.
"""
import numpy as np
import pandas as pd
import pytest

from casmi.channels.library import LibraryConfig, LibraryIndex


class FakeStore:
    def __init__(self):
        self.meta = pd.DataFrame({
            "ppm": [0.0, 0.0, 0.0],
            "n_peaks": [3, 3, 3],
            "is_timstof": [True, False, True],
            "ionization_mode": ["positive", "negative", "positive"],
            "adduct": ["[M+H]+", "[M-H]-", "[M+H]+"],
            "precursor_mz": [300.0, 302.0, 304.0],
        })
        self.spectra = [
            (np.array([100., 140., 180.]), np.array([1., 0.6, 0.2])),
            (np.array([110., 150., 190.]), np.array([1., 0.5, 0.2])),
            (np.array([125., 165., 205.]), np.array([1., 0.5, 0.2])),
        ]

    def peaks(self, row):
        return self.spectra[row]


def make_index(monkeypatch, cfg=None):
    store = FakeStore()
    idx = LibraryIndex(
        store, np.array(["K", "K", "K"]),
        cfg or LibraryConfig(max_refs_per_candidate=1),
    )
    idx._clean = lambda row: store.peaks(row)
    monkeypatch.setattr(
        "casmi.channels.library.entropy_similarity",
        lambda q_mz, q_it, r_mz, r_it, tol, shift:
            float(np.array_equal(q_mz, r_mz)),
    )
    return idx


def queries():
    return [
        {
            "mode": "positive", "adduct": "[M+H]+",
            "mz_c": np.array([100., 140., 180.]),
            "it_c": np.array([1., 0.6, 0.2]),
        },
        {
            "mode": "negative", "adduct": "[M-H]-",
            "mz_c": np.array([110., 150., 190.]),
            "it_c": np.array([1., 0.5, 0.2]),
        },
    ]


def test_mixed_modes_keep_correct_refs_even_with_one_ref_per_mode(monkeypatch):
    idx = make_index(monkeypatch)
    best, mean, used = idx.score(queries(), "K", excluded=None)
    assert best == 1.0 and mean == 1.0
    assert used == 2

    reversed_result = idx.score(queries()[::-1], "K", excluded=None)
    assert reversed_result == (best, mean, used)


def test_exclusions_apply_separately_to_each_query_mode(monkeypatch):
    idx = make_index(monkeypatch)
    assert idx.score(queries(), "K", np.array([1])) == (1.0, 0.5, 1)
    # The second positive-mode reference remains eligible but has no match;
    # the negative reference matches, so both are counted as used.
    assert idx.score(queries(), "K", np.array([0])) == (1.0, 0.5, 2)
    assert idx.score(queries(), "K", np.array([0, 2])) == (1.0, 0.5, 1)
    assert idx.score(queries(), "K", np.array([0, 1, 2])) == (0.0, 0.0, 0)


def test_unmatched_mode_gets_zero_not_wrong_mode_similarity(monkeypatch):
    idx = make_index(monkeypatch)
    q = queries()
    q[1] = dict(q[1], mode="unsupported-mode", mz_c=q[0]["mz_c"])
    assert idx.score(q, "K", excluded=None) == (1.0, 0.5, 1)


def test_single_mode_and_empty_queries_preserve_prior_behavior(monkeypatch):
    idx = make_index(monkeypatch)
    assert idx.score(queries()[:1], "K", excluded=None) == (1.0, 1.0, 1)
    assert idx.score([], "K", excluded=None) == (0.0, 0.0, 0)
    assert idx.score(queries(), "not-found", excluded=None) == (0.0, 0.0, 0)
