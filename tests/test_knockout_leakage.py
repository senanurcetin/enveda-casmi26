"""C1/C2 leakage regressions: synthetic stores, no competition data required."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from casmi.validate.knockout import (
    build_queries, include_equivalent_source_rows, knockout_for,
)
from casmi.validate.run_alignment import SIGNATURE_VERSION, query_signature


class SyntheticStore:
    def __init__(self):
        self.meta = pd.DataFrame({
            "ingest_lib": ["querylib"] * 4 + ["otherlib"],
            "test_adduct": [True, True, False, False, True],
            "n_peaks": [7, 7, 7, 7, 7],
            "precursor_mz": [300.1] * 5,
            "adduct": ["[M+H]+"] * 5,
            "ce_ev": [20.0] * 5,
            "ionization_mode": ["positive"] * 5,
        })

    def peaks(self, index):
        return np.array([50.0, 100.0, 125.0]), np.array([0.4, 1.0, 0.2])


KEYS = np.asarray(["A", "A", "A", "A-tautomer", "A"])


def test_c1_knockout_removes_unselected_same_key_source_spectra():
    store = SyntheticStore()
    query, = build_queries(store, KEYS, library="querylib", max_spectra=1, seed=3)
    assert len(query.spectra) == 1
    assert len(query.spectrum_rows) == 1  # original sampled query input
    np.testing.assert_array_equal(query.source_library_rows, [0, 1, 2])
    # Source-library same structure has two additional reference acquisitions,
    # including one with an adduct not used in query construction.
    own = knockout_for(query, "C1", {"A": np.array([0, 1, 2, 4])}, {"A"})
    np.testing.assert_array_equal(own.excluded_rows, [0, 1, 2])
    assert 4 not in own.excluded_rows  # Cross-library D1 evidence remains available.


def test_equivalent_tautomer_variants_are_excluded_in_all_classes():
    store = SyntheticStore()
    query, = build_queries(store, KEYS, library="querylib", max_spectra=1)
    original_digest = query_signature(query)
    source_idx = np.flatnonzero((store.meta.ingest_lib == "querylib").to_numpy())
    include_equivalent_source_rows(query, {"A", "A-tautomer"}, source_idx, KEYS)
    np.testing.assert_array_equal(query.source_library_rows, [0, 1, 2, 3])
    assert query_signature(query) != original_digest
    assert len(query.spectra) == 1  # references removed; query unchanged
    refs = {"A": np.array([0, 1, 2, 4]), "A-tautomer": np.array([3])}
    c1 = knockout_for(query, "C1", refs, {"A", "A-tautomer"})
    c2 = knockout_for(query, "C2", refs, {"A", "A-tautomer"})
    c3 = knockout_for(query, "C3", refs, {"A", "A-tautomer"})
    np.testing.assert_array_equal(c1.excluded_rows, [0, 1, 2, 3])
    np.testing.assert_array_equal(c2.excluded_rows, [0, 1, 2, 3, 4])
    np.testing.assert_array_equal(c3.excluded_rows, [0, 1, 2, 3, 4])
    assert c1.excluded_keys == c2.excluded_keys == set()
    assert c3.excluded_keys == {"A", "A-tautomer"}
    assert SIGNATURE_VERSION == "casmi-query-spectra-v2"


def test_source_scope_is_excluded_from_other_library_references():
    store = SyntheticStore()
    query, = build_queries(store, KEYS, library="querylib", max_spectra=1)
    own = np.array([0, 1, 2, 3])
    include_equivalent_source_rows(query, {"A", "A-tautomer"}, own, KEYS)
    c1 = knockout_for(query, "C1", {"A": np.array([0, 1, 2, 4])}, {"A"})
    assert 4 not in c1.excluded_rows


def test_reference_rows_outside_store_are_rejected():
    store = SyntheticStore()
    query, = build_queries(store, KEYS, library="querylib", max_spectra=1)
    with pytest.raises(ValueError, match="outside store"):
        include_equivalent_source_rows(query, {"A"}, [1000], KEYS)
    with pytest.raises(ValueError, match="Expected at least"):
        include_equivalent_source_rows(query, set(), [0, 1], KEYS)
    with pytest.raises(ValueError, match="Unknown knockout class"):
        knockout_for(query, "X", {}, {"A"})
