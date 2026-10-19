"""Construct real SpectrumStore instances with synthetic raw arrays.

Mocks only Parquet metadata retrieval; actual byte layout validation and
read-only NumPy memory maps are exercised. No private source data needed.
"""
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from casmi.spectra.store import SpectrumStore


def write_arrays(root, offsets=(0, 2, 2, 4), peaks=(50, 90, 100, 120)):
    np.asarray(offsets, dtype=np.int64).tofile(root / "offsets.i64")
    np.asarray(peaks, dtype=np.float32).tofile(root / "mz.f32")
    np.asarray(peaks, dtype=np.float32).tofile(root / "it.f32")


def test_store_rejects_invalid_metadata_before_returning_peaks():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        write_arrays(root)
        with patch("casmi.spectra.store.pd.read_parquet",
                   return_value=pd.DataFrame({"n_peaks": [2, 2, 0]})):
            with pytest.raises(ValueError, match="n_peaks disagrees"):
                SpectrumStore(root)
        with patch("casmi.spectra.store.pd.read_parquet",
                   return_value=pd.DataFrame({"n_peaks": [2, 0]})):
            with pytest.raises(ValueError, match="row count"):
                SpectrumStore(root)
        with patch("casmi.spectra.store.pd.read_parquet",
                   return_value=pd.DataFrame({"missing": [2, 0, 2]})):
            with pytest.raises(ValueError, match="lacks per-row n_peaks"):
                SpectrumStore(root)


def test_store_peaks_are_correctly_sliced_and_bounds_checked():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        write_arrays(root)
        with patch("casmi.spectra.store.pd.read_parquet",
                   return_value=pd.DataFrame({"n_peaks": [2, 0, 2]})):
            store = SpectrumStore(root)
        with store:
            assert len(store) == 3
            np.testing.assert_array_equal(store.peaks(0)[0], [50, 90])
            np.testing.assert_array_equal(store.peaks(1)[0], [])
            np.testing.assert_array_equal(store.peaks(np.int64(2))[0], [100, 120])
            for row in (-1, 3, 4):
                with pytest.raises(IndexError):
                    store.peaks(row)
            with pytest.raises(TypeError):
                store.peaks(True)
            with pytest.raises(TypeError):
                store.peaks(1.5)
            assert not store.mz.flags.writeable
            assert not store.it.flags.writeable
        with pytest.raises(ValueError, match="closed"):
            store.peaks(0)
        store.close()  # idempotent



def test_store_handles_valid_zero_peak_files_without_memmap_exception():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        write_arrays(root, offsets=(0, 0), peaks=())
        with patch("casmi.spectra.store.pd.read_parquet",
                   return_value=pd.DataFrame({"n_peaks": [0]})):
            store = SpectrumStore(root)
        assert len(store) == 1
        assert store.layout.empty_spectra == 1
        assert store.peaks(0)[0].size == 0
        assert store.peaks(0)[1].size == 0
