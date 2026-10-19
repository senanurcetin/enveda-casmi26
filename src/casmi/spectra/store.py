"""A1 — on-disk spectrum store: cleaned peaks as flat float32 memmaps + a metadata table.

Layout (artifacts/store/):
  mz.f32, it.f32   all peaks back to back
  offsets.i64      spectrum i owns peaks[offsets[i]:offsets[i+1]]
  meta.parquet     one row per spectrum (library, structure, adduct, masses, ppm, CE, ...)
"""
from pathlib import Path

import operator

import numpy as np
import pandas as pd

from casmi.spectra.store_integrity import inspect_store_layout


class SpectrumStore:
    def __init__(self, root: Path):
        root = Path(root)
        self.meta = pd.read_parquet(root / "meta.parquet")
        if "n_peaks" not in self.meta:
            raise ValueError("Spectrum metadata lacks per-row n_peaks counts")
        layout = inspect_store_layout(
            root, expected_spectra=len(self.meta),
            metadata_n_peaks=self.meta["n_peaks"].to_numpy(),
        )
        self.offsets = layout.offsets
        self.layout = layout
        if layout.peak_count:
            self.mz = np.memmap(root / "mz.f32", dtype=np.float32, mode="r")
            self.it = np.memmap(root / "it.f32", dtype=np.float32, mode="r")
        else:
            # A valid empty peak store cannot be memmapped by NumPy.
            self.mz = np.empty(0, dtype=np.float32)
            self.it = np.empty(0, dtype=np.float32)

    def close(self):
        """Release memory-mapped peak files (important for Windows file locks).

        Any arrays returned by peaks() must not be accessed after close().
        Calling close() repeatedly is safe. The source files are not changed.
        """
        for name in ("mz", "it"):
            array = getattr(self, name, None)
            mapped = getattr(array, "_mmap", None)
            if mapped is not None and not mapped.closed:
                mapped.close()
        self.mz = None
        self.it = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def __len__(self):
        return len(self.meta)

    def peaks(self, i):
        if self.mz is None or self.it is None:
            raise ValueError("SpectrumStore is closed")
        if isinstance(i, (bool, np.bool_)):
            raise TypeError("Spectrum index must be an integer row ID")
        row = operator.index(i)
        if row < 0 or row >= len(self.meta):
            raise IndexError(f"Spectrum row {row} is outside the store")
        a, b = int(self.offsets[row]), int(self.offsets[row + 1])
        return np.asarray(self.mz[a:b]), np.asarray(self.it[a:b])
