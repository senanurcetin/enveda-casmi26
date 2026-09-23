"""A1 — on-disk spectrum store: cleaned peaks as flat float32 memmaps + a metadata table.

Layout (artifacts/store/):
  mz.f32, it.f32   all peaks back to back
  offsets.i64      spectrum i owns peaks[offsets[i]:offsets[i+1]]
  meta.parquet     one row per spectrum (library, structure, adduct, masses, ppm, CE, ...)
"""
from pathlib import Path

import numpy as np
import pandas as pd


class SpectrumStore:
    def __init__(self, root: Path):
        root = Path(root)
        self.meta = pd.read_parquet(root / "meta.parquet")
        self.offsets = np.fromfile(root / "offsets.i64", dtype=np.int64)
        self.mz = np.memmap(root / "mz.f32", dtype=np.float32, mode="r")
        self.it = np.memmap(root / "it.f32", dtype=np.float32, mode="r")
        assert len(self.offsets) == len(self.meta) + 1, "store is inconsistent"

    def __len__(self):
        return len(self.meta)

    def peaks(self, i):
        a, b = self.offsets[i], self.offsets[i + 1]
        return np.asarray(self.mz[a:b]), np.asarray(self.it[a:b])
