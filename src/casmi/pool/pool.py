"""A3 — candidate pool: unique structures (metric key) sorted by neutral monoisotopic mass.

Sources are tracked as a bitmask so the ranker can learn source priors (e.g. that structures known
only from the drug-like enveda-180 library are unlikely answers for a natural-product test set)
and so validation can knock out a structure from the pool (Class-3 simulation).
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

SOURCE_BITS = {"train": 1, "enveda180_only": 2, "coconut": 4, "chebi": 8, "lipidmaps": 16, "pubchem": 32}


@dataclass
class Window:
    idx: np.ndarray          # pool row indices
    ppm: np.ndarray          # signed ppm of each candidate relative to the query mass


class CandidatePool:
    def __init__(self, table: pd.DataFrame):
        """table: pool_key, smiles, mass, sources (bitmask) [+ any extra per-structure columns]."""
        t = table.dropna(subset=["pool_key", "mass"]).sort_values("mass", kind="stable")
        t = t.drop_duplicates("pool_key", keep="first").reset_index(drop=True)
        self.table = t
        self.mass = t["mass"].to_numpy(np.float64)
        self.key = t["pool_key"].to_numpy()
        self.smiles = t["smiles"].to_numpy()
        self.row_of = dict(zip(self.key, range(len(t))))

    def smiles_of(self, key):
        return self.smiles[self.row_of[key]]

    def __len__(self):
        return len(self.mass)

    def window(self, query_mass, tol_ppm=10.0, exclude=None) -> Window:
        """All candidates with |ppm(candidate vs query)| <= tol_ppm, minus excluded keys."""
        d = query_mass * tol_ppm * 1e-6
        lo = np.searchsorted(self.mass, query_mass - d, side="left")
        hi = np.searchsorted(self.mass, query_mass + d, side="right")
        idx = np.arange(lo, hi)
        if exclude:
            idx = idx[~np.isin(self.key[idx], list(exclude))]
        return Window(idx, (query_mass - self.mass[idx]) / self.mass[idx] * 1e6)
