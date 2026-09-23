"""Per-molecule engine: neutral mass -> candidates -> channel features -> ranking.

The baseline ranks by D1 (direct library similarity), breaking ties by mass error. Later stages
replace `rank` with the learned ranker; the candidate/feature interface stays the same, so every
channel is validated through exactly the code path the submission uses.
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from casmi.channels.library import LibraryConfig, LibraryIndex
from casmi.chem.masses import neutral_mass
from casmi.spectra.clean import clean_peaks_fast


@dataclass(frozen=True)
class EngineConfig:
    tol_ppm: float = 10.0
    mass_offset_ppm: float = 0.0      # calibrated instrument offset (measured +1.5 ppm on np-examples)
    top_n: int = 25
    library: LibraryConfig = field(default_factory=LibraryConfig)


def molecule_mass(spectra, offset_ppm=0.0):
    masses = [float(neutral_mass(s["precursor_mz"], s["adduct"])) for s in spectra]
    masses = [m for m in masses if np.isfinite(m)]
    m = float(np.median(masses)) if masses else np.nan
    return m / (1 + offset_ppm * 1e-6)


class BaselineEngine:
    def __init__(self, pool, library: LibraryIndex, cfg: EngineConfig = EngineConfig()):
        self.pool, self.library, self.cfg = pool, library, cfg

    def features(self, spectra, excluded_rows=None, excluded_keys=None) -> pd.DataFrame:
        for s in spectra:   # clean once per query, with the library's cleaning settings
            s["mz_c"], s["it_c"] = clean_peaks_fast(s["mz"], s["it"], s["precursor_mz"], self.library.cfg.clean)
        mass = molecule_mass(spectra, self.cfg.mass_offset_ppm)
        win = self.pool.window(mass, self.cfg.tol_ppm, excluded_keys)
        keys = self.pool.key[win.idx]
        d1 = np.array([self.library.score(spectra, k, excluded_rows) for k in keys]).reshape(-1, 3)
        return pd.DataFrame({"pool_key": keys, "ppm": win.ppm, "d1_max": d1[:, 0], "d1_mean": d1[:, 1],
                             "d1_nref": d1[:, 2].astype(int), "n_cand": len(keys)})

    def rank(self, feats: pd.DataFrame, depth=60):
        """Ranked pool keys (deeper than 25, so output deduplication can refill the list)."""
        if feats.empty:
            return []
        order = np.lexsort((np.abs(feats.ppm.to_numpy()), -feats.d1_mean.to_numpy(), -feats.d1_max.to_numpy()))
        return feats.pool_key.to_numpy()[order][:depth].tolist()
