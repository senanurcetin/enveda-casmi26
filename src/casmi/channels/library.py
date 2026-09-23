"""D1 — direct library match: compare the query's spectra with each candidate's own reference spectra.

For every candidate structure we look up its library spectra (by metric key), keep the most
comparable ones (same ionisation mode, test adduct, plausible mass), and score the best
query × reference entropy similarity. Knocked-out spectra are never looked at.
"""
from dataclasses import dataclass

import numpy as np

from casmi.spectra.clean import CleanConfig, clean_peaks_fast
from casmi.spectra.similarity import entropy_similarity


@dataclass(frozen=True)
class LibraryConfig:
    tol_da: float = 0.02
    max_ppm: float = 20.0              # reference spectra with worse precursor error are not trusted
    max_refs_per_candidate: int = 40   # cap work per candidate; timsTOF + same adduct first
    clean: CleanConfig = CleanConfig(min_rel_intensity=0.01, top_k_per_window=6, window_da=50, max_peaks=64)


class LibraryIndex:
    """Store rows grouped by metric key, pre-filtered for trustworthiness."""

    def __init__(self, store, key_of_row, cfg: LibraryConfig = LibraryConfig()):
        self.store, self.cfg = store, cfg
        meta = store.meta
        ok = (np.abs(meta.ppm.to_numpy()) <= cfg.max_ppm) & (meta.n_peaks.to_numpy() >= 3)
        rows = np.flatnonzero(ok)
        keys = key_of_row[rows]
        order = np.argsort(keys.astype(str), kind="stable")
        rows, keys = rows[order], keys[order].astype(str)
        starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1]) + 1]
        ends = np.r_[starts[1:], len(keys)]
        self.rows_of = {keys[s]: rows[s:e] for s, e in zip(starts, ends)}
        # priority for picking references: timsTOF first, then same-adduct (applied per query)
        self.is_timstof = meta.is_timstof.to_numpy()
        self.adduct = meta.adduct.astype(str).to_numpy()
        self.mode = meta.ionization_mode.astype(str).to_numpy()
        self._cache = {}

    def _clean(self, row):
        hit = self._cache.get(row)
        if hit is None:
            mz, it = self.store.peaks(row)
            hit = clean_peaks_fast(mz, it, float(self.store.meta.precursor_mz.iat[row]), self.cfg.clean)
            if len(self._cache) < 400_000:
                self._cache[row] = hit
        return hit

    def references(self, key, mode, adducts, excluded):
        rows = self.rows_of.get(key)
        if rows is None:
            return rows
        if excluded is not None and len(excluded):
            rows = rows[~np.isin(rows, excluded)]
        rows = rows[self.mode[rows] == mode]
        if len(rows) > self.cfg.max_refs_per_candidate:
            prio = self.is_timstof[rows].astype(int) * 2 + np.isin(self.adduct[rows], list(adducts)).astype(int)
            rows = rows[np.argsort(-prio, kind="stable")[: self.cfg.max_refs_per_candidate]]
        return rows

    def score(self, query_spectra, key, excluded):
        """Best, mean-of-best-per-spectrum similarity, and number of references used."""
        mode = query_spectra[0]["mode"]
        rows = self.references(key, mode, {s["adduct"] for s in query_spectra}, excluded)
        if rows is None or len(rows) == 0:
            return 0.0, 0.0, 0
        refs = [self._clean(r) for r in rows]
        per_spec = []
        for q in query_spectra:
            per_spec.append(max(entropy_similarity(q["mz_c"], q["it_c"], rm, ri, self.cfg.tol_da, 0.0)
                                for rm, ri in refs))
        return float(max(per_spec)), float(np.mean(per_spec)), len(rows)
