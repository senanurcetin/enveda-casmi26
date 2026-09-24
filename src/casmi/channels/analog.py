"""D2 — mass-shifted analog propagation.

A Class-2 molecule has no reference spectrum, but its structural relatives usually do. Relatives
that differ by one substituent fragment into largely the same ions, either at the same m/z (the
substituent was lost) or shifted by the precursor difference (it was retained). We find such
analog spectra in two stages and propagate their structures onto the candidates:

  1. prefilter: every library spectrum is a sparse vector over fragment-m/z bins and
     neutral-loss bins; one sparse mat-vec scores a query spectrum against the whole library;
  2. refine: the top hits are rescored with mass-shifted entropy similarity;
  3. propagate: candidate score = max over analogs a of  sim(a)^p · Tanimoto(fp_c, fp_a).
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from casmi.spectra.clean import CleanConfig, clean_peaks_fast
from casmi.spectra.similarity import entropy_similarity

BIN_DA = 0.05
MAX_MZ = 1300.0
N_BINS = int(MAX_MZ / BIN_DA) + 1


@dataclass(frozen=True)
class AnalogConfig:
    max_shift_da: float = 200.0
    prefilter_top: int = 200
    keep_structures: int = 100
    per_structure_adduct: int = 3
    tol_da: float = 0.02
    powers: tuple = (1.0, 2.0, 4.0)
    clean: CleanConfig = CleanConfig(min_rel_intensity=0.01, top_k_per_window=6, window_da=50, max_peaks=64)


def _vector(mz, it, precursor_mz):
    """Sparse (indices, values) over [fragment bins | neutral-loss bins], L2-normalised."""
    w = np.sqrt(it.astype(np.float64))
    frag = np.minimum((mz / BIN_DA).astype(np.int64), N_BINS - 1)
    loss_mz = precursor_mz - mz
    ok = loss_mz > 0.5
    loss = N_BINS + np.minimum((loss_mz[ok] / BIN_DA).astype(np.int64), N_BINS - 1)
    idx = np.concatenate([frag, loss])
    val = np.concatenate([w, w[ok]])
    # duplicate bins (two peaks in one bin) are summed by the CSR constructor
    norm = np.sqrt((val ** 2).sum())
    return idx, val / norm if norm > 0 else val


class AnalogIndex:
    def __init__(self, store, key_of_row, cfg: AnalogConfig = AnalogConfig(), cache: Path | None = None):
        self.store, self.cfg = store, cfg
        self._clean_cache = {}
        meta = store.meta
        if cache is not None and Path(cache).exists():
            z = np.load(cache, allow_pickle=True)
            self.rows = z["rows"]
            self.matrix = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        else:
            self.rows = self._select(meta, key_of_row)
            self.matrix = self._build(self.rows)
            if cache is not None:
                np.savez(cache, rows=self.rows, data=self.matrix.data, indices=self.matrix.indices,
                         indptr=self.matrix.indptr, shape=np.array(self.matrix.shape))
        self.key = key_of_row[self.rows]
        self.prec = meta.precursor_mz.to_numpy()[self.rows]
        self.adduct = meta.adduct.astype(str).to_numpy()[self.rows]
        self.nmass = meta.neutral_mass.to_numpy()[self.rows]

    def _select(self, meta, key_of_row):
        ok = (meta.test_adduct.to_numpy() & (np.abs(meta.ppm.to_numpy()) <= 20)
              & (meta.n_peaks.to_numpy() >= 5))
        t = pd.DataFrame({"row": np.flatnonzero(ok)})
        t["key"] = key_of_row[t.row]
        t["adduct"] = meta.adduct.to_numpy()[t.row]
        t["prio"] = meta.is_timstof.to_numpy()[t.row] * 1000 + meta.n_peaks.to_numpy()[t.row].clip(0, 999)
        t = t.sort_values("prio", ascending=False, kind="stable")
        t = t.groupby(["key", "adduct"], observed=True, sort=False).head(self.cfg.per_structure_adduct)
        return np.sort(t.row.to_numpy())

    def _clean(self, row):
        hit = self._clean_cache.get(row)
        if hit is None:
            mz, it = self.store.peaks(row)
            hit = clean_peaks_fast(mz, it, float(self.store.meta.precursor_mz.iat[row]), self.cfg.clean)
            if len(self._clean_cache) < 300_000:
                self._clean_cache[row] = hit
        return hit

    def _build(self, rows):
        indptr, indices, data = [0], [], []
        prec = self.store.meta.precursor_mz.to_numpy()
        for r in rows:
            mz, it = self._clean(r)
            i, v = _vector(mz, it, prec[r]) if mz.size else (np.empty(0, np.int64), np.empty(0))
            indices.append(i)
            data.append(v.astype(np.float32))
            indptr.append(indptr[-1] + i.size)
        self._clean_cache.clear()
        m = sp.csr_matrix((np.concatenate(data), np.concatenate(indices), np.array(indptr)),
                          shape=(len(rows), 2 * N_BINS))
        m.sum_duplicates()
        return m

    def search(self, query_spectra, query_mass, excluded_rows=None):
        """Best refined similarity per analog structure: dict key -> (sim, analog neutral mass)."""
        best = {}
        excluded = set(excluded_rows.tolist()) if excluded_rows is not None and len(excluded_rows) else set()
        mass_ok = np.abs(self.nmass - query_mass) <= self.cfg.max_shift_da
        for q in query_spectra:
            if q["mz_c"].size == 0:
                continue
            i, v = _vector(q["mz_c"], q["it_c"], q["precursor_mz"])
            # dense query vector: a CSR mat-vec needs no nnz-sized temporary (a sparse x sparse
            # product allocated ~260 MB per call and ran out of memory); same scores, ~5x faster
            dq = np.zeros(2 * N_BINS, np.float32)
            np.add.at(dq, i, v.astype(np.float32))
            scores = self.matrix @ dq
            scores[~(mass_ok & (self.adduct == q["adduct"]))] = -1
            top = np.argpartition(-scores, min(self.cfg.prefilter_top, len(scores) - 1))[: self.cfg.prefilter_top]
            for t in top:
                if scores[t] <= 0 or self.rows[t] in excluded:
                    continue
                rm, ri = self._clean(self.rows[t])
                sim = entropy_similarity(q["mz_c"], q["it_c"], rm, ri, self.cfg.tol_da,
                                         float(q["precursor_mz"] - self.prec[t]))
                k = self.key[t]
                if sim > best.get(k, (0.0, 0.0))[0]:
                    best[k] = (sim, float(self.nmass[t]))
        items = sorted(best.items(), key=lambda kv: -kv[1][0])[: self.cfg.keep_structures]
        return dict(items)


def tanimoto_matrix(fa, fb):
    """Tanimoto between packed-bit fingerprint rows: fa [n, 256] vs fb [m, 256] -> [n, m]."""
    a = np.unpackbits(fa, axis=1).astype(np.float32)
    b = np.unpackbits(fb, axis=1).astype(np.float32)
    inter = a @ b.T
    return inter / (a.sum(1)[:, None] + b.sum(1)[None, :] - inter + 1e-9)


def propagate(analogs: dict, cand_rows, cand_mass, fps, row_of, powers):
    """Per-candidate analog features. `analogs`: key -> (sim, analog neutral mass)."""
    n = len(cand_rows)
    feats = {f"an_p{int(p)}": np.zeros(n) for p in powers}
    feats.update(an_sim=np.zeros(n), an_tani=np.zeros(n), an_dmass=np.full(n, 999.0),
                 an_n_close=np.zeros(n), an_self=np.zeros(n))
    keys = [k for k in analogs if k in row_of]
    if not keys or n == 0:
        return feats
    sims = np.array([analogs[k][0] for k in keys])
    amass = np.array([analogs[k][1] for k in keys])
    T = tanimoto_matrix(fps[cand_rows], fps[[row_of[k] for k in keys]])
    for p in powers:
        feats[f"an_p{int(p)}"] = (T * sims[None, :] ** p).max(axis=1)
    best = (T * sims[None, :] ** 4).argmax(axis=1)
    feats["an_sim"] = sims[best]
    feats["an_tani"] = T[np.arange(n), best]
    feats["an_dmass"] = cand_mass - amass[best]
    feats["an_n_close"] = ((T >= 0.7) & (sims[None, :] >= 0.5)).sum(axis=1)
    feats["an_self"] = np.where(T >= 0.999, sims[None, :], 0).max(axis=1)
    return feats
