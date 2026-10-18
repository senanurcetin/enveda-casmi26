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
from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib
import json
import zipfile

import numpy as np
import pandas as pd
import scipy.sparse as sp

from casmi.spectra.clean import CleanConfig, clean_peaks_fast
from casmi.channels.prefilter import eligible_top_k
from casmi.channels.reference_selection import select_reference_rows, SELECTION_VERSION
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
    reference_policy: str = "legacy"  # opt-in: source_balanced (requires new index cache)
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
        expected = None
        stamp = None
        if cfg.reference_policy != "legacy":
            # A policy-specific cache never overwrites the frozen leaderboard
            # baseline's analog_index.npz. Nonlegacy loads are fail-closed.
            expected = self._select(meta, key_of_row)
            stamp = self._cache_stamp(store, key_of_row, expected)
            if cache is not None:
                path = Path(cache)
                cache = path.with_name(
                    f"{path.stem}_{cfg.reference_policy}_v1{path.suffix}"
                )
        loaded = False
        if cache is not None and Path(cache).exists():
            try:
                with np.load(cache, allow_pickle=False) as z:
                    if expected is not None and (
                        str(z["selection_version"].item()) != SELECTION_VERSION
                        or str(z["selection_digest"].item()) != stamp
                        or not np.array_equal(z["rows"], expected)
                    ):
                        raise ValueError("Cached D2 selection provenance mismatch")
                    cached_rows = z["rows"]
                    cached_matrix = sp.csr_matrix(
                        (z["data"], z["indices"], z["indptr"]),
                        shape=tuple(z["shape"]),
                    )
                    if cached_matrix.shape != (len(cached_rows), 2 * N_BINS):
                        raise ValueError("Cached D2 matrix shape is invalid")
                    if len(cached_rows) and (cached_rows.min() < 0 or cached_rows.max() >= len(meta)):
                        raise ValueError("Cached D2 reference row index out of range")
                    self.rows, self.matrix = cached_rows, cached_matrix
                    loaded = True
            except (KeyError, OSError, ValueError, zipfile.BadZipFile):
                if expected is None:
                    raise  # fail closed on malformed frozen legacy caches
        if not loaded:
            self.rows = expected if expected is not None else self._select(meta, key_of_row)
            self.matrix = self._build(self.rows)
            if cache is not None:
                output = dict(rows=self.rows, data=self.matrix.data,
                              indices=self.matrix.indices, indptr=self.matrix.indptr,
                              shape=np.array(self.matrix.shape))
                if expected is not None:
                    output.update(selection_version=np.array(SELECTION_VERSION),
                                  selection_digest=np.array(stamp))
                np.savez(cache, **output)
        self.key = key_of_row[self.rows]
        self.prec = meta.precursor_mz.to_numpy()[self.rows]
        self.adduct = meta.adduct.astype(str).to_numpy()[self.rows]
        self.nmass = meta.neutral_mass.to_numpy()[self.rows]

    def _select(self, meta, key_of_row):
        return select_reference_rows(
            meta, key_of_row,
            per_structure_adduct=self.cfg.per_structure_adduct,
            policy=self.cfg.reference_policy,
        )

    def _cache_stamp(self, store, keys, rows):
        """Hash selected metadata, source-array file stamps and index config.

        The selected rows are recomputed on every nonlegacy cache load, so
        updated eligibility or changes in available libraries invalidate the
        cached index even if the number of spectra is unchanged.
        """
        digest = hashlib.sha256()
        digest.update(SELECTION_VERSION.encode())
        digest.update(json.dumps(asdict(self.cfg), sort_keys=True).encode())
        digest.update(np.asarray(rows, dtype=np.int64).tobytes())
        selected = store.meta.iloc[rows][[
            "ingest_lib", "adduct", "precursor_mz", "neutral_mass",
            "ppm", "n_peaks", "is_timstof", "test_adduct",
        ]]
        digest.update(pd.util.hash_pandas_object(selected, index=False).values.tobytes())
        digest.update(pd.util.hash_pandas_object(
            pd.Series(np.asarray(keys)[rows], dtype="string"), index=False,
        ).values.tobytes())
        for field in ("mz", "it", "offsets"):
            ref = getattr(store, field, None)
            filename = getattr(ref, "filename", None)
            if filename is None:
                continue
            stat = Path(filename).stat()
            digest.update(
                f"{field}:{stat.st_size}:{stat.st_mtime_ns}".encode()
            )
        return digest.hexdigest()

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
        # A structure-level knockout must be applied BEFORE selecting the top
        # analog hits. Previously self-references displaced allowed hits in the
        # shortlist and were discarded only afterwards.
        excluded = (np.asarray(excluded_rows, dtype=np.int64)
                    if excluded_rows is not None else np.empty(0, np.int64))
        not_excluded = (~np.isin(self.rows, excluded)
                        if excluded.size else np.ones(len(self.rows), dtype=bool))
        mass_ok = np.isfinite(self.nmass) & np.isfinite(query_mass) & (
            np.abs(self.nmass - query_mass) <= self.cfg.max_shift_da
        )
        for q in query_spectra:
            if q["mz_c"].size == 0:
                continue
            i, v = _vector(q["mz_c"], q["it_c"], q["precursor_mz"])
            # Dense CSR mat-vec avoids a large sparse x sparse temporary.
            dq = np.zeros(2 * N_BINS, np.float32)
            np.add.at(dq, i, v.astype(np.float32))
            scores = self.matrix @ dq
            allowed = not_excluded & mass_ok & (self.adduct == q["adduct"])
            top = eligible_top_k(scores, allowed, self.cfg.prefilter_top)
            for t in top:
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
