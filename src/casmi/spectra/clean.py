"""Peak cleaning. One function is used for library spectra (offline) and query spectra (online),
so the two can never drift apart.

`clean_peaks` is the readable NumPy reference; `clean_peaks_fast` is the numba version used on
millions of spectra. tests/test_spectra.py asserts they agree on random spectra.
"""
from dataclasses import dataclass

import numba as nb
import numpy as np


@dataclass(frozen=True)
class CleanConfig:
    max_above_precursor: float = 2.0   # Da; a fragment cannot outweigh its singly charged precursor
    min_rel_intensity: float = 0.001   # relative to base peak
    window_da: float = 50.0            # windowed top-k keeps low-mass fragments a global top-N drops
    top_k_per_window: int = 0          # 0 disables windowed filtering
    max_peaks: int = 150               # global cap after the filters above
    merge_tol_da: float = 0.005        # peaks closer than this are merged (centroiding artefacts)


def clean_peaks(mz, intensity, precursor_mz, cfg: CleanConfig = CleanConfig()):
    """Return (mz, intensity) sorted by m/z, float32, intensity normalised to base peak = 1."""
    mz = np.asarray(mz, dtype=np.float64)
    it = np.asarray(intensity, dtype=np.float64)
    keep = np.isfinite(mz) & np.isfinite(it) & (it > 0)
    if precursor_mz is not None and np.isfinite(precursor_mz):
        keep &= mz <= precursor_mz + cfg.max_above_precursor
    mz, it = mz[keep], it[keep]
    if mz.size == 0:
        return np.empty(0, np.float32), np.empty(0, np.float32)

    order = np.argsort(mz, kind="stable")
    mz, it = mz[order], it[order]
    if cfg.merge_tol_da > 0 and mz.size > 1:
        # merge runs of peaks closer than merge_tol_da: intensity-weighted m/z, summed intensity
        group = np.concatenate([[0], np.cumsum(np.diff(mz) > cfg.merge_tol_da)])
        if group[-1] + 1 < mz.size:
            w = np.bincount(group, weights=it)
            mz = np.bincount(group, weights=mz * it) / w
            it = w

    it = it / it.max()
    keep = it >= cfg.min_rel_intensity
    mz, it = mz[keep], it[keep]

    if cfg.top_k_per_window and mz.size > cfg.top_k_per_window:
        win = np.floor(mz / cfg.window_da).astype(np.int64)
        # rank peaks within each window by intensity (descending)
        order = np.lexsort((-it, win))
        rank = np.empty(mz.size, np.int64)
        starts = np.r_[0, np.flatnonzero(np.diff(win[order])) + 1]
        run = np.arange(mz.size) - np.repeat(starts, np.diff(np.r_[starts, mz.size]))
        rank[order] = run
        keep = rank < cfg.top_k_per_window
        mz, it = mz[keep], it[keep]

    if mz.size > cfg.max_peaks:
        top = np.argpartition(-it, cfg.max_peaks - 1)[:cfg.max_peaks]
        top.sort()
        mz, it = mz[top], it[top]
    return mz.astype(np.float32), it.astype(np.float32)


@nb.njit(cache=True)
def _clean_nb(mz, it, precursor_mz, max_above, min_rel, window_da, top_k, max_peaks, merge_tol):
    n = mz.size
    keep = np.empty(n, np.bool_)
    lim = precursor_mz + max_above if np.isfinite(precursor_mz) else np.inf
    m = 0
    for i in range(n):
        ok = np.isfinite(mz[i]) and np.isfinite(it[i]) and it[i] > 0 and mz[i] <= lim
        keep[i] = ok
        m += ok
    mz, it = mz[keep], it[keep]
    if m == 0:
        return np.empty(0, np.float32), np.empty(0, np.float32)
    order = np.argsort(mz, kind="mergesort")
    mz, it = mz[order], it[order]

    # chain-merge consecutive peaks closer than merge_tol
    out_mz = np.empty(m, np.float64)
    out_it = np.empty(m, np.float64)
    g = 0
    wsum, isum = mz[0] * it[0], it[0]
    for i in range(1, m):
        if merge_tol > 0 and mz[i] - mz[i - 1] <= merge_tol:
            wsum += mz[i] * it[i]
            isum += it[i]
        else:
            out_mz[g], out_it[g] = wsum / isum, isum
            g += 1
            wsum, isum = mz[i] * it[i], it[i]
    out_mz[g], out_it[g] = wsum / isum, isum
    g += 1
    mz, it = out_mz[:g], out_it[:g]

    it = it / it.max()
    sel = it >= min_rel
    mz, it = mz[sel], it[sel]

    if top_k > 0 and mz.size > top_k:
        win = np.floor(mz / window_da)
        keepw = np.zeros(mz.size, np.bool_)
        s = 0
        while s < mz.size:
            e = s
            while e < mz.size and win[e] == win[s]:
                e += 1
            seg = np.argsort(-it[s:e], kind="mergesort")
            for j in range(min(top_k, e - s)):
                keepw[s + seg[j]] = True
            s = e
        mz, it = mz[keepw], it[keepw]

    if mz.size > max_peaks:
        top = np.sort(np.argsort(-it, kind="mergesort")[:max_peaks])
        mz, it = mz[top], it[top]
    return mz.astype(np.float32), it.astype(np.float32)


def clean_peaks_fast(mz, intensity, precursor_mz, cfg: CleanConfig = CleanConfig()):
    """numba equivalent of `clean_peaks` (~50x faster)."""
    p = np.nan if precursor_mz is None else float(precursor_mz)
    return _clean_nb(np.asarray(mz, np.float64), np.asarray(intensity, np.float64), p,
                     cfg.max_above_precursor, cfg.min_rel_intensity, cfg.window_da,
                     cfg.top_k_per_window, cfg.max_peaks, cfg.merge_tol_da)
