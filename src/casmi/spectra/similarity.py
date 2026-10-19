"""Spectral similarities, compiled with numba.

Entropy similarity (Li et al., Nature Methods 2021) with intensity weighting:
    sim = sum over matched peak pairs [ (a+b)ln(a+b) - a ln a - b ln b ] / ln 4
with each spectrum's intensities normalised to sum 1. This closed form follows from
1 - (2 S_AB - S_A - S_B) / ln 4 and needs only the matched pairs.

`shift` enables the mass-shifted ("modified") variant used for analog search: a peak may match
either at the same m/z or at m/z + shift, where shift = precursor difference between the spectra.
"""
import numba as nb
import numpy as np

LN4 = np.log(4.0)


@nb.njit(cache=True)
def _entropy(p):
    s = 0.0
    for x in p:
        if x > 0:
            s -= x * np.log(x)
    return s


@nb.njit(cache=True)
def _valid_spectrum(mz, it):
    """Reject invalid numerical peak arrays before any logarithm or square root."""
    if mz.ndim != 1 or it.ndim != 1 or mz.size != it.size:
        raise ValueError("Spectrum m/z and intensity must be aligned 1D arrays")
    for i in range(mz.size):
        if (not np.isfinite(mz[i]) or mz[i] < 0.0
                or not np.isfinite(it[i]) or it[i] < 0.0):
            return False
    return True


@nb.njit(cache=True)
def weight_intensity(it):
    """Normalise positive intensity; return zeros for zero/invalid total."""
    total = it.sum()
    if not np.isfinite(total) or total <= 0.0:
        return np.zeros(it.size, dtype=np.float64)
    p = it / total
    s = _entropy(p)
    if s < 3.0:
        w = 0.25 + 0.25 * s
        p = p ** w
        p = p / p.sum()
    return p


@nb.njit(cache=True)
def _match(mz_a, mz_b, tol_da, shift):
    """Greedy one-to-one matching of peaks within tol_da, allowing an optional mass shift.
    Returns index pairs sorted by m/z proximity (closest pairs claimed first)."""
    na, nb_ = mz_a.size, mz_b.size
    cand_i, cand_j, cand_d = [], [], []
    for s_idx in range(2 if shift != 0.0 else 1):
        off = shift if s_idx == 1 else 0.0
        j0 = 0
        for i in range(na):
            lo = mz_a[i] + off - tol_da
            while j0 < nb_ and mz_b[j0] < lo:
                j0 += 1
            j = j0
            while j < nb_ and mz_b[j] <= mz_a[i] + off + tol_da:
                cand_i.append(i)
                cand_j.append(j)
                cand_d.append(abs(mz_b[j] - mz_a[i] - off))
                j += 1
    n = len(cand_i)
    order = np.argsort(np.array(cand_d)) if n else np.empty(0, np.int64)
    used_a = np.zeros(na, np.bool_)
    used_b = np.zeros(nb_, np.bool_)
    pi, pj = [], []
    for k in order:
        i, j = cand_i[k], cand_j[k]
        if not used_a[i] and not used_b[j]:
            used_a[i] = True
            used_b[j] = True
            pi.append(i)
            pj.append(j)
    return pi, pj


@nb.njit(cache=True)
def entropy_similarity(mz_a, it_a, mz_b, it_b, tol_da=0.02, shift=0.0):
    if not _valid_spectrum(mz_a, it_a) or not _valid_spectrum(mz_b, it_b):
        return 0.0
    if mz_a.size == 0 or mz_b.size == 0:
        return 0.0
    pa = weight_intensity(it_a.astype(np.float64))
    pb = weight_intensity(it_b.astype(np.float64))
    if pa.sum() <= 0.0 or pb.sum() <= 0.0:
        return 0.0
    pi, pj = _match(mz_a.astype(np.float64), mz_b.astype(np.float64), tol_da, shift)
    s = 0.0
    for k in range(len(pi)):
        a, b = pa[pi[k]], pb[pj[k]]
        if a <= 0.0 or b <= 0.0:
            continue  # 0*log(0) is mathematically zero, not NaN.
        s += (a + b) * np.log(a + b) - a * np.log(a) - b * np.log(b)
    return s / LN4


@nb.njit(cache=True)
def cosine_similarity(mz_a, it_a, mz_b, it_b, tol_da=0.02, shift=0.0):
    """Greedy (modified, if shift != 0) cosine on sqrt intensities."""
    if not _valid_spectrum(mz_a, it_a) or not _valid_spectrum(mz_b, it_b):
        return 0.0
    if mz_a.size == 0 or mz_b.size == 0:
        return 0.0
    a = np.sqrt(it_a.astype(np.float64))
    b = np.sqrt(it_b.astype(np.float64))
    pi, pj = _match(mz_a.astype(np.float64), mz_b.astype(np.float64), tol_da, shift)
    s = 0.0
    for k in range(len(pi)):
        s += a[pi[k]] * b[pj[k]]
    norm = np.sqrt((a * a).sum()) * np.sqrt((b * b).sum())
    return s / norm if norm > 0.0 and np.isfinite(norm) else 0.0
