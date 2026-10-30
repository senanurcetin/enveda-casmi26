"""Bounded optimal one-to-one cosine baseline; no learned parameters.

The assignment maximizes intensity-product weight, unlike proximity-greedy
matching. SciPy is imported only when the optional baseline is requested.
"""
from __future__ import annotations

import numpy as np


def optimal_cosine(mz_a, it_a, mz_b, it_b, tol_da=0.02, shift=0.0,
                   *, intensity_power=1.0, max_peaks=256):
    """Cosine with exact peak assignment and optional precursor mass shift.

    All peaks enter the denominator, including unmatched peaks. Direct and
    shifted edges share the same assignment, preventing double counting.
    Unsorted input is supported. Invalid inputs fail rather than become scores.
    """
    if (not np.isfinite(tol_da) or tol_da < 0 or not np.isfinite(shift)
            or not np.isfinite(intensity_power) or intensity_power <= 0):
        raise ValueError('Finite nonnegative tolerance, finite shift and positive power required')
    if type(max_peaks) is not int or not 1 <= max_peaks <= 1024:
        raise ValueError('max_peaks must be an integer from 1 through 1024')
    arrays = [np.asarray(v, dtype=np.float64) for v in (mz_a, it_a, mz_b, it_b)]
    for mz, intensity in (arrays[:2], arrays[2:]):
        if mz.ndim != 1 or intensity.ndim != 1 or mz.shape != intensity.shape:
            raise ValueError('Aligned one-dimensional peak arrays required')
        if len(mz) > max_peaks:
            raise ValueError('Peak budget exceeded; clean spectra explicitly before scoring')
        if not np.isfinite(mz).all() or not np.isfinite(intensity).all() or (mz < 0).any() or (intensity < 0).any():
            raise ValueError('Finite nonnegative peak arrays required')
    ma, a, mb, b = arrays
    if not len(a) or not len(b) or a.max() == 0 or b.max() == 0:
        return 0.0
    # Scale before powers/norms so large finite acquisition units cannot overflow.
    a, b = (a / a.max()) ** intensity_power, (b / b.max()) ** intensity_power
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    distance = mb[None, :] - ma[:, None]
    allowed = (np.abs(distance) <= tol_da) | (np.abs(distance - shift) <= tol_da)
    weights = np.outer(a, b) * allowed
    from scipy.optimize import linear_sum_assignment
    rows, cols = linear_sum_assignment(weights, maximize=True)
    return float(np.clip(weights[rows, cols].sum(), 0.0, 1.0))
