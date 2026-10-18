"""Numerically safe, deterministic top-k prefilter for MS/MS analog retrieval.

All eligibility restrictions, especially structure-level knockout rows, must
be enforced *before* top-k selection. Otherwise withheld/self-references can
consume the entire shortlist and silently suppress genuine analogs.

This module requires only NumPy, so offline numerical regression tests need
neither the full spectrum store nor GPU/chemistry dependencies.
"""
from __future__ import annotations

import numpy as np


def eligible_top_k(scores, allowed, k: int) -> np.ndarray:
    """Indices of up to k highest-scoring eligible, positive, finite entries.

    Ties are resolved by ascending original index (stable, deterministic);
    inputs are never mutated. The order is descending score, then index.
    """
    scores = np.asarray(scores)
    allowed = np.asarray(allowed)
    if scores.ndim != 1 or allowed.ndim != 1 or scores.shape != allowed.shape:
        raise ValueError("scores and allowed must be aligned 1D vectors")
    if not np.issubdtype(scores.dtype, np.number):
        raise ValueError("scores must be numeric")
    if allowed.dtype != np.dtype(bool):
        raise ValueError("allowed must be a boolean vector")
    if type(k) is not int or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if k == 0 or scores.size == 0:
        return np.empty(0, dtype=np.int64)
    eligible = np.flatnonzero(allowed & np.isfinite(scores) & (scores > 0))
    if eligible.size <= k:
        return eligible[np.lexsort((eligible, -scores[eligible]))]
    # Use the score of the kth-best candidate as a cutoff, then break ties
    # deterministically by index without sorting the whole reference library.
    local_values = scores[eligible]
    cutoff = np.partition(local_values, local_values.size - k)[local_values.size - k]
    better = eligible[local_values > cutoff]
    tied = eligible[local_values == cutoff]
    chosen = np.concatenate((better, tied[: k - better.size]))
    return chosen[np.lexsort((chosen, -scores[chosen]))]
