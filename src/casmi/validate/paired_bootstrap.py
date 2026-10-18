"""Minimal numerical bootstrap for CASMI ranking comparisons.

Separate from RDKit/pandas scoring so synthetic CI regression tests and
experiment comparisons need only NumPy, not the entire inference stack.
"""
from __future__ import annotations

import numpy as np


def paired_difference(ranks_a, ranks_b, n_boot=4000, seed=0):
    """Molecule-level bootstrap CI of MRR(B)-MRR(A).

    Input sequences MUST already have been aligned on verified query IDs and
    spectra. This routine does not infer or validate pairing.
    """
    a, b = np.asarray(ranks_a), np.asarray(ranks_b)
    if a.ndim != 1 or b.ndim != 1 or len(a) == 0 or a.shape != b.shape:
        raise ValueError("Paired ranks must be nonempty, one-dimensional and equal length")
    if n_boot < 100:
        raise ValueError("At least 100 bootstrap draws are required")
    to_rr = lambda r: np.where(r > 0, 1.0 / np.maximum(r, 1), 0.0)
    d = to_rr(b) - to_rr(a)
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(0, len(d), size=(n_boot, len(d)))].mean(axis=1)
    return float(d.mean()), float(np.quantile(boots, .025)), float(np.quantile(boots, .975))
