"""Shared numerical contracts for fixed MRR@25 research summaries."""
from __future__ import annotations

import numpy as np


def validated_ranks(values):
    ranks = np.asarray(values)
    if ranks.ndim != 1 or not len(ranks):
        raise ValueError("Ranks must be a nonempty one-dimensional integer sequence")
    if not np.issubdtype(ranks.dtype, np.integer):
        raise ValueError("Ranks must be integers from 0 to 25; 0 denotes a miss")
    if np.any((ranks < 0) | (ranks > 25)):
        raise ValueError("Ranks must lie between 0 and 25; out-of-list hits are encoded as 0")
    return ranks


def validated_bootstrap(n_boot, seed):
    if not isinstance(n_boot, (int, np.integer)) or isinstance(n_boot, (bool, np.bool_)) or n_boot < 100:
        raise ValueError("At least 100 integer bootstrap draws are required")
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)) or seed < 0:
        raise ValueError("Bootstrap seed must be a nonnegative integer")
