"""Simultaneous paired-bootstrap uncertainty for CASMI novelty mixtures.

The old comparator calculated 12 pointwise 95% intervals and inspected
all of them. They are not a joint 95% confidence statement.
Each bootstrap draw here reuses the same molecules across C1/C2/C3;
a max-deviation statistic gives approximately simultaneous intervals
over the predeclared novelty-share grid.

These bounds capture validation sampling variation, not Kaggle domain shift,
hyperparameter search multiplicity, model bias, or data leakage.
"""
from __future__ import annotations

import numpy as np

DEFAULT_GRID = tuple(
    (f1, f3) for f1 in (0.15, 0.20, 0.25, 0.30)
    for f3 in (0.2, 0.3, 0.4)
)


def _ranks_to_rr(values, cls, side):
    ranks = np.asarray(values)
    if ranks.ndim != 1 or not len(ranks):
        raise ValueError(f"{side} {cls}: ranks must be nonempty 1D")
    if not np.issubdtype(ranks.dtype, np.integer):
        raise ValueError(f"{side} {cls}: ranks must be integers")
    if np.any((ranks < 0) | (ranks > 25)):
        raise ValueError(f"{side} {cls}: invalid rank (valid: 0 to 25)")
    return np.where(ranks > 0, 1.0 / np.maximum(ranks, 1), 0.0)


def paired_novelty_deltas(a: dict, b: dict) -> np.ndarray:
    """Molecule x class effect matrix, with preverified aligned query ordering."""
    columns = []
    n = None
    for cls in ("C1", "C2", "C3"):
        if cls not in a or cls not in b:
            raise ValueError(f"Missing novelty class {cls}")
        x, y = _ranks_to_rr(a[cls], cls, "A"), _ranks_to_rr(b[cls], cls, "B")
        if x.shape != y.shape or (n is not None and len(x) != n):
            raise ValueError("Unpaired class ranks, unequal lengths")
        n = len(x)
        columns.append(y - x)
    return np.stack(columns, axis=1)


def compare_mixtures(
    a: dict,
    b: dict,
    *,
    grid=DEFAULT_GRID,
    n_boot: int = 4000,
    seed: int = 0,
    confidence: float = 0.95,
) -> dict:
    """Approximate simultaneous CI over a fixed mixture grid.

    Must follow align_runs; this function is not an identity check.
    Bootstrap rows remain paired across all three novelty classes.
    """
    if type(n_boot) is not int or n_boot < 100:
        raise ValueError("Use at least 100 bootstrap resamples")
    if not 0 < confidence < 1:
        raise ValueError("Confidence must lie strictly between 0 and 1")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Seed must be a nonnegative integer")
    shares = np.asarray(grid, dtype=float)
    if shares.ndim != 2 or shares.shape[1] != 2 or len(shares) == 0:
        raise ValueError("Grid must contain (f1, f3) pairs")
    if not np.isfinite(shares).all():
        raise ValueError("Mixture shares must be finite")
    f1, f3 = shares[:, 0], shares[:, 1]
    f2 = 1.0 - f1 - f3
    if np.any(f1 < 0) or np.any(f3 < 0) or np.any(f2 < -1e-12):
        raise ValueError("Mixture fractions must be nonnegative and sum to 1")
    f2 = np.maximum(f2, 0.0)
    weights = np.stack([f1, f2, f3], axis=1)

    effects = paired_novelty_deltas(a, b)
    mean = effects.mean(axis=0) @ weights.T
    rng = np.random.default_rng(int(seed))
    indices = rng.integers(0, effects.shape[0], size=(n_boot, effects.shape[0]))
    sampled = effects[indices].mean(axis=1) @ weights.T
    max_abs_deviation = np.max(np.abs(sampled - mean), axis=1)
    radius = float(np.quantile(max_abs_deviation, confidence))
    lower, upper = mean - radius, mean + radius
    verdict = ("ACCEPT" if np.all(lower > 0) else
               "REJECT" if np.all(upper < 0) else "NOT PROVEN")
    scenarios = [
        {"f1": float(f1[i]), "f2": float(f2[i]), "f3": float(f3[i]),
         "delta_mrr": float(mean[i]), "ci_low": float(lower[i]),
         "ci_high": float(upper[i])}
        for i in range(len(mean))
    ]
    return {
        "method": "paired bootstrap; maximum absolute deviation over fixed mixture grid",
        "confidence": confidence,
        "familywise": "approximate for a fixed prespecified grid; not post-selection corrected",
        "n_molecules": int(len(effects)),
        "n_boot": n_boot,
        "seed": int(seed),
        "half_width": radius,
        "verdict": verdict,
        "scenarios": scenarios,
    }
