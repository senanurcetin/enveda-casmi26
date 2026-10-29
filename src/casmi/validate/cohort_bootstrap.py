"""Joint fixed-family uncertainty with independent resampling across source cohorts."""
from __future__ import annotations

import numpy as np

from .rank_inputs import validated_bootstrap


def joint_cohort_intervals(effects_by_cohort, *, n_boot=4000, seed=0):
    validated_bootstrap(n_boot, seed)
    if not effects_by_cohort:
        raise ValueError("At least one declared source cohort required")
    groups, names_seen = {}, set()
    for cohort, effects in sorted(effects_by_cohort.items()):
        if not effects or set(effects) & names_seen:
            raise ValueError("Empty or duplicate secondary contrast family")
        names_seen.update(effects)
        arrays = [np.asarray(value, dtype=float) for value in effects.values()]
        if any(a.ndim != 1 or len(a) < 2 or not np.isfinite(a).all() or np.any(np.abs(a) > 1) for a in arrays):
            raise ValueError("Finite bounded paired effects from at least two molecules required")
        if len({len(a) for a in arrays}) != 1:
            raise ValueError("Every contrast in one cohort must share paired molecules")
        groups[cohort] = (list(effects), np.column_stack(arrays))
    maximum = np.zeros(n_boot)
    outcomes, counts = {}, {}
    for child_seed, (cohort, (names, values)) in zip(np.random.SeedSequence(seed).spawn(len(groups)), groups.items()):
        mean = values.mean(axis=0)
        counts[cohort] = len(values)
        rng = np.random.default_rng(child_seed)
        for start in range(0, n_boot, 100):
            ids = rng.integers(len(values), size=(min(100, n_boot - start), len(values)))
            deviations = np.max(np.abs(values[ids].mean(axis=1) - mean), axis=1)
            maximum[start:start + len(ids)] = np.maximum(maximum[start:start + len(ids)], deviations)
        outcomes.update({name: {"cohort": cohort, "delta": float(value)} for name, value in zip(names, mean)})
    radius = float(np.quantile(maximum, .95))
    for value in outcomes.values():
        value.update(ci_low=value["delta"] - radius, ci_high=value["delta"] + radius)
    return {"method": "Independent molecule resampling across cohorts; paired within cohort; fixed-family max-deviation bootstrap",
            "n_bootstrap": n_boot, "seed": seed, "family_size": len(outcomes), "cohort_counts": counts,
            "joint_radius": radius, "outcomes": outcomes}
