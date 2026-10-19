"""Fail-closed reconciliation before merging CASMI candidate source records.

Same first-block InChIKey indicates matched connectivity but does not justify
joining arbitrarily different neutral masses. Some source rows lack complete
identity fields; those may be excluded, but an impossible *present* mass and
a material same-key mass contradiction must never be silently accepted.

This is a source-identity quality gate, not a proof that source providers
have correct chemical annotations.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

REQUIRED = ("pool_key", "smiles", "mass")
DEFAULT_MAX_MASS_SPREAD_DA = 0.01
DEFAULT_MAX_MASS_SPREAD_PPM = 20.0


def reconcile_source_rows(
    records: pd.DataFrame, *,
    max_spread_da: float = DEFAULT_MAX_MASS_SPREAD_DA,
    max_spread_ppm: float = DEFAULT_MAX_MASS_SPREAD_PPM,
) -> tuple[pd.DataFrame, dict]:
    """Return whole, complete source rows after auditing duplicate key masses.

    The intended source order (training library before external providers)
    is preserved. Crucially, the first source's SMILES and mass must both
    exist on the SAME row before groupby.first() is used by the caller.
    """
    if not isinstance(records, pd.DataFrame):
        raise TypeError("Candidate source records must be a DataFrame")
    absent = sorted(set(REQUIRED) - set(records.columns))
    if absent:
        raise ValueError("Candidate sources missing required fields: " + ", ".join(absent))
    if (not np.isfinite(max_spread_da) or max_spread_da < 0
            or not np.isfinite(max_spread_ppm) or max_spread_ppm < 0):
        raise ValueError("Mass disagreement tolerances must be finite nonnegative values")
    complete = records.dropna(subset=list(REQUIRED)).copy()
    omitted = len(records) - len(complete)
    if complete.empty:
        raise ValueError("No complete candidate source identity rows available")
    # A source key is the standard InChIKey connectivity block. Reject
    # malformed *present* keys instead of merging them as new structures.
    if (not complete["pool_key"].map(lambda value:
                                    isinstance(value, str) and len(value) == 14
                                    and value.isascii() and value.isalpha()
                                    and value.isupper()).all()):
        raise ValueError("Malformed candidate source InChIKey14")
    if not complete["smiles"].map(lambda value:
                                   isinstance(value, str) and bool(value.strip())).all():
        raise ValueError("Empty or malformed candidate source SMILES")
    try:
        mass = pd.to_numeric(complete["mass"], errors="raise").to_numpy(dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Non-numeric candidate source mass") from exc
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0):
        raise ValueError("Candidate source mass must be finite and positive")
    complete.loc[:, "mass"] = mass

    ranges = complete.groupby("pool_key", sort=False)["mass"].agg(["min", "max"])
    discrepancies = (ranges["max"] - ranges["min"]).to_numpy(dtype=float)
    thresholds = np.maximum(
        float(max_spread_da),
        ranges["min"].to_numpy(dtype=float) * float(max_spread_ppm) * 1e-6,
    )
    conflicts = int(np.count_nonzero(discrepancies > thresholds))
    if conflicts:
        # Aggregate count only: no licensed source identifiers, chemical
        # structures, or per-molecule annotations escape into log messages.
        raise ValueError(
            f"Conflicting neutral mass for {conflicts} shared InChIKey14 groups; "
            "quarantine or correct providers before building a new candidate pool"
        )
    summary = {
        "status": "SOURCE_ROWS_RECONCILED",
        "input_rows": int(len(records)),
        "complete_rows": int(len(complete)),
        "incomplete_rows_excluded": int(omitted),
        "distinct_pool_keys": int(len(ranges)),
        "shared_pool_keys": int((complete.groupby("pool_key").size() > 1).sum()),
        "mass_conflict_groups": 0,
        "max_source_mass_spread_da": float(discrepancies.max(initial=0.0)),
        "mass_threshold_da": float(max_spread_da),
        "mass_threshold_ppm": float(max_spread_ppm),
        "interpretation": "Pre-aggregation identity consistency only; not RDKit verification",
    }
    return complete, summary
