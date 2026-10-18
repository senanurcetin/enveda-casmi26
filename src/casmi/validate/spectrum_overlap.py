"""Exact-peak overlap helpers for CASMI validation leakage checks.

A precursor mass, number of peaks and first m/z are a *bucket*, not evidence
of a duplicate spectrum. Only full arrays, including intensities, can confirm
a same-spectrum match. This module deliberately has no numerical dependencies.
"""
from __future__ import annotations

import math
from collections.abc import Sequence


def coarse_spectrum_bucket(precursor_mz: float, ms2_mzs: Sequence) -> tuple | None:
    """Cheap candidate retrieval; NEVER use the bucket as a match on its own."""
    if ms2_mzs is None or len(ms2_mzs) == 0:
        return None
    try:
        precursor, first = float(precursor_mz), float(ms2_mzs[0])
    except (TypeError, ValueError):
        return None
    if not math.isfinite(precursor) or not math.isfinite(first):
        return None
    return round(precursor, 4), len(ms2_mzs), round(first, 4)


def full_spectrum_match(
    precursor_a, mz_a: Sequence, intensities_a: Sequence,
    precursor_b, mz_b: Sequence, intensities_b: Sequence,
) -> bool:
    """Require full precursor, fragment m/z and normalized intensity equality.

    Converts array scalar dtypes to Python float without rounding. An equal
    coarse bucket with shifted *any* peak or intensity is not a verified match.
    """
    if mz_a is None or mz_b is None or intensities_a is None or intensities_b is None:
        return False
    if len(mz_a) == 0 or len(mz_b) == 0:
        return False
    if len(mz_a) != len(mz_b) or len(intensities_a) != len(intensities_b):
        return False
    if len(mz_a) != len(intensities_a):
        return False
    try:
        values_a = [float(precursor_a)] + [float(x) for x in mz_a] + [float(x) for x in intensities_a]
        values_b = [float(precursor_b)] + [float(x) for x in mz_b] + [float(x) for x in intensities_b]
    except (TypeError, ValueError):
        return False
    if not all(math.isfinite(x) for x in values_a + values_b):
        return False
    return values_a == values_b
