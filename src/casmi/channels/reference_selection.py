"""Reference sampling for the CASMI D2 analog index.

The historical policy selects up to N spectra per (standard key, adduct),
preferring timsTOF rows. This can eliminate every independent reference
acquired in another library before a cross-library C1 knockout is applied.

The opt-in source_balanced policy retains the best row, reserves one
non-timsTOF source (if the best is timsTOF and N >= 2), then chooses distinct
libraries before filling unused slots by quality. It controls reference
coverage; it does NOT guarantee improved MRR or metric-key equivalence.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

POLICIES = ("legacy", "source_balanced")
SELECTION_VERSION = "casmi-d2-reference-selection-v1"


def select_reference_rows(meta: pd.DataFrame, key_of_row, *,
                          per_structure_adduct: int = 3,
                          policy: str = "legacy") -> np.ndarray:
    """Return sorted, unique spectrum row IDs for building the D2 index.

    The legacy branch intentionally reproduces AnalogIndex._select,
    including its global timsTOF and peak-count priority and stable ties.
    Source balancing does not raise the existing per-(key, adduct) cap.
    """
    if policy not in POLICIES:
        raise ValueError(f"Unknown selection policy {policy!r}")
    if type(per_structure_adduct) is not int or per_structure_adduct < 1:
        raise ValueError("per_structure_adduct must be a positive integer")
    keys = np.asarray(key_of_row)
    if len(meta) != len(keys):
        raise ValueError("Reference key and metadata lengths differ")
    needed = {"test_adduct", "ppm", "n_peaks", "is_timstof", "ingest_lib", "adduct"}
    if not needed.issubset(meta.columns):
        raise ValueError(f"Missing index metadata: {sorted(needed - set(meta.columns))}")
    ok = (meta.test_adduct.to_numpy(dtype=bool) &
          (np.abs(meta.ppm.to_numpy(dtype=float)) <= 20) &
          (meta.n_peaks.to_numpy(dtype=float) >= 5))
    rows = np.flatnonzero(ok)
    if len(rows) == 0:
        return np.empty(0, dtype=np.int64)
    t = pd.DataFrame({"row": rows})
    t["key"] = keys[rows]
    t["adduct"] = meta.adduct.to_numpy()[rows]
    t["prio"] = (
        meta.is_timstof.to_numpy(dtype=bool)[rows].astype(np.int64) * 1000
        + np.clip(meta.n_peaks.to_numpy(dtype=float)[rows], 0, 999)
    )
    t = t.sort_values("prio", ascending=False, kind="stable")
    if policy == "legacy":
        # Exact historical selection path.
        t = t.groupby(["key", "adduct"], observed=True, sort=False).head(per_structure_adduct)
        return np.sort(t.row.to_numpy(dtype=np.int64))

    t["lib"] = meta.ingest_lib.astype(str).to_numpy()[t.row.to_numpy()]
    t["timstof"] = meta.is_timstof.to_numpy(dtype=bool)[t.row.to_numpy()]
    selected = []
    for _, group in t.groupby(["key", "adduct"], observed=True, sort=False):
        group_rows = group.row.to_numpy(dtype=np.int64)
        if len(group_rows) <= per_structure_adduct:
            selected.extend(group_rows)
            continue
        libs = group.lib.to_numpy()
        timstof = group.timstof.to_numpy(dtype=bool)
        chosen = [0]  # Never displace the highest-quality reference.
        used_libs = {libs[0]}
        if per_structure_adduct > 1 and timstof[0]:
            for i in range(1, len(group_rows)):
                if not timstof[i] and libs[i] not in used_libs:
                    chosen.append(i)
                    used_libs.add(libs[i])
                    break
        for i in range(1, len(group_rows)):
            if len(chosen) == per_structure_adduct:
                break
            if libs[i] not in used_libs:
                chosen.append(i)
                used_libs.add(libs[i])
        if len(chosen) < per_structure_adduct:
            known = set(chosen)
            for i in range(1, len(group_rows)):
                if len(chosen) == per_structure_adduct:
                    break
                if i not in known:
                    chosen.append(i)
        selected.extend(group_rows[chosen])
    return np.sort(np.asarray(selected, dtype=np.int64))
