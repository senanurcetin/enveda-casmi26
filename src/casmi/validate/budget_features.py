"""Recompute list-relative features after selecting a declared candidate budget."""
from __future__ import annotations

import numpy as np
import pandas as pd

from casmi.rank.features import add_context


def budget_features(frame, policy):
    if policy not in ("legacy", "expanded"):
        raise ValueError("Unknown declared candidate-budget policy")
    if frame.empty:
        return frame.copy()
    if "audit_legacy" not in frame or not frame.audit_legacy.isin([0, 1]).all():
        raise ValueError("Candidate origin audit is missing or invalid")
    selected = frame if policy == "expanded" else frame[frame.audit_legacy == 1]
    selected = selected.drop(columns=[c for c in selected if c.startswith("audit_")]).reset_index(drop=True).copy()
    if selected.empty:
        return selected
    selected["n_cand"] = len(selected)
    selected = add_context(selected)
    for column in ("d3_ll", "d3_cos", "d4_int", "d4_pk", "d4x_int", "d4x_pk"):
        if column not in selected:
            continue
        values = selected[column].to_numpy(dtype=float)
        sd = values.std()
        selected[column + "_rank"] = pd.Series(values).rank(ascending=False, method="min").to_numpy()
        selected[column + "_gap"] = values.max() - values
        selected[column + "_z"] = (values - values.mean()) / sd if sd > 0 else np.zeros_like(values)
    return selected
