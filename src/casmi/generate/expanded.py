"""Explicitly budgeted proposal expansion preserving historical candidates/features."""
from __future__ import annotations

from .transforms import propose_agg


def expanded_proposals(analogs, query_mass, tol_ppm=10.0, max_analogs=30, cap=900, *, legacy_analogs=None):
    if type(max_analogs) is not int or max_analogs < 15 or type(cap) is not int or cap < 300:
        raise ValueError("Expanded policy must preserve the 15-seed / 300-candidate legacy budget")
    legacy = propose_agg(analogs[:15] if legacy_analogs is None else legacy_analogs, query_mass, tol_ppm, 15, 300)
    expanded = propose_agg(analogs, query_mass, tol_ppm, max_analogs, cap, ordered_pairs=True)
    out = {smiles: {**record, "audit_legacy": 1} for smiles, record in legacy.items()}
    for smiles, record in expanded.items():
        if smiles not in out and len(out) < cap:
            out[smiles] = {**record, "audit_legacy": 0}
    return out
