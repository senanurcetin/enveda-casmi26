"""Prespecified molecule-level diagnostic slices for CASMI research validation.

Groupings use only input metadata. All outcomes are descriptive, with small
cells suppressed and no post-hoc statistical ACCEPT/REJECT claims. This module
requires only the Python standard library and signed validation reports.
"""
from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping
from statistics import median

from casmi.validate.run_alignment import CLASSES, align_runs

FIELDS = ("adduct", "ionization_mode", "spectra_bin", "precursor_mz_bin")
MASS_BINS = ("<300", "300-600", ">=600", "unknown")


def describe_query(spectra) -> dict:
    """Compute one query's *original* non-label metadata, before perturbation."""
    if not isinstance(spectra, (list, tuple)) or not spectra:
        raise ValueError("Query must contain at least one selected spectrum")
    def shared(field):
        values = {str(s[field]) for s in spectra}
        return next(iter(values)) if len(values) == 1 else "mixed"

    n = len(spectra)
    masses = [float(s["precursor_mz"]) for s in spectra]
    masses = [m for m in masses if math.isfinite(m)]
    midpoint = median(masses) if masses else math.nan
    mass_bin = ("unknown" if not math.isfinite(midpoint) else
                "<300" if midpoint < 300 else
                "300-600" if midpoint < 600 else ">=600")
    return {
        "adduct": shared("adduct"),
        "ionization_mode": shared("mode"),
        "spectra_bin": "1" if n == 1 else "2-3" if n <= 3 else "4+",
        "precursor_mz_bin": mass_bin,
    }


def _checked_contexts(run, ids, label):
    context = run.get("query_context")
    if not isinstance(context, list) or len(context) != len(ids):
        raise ValueError(f"{label}: missing or incomplete query_context")
    by_id = {}
    for qid, item in zip(ids, context):
        if not isinstance(item, Mapping) or set(item) != set(FIELDS):
            raise ValueError(f"{label}: invalid query_context fields")
        row = dict(item)
        if any(not isinstance(row[field], str) or not row[field]
               for field in FIELDS):
            raise ValueError(f"{label}: empty or malformed query_context value")
        if row["spectra_bin"] not in ("1", "2-3", "4+"):
            raise ValueError(f"{label}: invalid spectra-count bin")
        if row["precursor_mz_bin"] not in MASS_BINS:
            raise ValueError(f"{label}: invalid precursor-mass bin")
        by_id[qid] = row
    return by_id


def _mrr(ranks, indices):
    return sum(1.0 / ranks[i] for i in indices if ranks[i]) / len(indices)


def diagnostic_slices(reference: Mapping, candidate: Mapping, *, min_n: int = 15) -> dict:
    """Compare known, query-aligned context slices; NEVER infer superiority."""
    if type(min_n) is not int or min_n < 2:
        raise ValueError("min_n must be an integer >= 2")
    aligned = align_runs(reference, candidate)
    ids = reference["query_ids"]
    b_ids = candidate["query_ids"]
    context_a = _checked_contexts(reference, ids, "Reference")
    context_b = _checked_contexts(candidate, b_ids, "Candidate")
    if context_a != context_b:
        raise ValueError("Query input metadata differ across validation runs")
    groups = {}
    for field in FIELDS:
        cells = defaultdict(list)
        for index, qid in enumerate(ids):
            cells[context_a[qid][field]].append(index)
        records = []
        for value in sorted(cells):
            subset = cells[value]
            row = {"value": value, "n_molecules": len(subset)}
            if len(subset) < min_n:
                row.update(status="SUPPRESSED_SMALL_CELL")
            else:
                classes = {}
                for cls in CLASSES:
                    a = _mrr(aligned["a"][cls], subset)
                    b = _mrr(aligned["b"][cls], subset)
                    classes[cls] = {
                        "reference_mrr": a,
                        "candidate_mrr": b,
                        "candidate_minus_reference": b - a,
                    }
                row.update(status="DESCRIPTIVE_ONLY", classes=classes)
            records.append(row)
        groups[field] = records
    return {
        "status": "DESCRIPTIVE_ONLY",
        "source_evidence": (
            "VERIFIED_REQUIRED_INPUTS" if aligned["decision_eligible"]
            else "QUERY_PAIRED_ONLY_SOURCE_UNVERIFIED"
        ),
        "n_molecules": aligned["n"],
        "min_cell_size": min_n,
        "slices": groups,
        "interpretation": (
            "Exploratory error characterization, not a model-improvement "
            "verdict. Multiple groups and novelty classes are inspected "
            "without multiplicity correction or independent holdout."
        ),
    }
