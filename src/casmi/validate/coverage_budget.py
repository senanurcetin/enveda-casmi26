"""Descriptive candidate-coverage and ranking-loss accounting from aggregates.

No per-query outcomes are inferred. Conditional means describe the observed
covered subset and are unsuitable for comparing models with different pools.
"""
from __future__ import annotations

import math


def coverage_budget(summary, truth_candidates_present):
    n = summary["n"]
    if type(n) is not int or n < 1 or type(truth_candidates_present) is not int or not 0 <= truth_candidates_present <= n:
        raise ValueError("Positive integer population and bounded integer coverage count required")
    values = {}
    for field in ("mrr", "top1", "top10", "recall_at_25"):
        value = summary[field]
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Finite bounded aggregate metrics required")
        values[field] = value
    counts = {}
    for field in ("top1", "top10", "recall_at_25"):
        count = round(values[field] * n)
        if not math.isclose(count, values[field] * n, abs_tol=1e-9, rel_tol=0):
            raise ValueError("Aggregate hit rate does not correspond to integer molecules")
        counts[field] = count
    if not 0 <= counts["top1"] <= counts["top10"] <= counts["recall_at_25"] <= truth_candidates_present:
        raise ValueError("Ranked hits exceed candidate coverage or violate nested cutoffs")
    ceiling = truth_candidates_present / n
    mrr = values["mrr"]
    # Non-top-1 hits contribute at most 1/2; top-10 and top-25
    # hit groups contribute at least 1/10 and 1/25 respectively.
    lower = (counts["top1"] + (counts["top10"] - counts["top1"]) / 10
             + (counts["recall_at_25"] - counts["top10"]) / 25) / n
    upper = (counts["top1"] + (counts["top10"] - counts["top1"]) / 2
             + (counts["recall_at_25"] - counts["top10"]) / 11) / n
    if not lower - 1e-12 <= mrr <= upper + 1e-12:
        raise ValueError("MRR contradicts the observed top-1/top-10/top-25 counts")
    return {
        "n_molecules": n,
        "truth_candidates_present": truth_candidates_present,
        "truth_candidates_absent": n - truth_candidates_present,
        "top1_hits": counts["top1"], "top10_hits": counts["top10"], "top25_hits": counts["recall_at_25"],
        "covered_outside_top25": truth_candidates_present - counts["recall_at_25"],
        "coverage_ceiling": ceiling, "observed_mrr": mrr,
        "coverage_loss": 1 - ceiling, "within_coverage_ranking_loss": max(0.0, ceiling - mrr),
        "conditional_mrr_given_coverage": mrr / ceiling if truth_candidates_present else None,
        "conditional_hit25_given_coverage": counts["recall_at_25"] / truth_candidates_present if truth_candidates_present else None,
    }
