"""MRR@25 with the official metric's semantics, plus the per-molecule detail we need for analysis.

Semantics reproduced from the official scoring notebook (`metric/casmi-mean-reciprocal-rank`):
  * rank = position in the list; an unparseable or duplicate guess still consumes its position,
  * only the first correct guess counts (reciprocal rank, not average precision),
  * a molecule absent from the submission scores 0.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from casmi.chem.canon import canonical_ik14, heavy_atom_key
from casmi.validate.paired_bootstrap import paired_difference
from casmi.validate.ranked_slots import split_ranked_guesses

MAX_RANK = 25


def reciprocal_rank(guesses, answer_ik14s, answer_heavy_keys, max_rank=MAX_RANK):
    for rank, guess in enumerate(guesses[:max_rank], 1):
        if heavy_atom_key(guess) not in answer_heavy_keys:
            continue
        if canonical_ik14(guess) in answer_ik14s:
            return 1.0 / rank
    return 0.0


def score_submission(solution: pd.DataFrame, submission: pd.DataFrame, max_rank=MAX_RANK) -> pd.Series:
    """Per-molecule reciprocal rank. `solution` has molecule_id + smiles (comma-separated if several
    structures are accepted); `submission` has molecule_id + semicolon-joined guesses."""
    if submission["molecule_id"].duplicated().any():
        raise ValueError("Each molecule_id must appear exactly once")
    # Never compress malformed empty slots: a blank candidate still occupies
    # its original rank, even though a valid competition submission forbids it.
    guesses = {m: split_ranked_guesses(str(s))
               for m, s in zip(submission["molecule_id"], submission["smiles"])}
    too_many = [m for m, g in guesses.items() if len(g) > max_rank]
    if too_many:
        raise ValueError(f"More than {max_rank} guesses for {too_many[:3]}")
    out = {}
    for m, answer in zip(solution["molecule_id"], solution["smiles"]):
        accepted = [s for s in str(answer).split(",") if s]
        ik14s = {k for k in map(canonical_ik14, accepted) if k}
        heavy = {k for k in map(heavy_atom_key, accepted) if k}
        if not ik14s:
            raise ValueError(f"Unparseable solution SMILES for {m}")
        out[m] = reciprocal_rank(guesses.get(m, []), ik14s, heavy, max_rank)
    return pd.Series(out, name="rr")


def mrr(solution, submission, max_rank=MAX_RANK) -> float:
    return float(score_submission(solution, submission, max_rank).mean())


# ----------------------------------------------------------------------------- key-level scoring
def rank_of_truth(ranked_ik14s, true_ik14, max_rank=MAX_RANK):
    """1-based rank of the true key in an already canonicalised, deduplicated list (0 = absent).
    Validation runs on keys directly, which is exact once the output writer guarantees
    canonical, deduplicated lists — and ~1000x faster than re-parsing SMILES."""
    for rank, key in enumerate(ranked_ik14s[:max_rank], 1):
        if key == true_ik14:
            return rank
    return 0


@dataclass
class Summary:
    mrr: float
    ci_low: float
    ci_high: float
    top1: float
    top10: float
    recall_at_25: float
    n: int

    def __str__(self):
        return (f"MRR {self.mrr:.4f} [{self.ci_low:.4f}, {self.ci_high:.4f}]  top1 {self.top1:.3f}  "
                f"top10 {self.top10:.3f}  hit@25 {self.recall_at_25:.3f}  n={self.n}")


def summarize(ranks, n_boot=2000, seed=0) -> Summary:
    """MRR with a molecule-level bootstrap 95% CI, from 1-based ranks (0 = not in top 25)."""
    ranks = np.asarray(ranks)
    rr = np.where(ranks > 0, 1.0 / np.maximum(ranks, 1), 0.0)
    rng = np.random.default_rng(seed)
    boots = rr[rng.integers(0, len(rr), size=(n_boot, len(rr)))].mean(axis=1)
    return Summary(float(rr.mean()), float(np.quantile(boots, .025)), float(np.quantile(boots, .975)),
                   float((ranks == 1).mean()), float(((ranks > 0) & (ranks <= 10)).mean()),
                   float((ranks > 0).mean()), len(rr))


