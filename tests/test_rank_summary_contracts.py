"""Reject malformed rank arrays before they become scientific score evidence."""
import numpy as np
import pytest

from casmi.validate.metric import summarize
from casmi.validate.paired_bootstrap import paired_difference


@pytest.mark.parametrize("ranks", [[], [[1, 2]], [26], [-1], [1.5], [1.0], [np.nan], [np.inf], [True], ["1"]])
def test_invalid_ranks_are_rejected_by_both_summary_paths(ranks):
    with pytest.raises(ValueError):
        summarize(ranks, n_boot=100)
    with pytest.raises(ValueError):
        paired_difference(ranks, [1], n_boot=100)
    with pytest.raises(ValueError):
        paired_difference([1], ranks, n_boot=100)


@pytest.mark.parametrize("kwargs", [{"n_boot": 99}, {"n_boot": 100.5}, {"n_boot": True}, {"seed": -1}, {"seed": True}, {"seed": 1.5}])
def test_invalid_bootstrap_parameters_fail_explicitly(kwargs):
    with pytest.raises(ValueError):
        summarize([1, 2], **kwargs)
    with pytest.raises(ValueError):
        paired_difference([1, 2], [1, 2], **kwargs)


def test_valid_fixed_seed_summaries_match_direct_reciprocal_rank_bootstrap():
    ranks = np.array([1, 2, 25, 0], dtype=np.int64)
    rr = np.array([1, .5, .04, 0])
    ids = np.random.default_rng(7).integers(4, size=(200, 4))
    expected = np.quantile(rr[ids].mean(axis=1), [.025, .975])
    result = summarize(ranks, n_boot=200, seed=7)
    assert result.mrr == pytest.approx(rr.mean())
    assert (result.ci_low, result.ci_high) == pytest.approx(expected)
    assert (result.top1, result.top10, result.recall_at_25) == (.25, .5, .75)
    assert paired_difference(ranks, ranks, n_boot=200, seed=7) == (0, 0, 0)
    with pytest.raises(ValueError, match="equal length"):
        paired_difference([1], [1, 2], n_boot=200)
