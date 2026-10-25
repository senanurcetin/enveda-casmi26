"""Candidate quarantine must not delete labels or independent reference exclusions."""
import numpy as np
import pandas as pd

from casmi.chem.canon import canonical_ik14
from casmi.chem.keys import MetricKeys
from casmi.pool.pool import CandidatePool
from casmi.validate.knockout import equivalent_keys, knockout_for, Query


def test_absent_candidate_truth_still_has_metric_label_and_reference_exclusion(tmp_path):
    truth = canonical_ik14("CCO")
    pool = CandidatePool(pd.DataFrame({"pool_key": ["decoy"], "smiles": ["COC"], "mass": [46.041865]}))
    keys = MetricKeys(tmp_path / "metric.parquet")
    label, retained = equivalent_keys(truth, 46.041865, pool, keys, truth_smiles="CCO")
    assert label == truth
    assert retained == set()
    # The runner unions independent training equivalents with retained pool
    # equivalents; the absent truth must still hide its cross-library references.
    same = {truth} | retained
    query = Query("synthetic", truth, [], np.array([0]), np.array([0, 1]))
    ko = knockout_for(query, "C2", {truth: np.array([0, 1, 2])}, same)
    np.testing.assert_array_equal(ko.excluded_rows, [0, 1, 2])
    assert truth not in pool.row_of
