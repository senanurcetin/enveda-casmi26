"""Ensure validation uses the actual prediction path without losing knockout scope."""
import numpy as np
import pandas as pd

from casmi.predict import Predictor


class Pool:
    row_of = {"candidate": 0}
    table = pd.DataFrame({"sources": [1], "n_spectra": [100],
                          "np_likeness": [0.5], "has_sugar": [False]})
    def smiles_of(self, key):
        return "CCO"


class Keys:
    def dedupe_ranked(self, ranked, smiles_of, top_n):
        return [(key, "metric-key") for key in ranked[:top_n]]


class Engine:
    def features(self, spectra, excluded_rows=None, excluded_keys=None):
        self.received = (excluded_rows, excluded_keys)
        return pd.DataFrame({"pool_key": ["candidate"], "ppm": [0.0],
                             "d1_max": [0.5], "d1_mean": [0.5], "d1_nref": [1]})


class Ranker:
    columns = ["prior_log_spectra"]
    def predict(self, features):
        self.received = features
        return np.ones(len(features))


def test_prediction_path_forwards_knockout_and_clears_withheld_reference_priors():
    engine, ranker = Engine(), Ranker()
    predictor = Predictor(engine, Pool(), Keys(), ranker=ranker)
    excluded_rows = np.array([2, 5])
    excluded_keys = {"equivalent"}
    result = predictor.rank([], excluded_rows=excluded_rows,
                            excluded_keys=excluded_keys, knocked_keys={"candidate"})
    assert engine.received[0] is excluded_rows
    assert engine.received[1] is excluded_keys
    assert ranker.received.prior_log_spectra.iat[0] == 0
    assert ranker.received.prior_in_train.iat[0] == 0
    assert result == [("candidate", "metric-key")]


def test_default_prediction_retains_existing_priors_and_has_no_exclusions():
    engine, ranker = Engine(), Ranker()
    Predictor(engine, Pool(), Keys(), ranker=ranker).rank([])
    assert engine.received == (None, None)
    assert ranker.received.prior_log_spectra.iat[0] > 0
    assert ranker.received.prior_in_train.iat[0] == 1
