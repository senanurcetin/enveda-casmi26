"""E — LightGBM LambdaRank over per-molecule candidate lists, seed-bagged.

Grouping is by (query, simulated class): each group is one molecule's candidate list with exactly
one relevant item (or none, for recall misses — those groups carry no gradient and are dropped).
"""
from dataclasses import dataclass, field

import lightgbm as lgb
import numpy as np
import pandas as pd

MONOTONE_UP = {"d1_max", "d1_mean", "an_p1", "an_p2", "an_p4", "an_self"}


@dataclass
class RankerConfig:
    seeds: tuple = (0, 1, 2, 3, 4)
    params: dict = field(default_factory=lambda: dict(
        objective="lambdarank", metric="map", eval_at=[25], lambdarank_truncation_level=25,
        learning_rate=0.05, num_leaves=31, min_data_in_leaf=50, feature_fraction=0.8,
        bagging_fraction=0.8, bagging_freq=1, lambda_l2=1.0, verbose=-1))
    num_boost_round: int = 2000
    early_stopping: int = 100


class Ranker:
    def __init__(self, cfg: RankerConfig = RankerConfig()):
        self.cfg, self.models, self.columns = cfg, [], None

    @staticmethod
    def _groups(df):
        return df.groupby(["qid", "cls"], sort=False).size().to_numpy()

    def fit(self, train: pd.DataFrame, valid: pd.DataFrame, columns, class_weight=None):
        self.columns = list(columns)
        keep = lambda d: d[d.groupby(["qid", "cls"]).label.transform("max") > 0]
        train, valid = keep(train), keep(valid)
        mono = [1 if c in MONOTONE_UP else 0 for c in self.columns]
        w = train.cls.map(class_weight).to_numpy() if class_weight else None
        dtr = lgb.Dataset(train[self.columns], train.label, group=self._groups(train), weight=w, free_raw_data=False)
        dva = lgb.Dataset(valid[self.columns], valid.label, group=self._groups(valid), reference=dtr)
        self.models, self.best_iters = [], []
        for seed in self.cfg.seeds:
            params = {**self.cfg.params, "seed": seed, "monotone_constraints": mono}
            m = lgb.train(params, dtr, self.cfg.num_boost_round, valid_sets=[dva],
                          callbacks=[lgb.early_stopping(self.cfg.early_stopping, verbose=False)])
            self.models.append(m)
            self.best_iters.append(m.best_iteration)
        return self

    def predict(self, feats: pd.DataFrame) -> np.ndarray:
        X = feats[self.columns]
        return np.mean([m.predict(X, num_iteration=m.best_iteration) for m in self.models], axis=0)

    def importance(self) -> pd.Series:
        imp = np.mean([m.feature_importance("gain") for m in self.models], axis=0)
        return pd.Series(imp, index=self.columns).sort_values(ascending=False)
