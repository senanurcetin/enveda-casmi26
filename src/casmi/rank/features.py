"""Ranker features: channel outputs + knockout-consistent priors + within-molecule list context.

Two leakage rules are enforced here, because the ranker would otherwise learn them:
  * priors derived from training-set membership (has spectra, spectrum count, which libraries)
    are recomputed after the knockout — a simulated Class-2 answer must look exactly like a real
    one, i.e. like a structure with no reference spectra;
  * list-context features are computed per molecule over the candidates the engine actually saw.
"""
import numpy as np
import pandas as pd

from casmi.pool.pool import SOURCE_BITS

CHANNEL_COLS = ["d1_max", "d1_mean", "an_p1", "an_p2", "an_p4", "an_sim", "an_tani", "an_self"]
CONTEXT_OF = ["d1_max", "an_p4", "an_p2", "an_tani"]


def add_priors(feats: pd.DataFrame, pool, knocked_keys=frozenset()) -> pd.DataFrame:
    t = pool.table
    rows = np.fromiter((pool.row_of.get(k, -1) for k in feats.pool_key), np.int64, len(feats))
    gen = rows < 0                                     # generated structures: in no source
    knocked = feats.pool_key.isin(knocked_keys).to_numpy() | gen
    rows = np.where(gen, 0, rows)
    src = np.where(gen, 0, t.sources.to_numpy()[rows])
    in_train = ((src & SOURCE_BITS["train"]) > 0) & ~knocked
    feats = feats.assign(
        prior_in_train=in_train.astype(np.int8),
        prior_in_coconut=((src & SOURCE_BITS["coconut"]) > 0).astype(np.int8),
        prior_enveda180_only=(((src & SOURCE_BITS["enveda180_only"]) > 0) & ~knocked).astype(np.int8),
        prior_log_spectra=np.where(knocked, 0.0, np.log1p(t.n_spectra.to_numpy()[rows])),
        prior_np_likeness=np.where(gen, np.nan, t.np_likeness.to_numpy(dtype=float)[rows]),
        prior_has_sugar=np.where(gen, np.nan, t.has_sugar.map({True: 1.0, False: 0.0}).to_numpy(dtype=float)[rows]),
        abs_ppm=np.abs(feats.ppm.to_numpy()),
    )
    return feats


def add_context(feats: pd.DataFrame) -> pd.DataFrame:
    """Rank, gap-to-best and z-score of key channels within the molecule's candidate list."""
    out = {}
    for c in CONTEXT_OF:
        if c not in feats:
            continue
        v = feats[c].to_numpy(dtype=float)
        best = v.max() if len(v) else 0.0
        sd = v.std() if len(v) > 1 else 0.0
        out[f"{c}_rank"] = pd.Series(v).rank(ascending=False, method="min").to_numpy()
        out[f"{c}_gap"] = best - v
        out[f"{c}_z"] = (v - v.mean()) / sd if sd > 0 else np.zeros_like(v)
    out["log_n_cand"] = np.full(len(feats), np.log1p(len(feats)))
    return feats.assign(**out)


def feature_columns(feats: pd.DataFrame):
    drop = {"pool_key", "ppm", "n_cand", "label", "qid", "cls", "gen_smiles", "stratum"}
    return [c for c in feats.columns if c not in drop]
