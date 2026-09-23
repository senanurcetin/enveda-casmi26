"""End-to-end prediction for one molecule or a whole test file — the code path the submission runs.

    engine features (D1, D2, D3) -> priors + list context -> ranker score -> metric-key dedupe
    -> 25 SMILES, best first

Every molecule is wrapped in a fallback: if anything fails, it is answered with its mass-window
candidates ordered by mass error, so a single bad molecule can never cost a row (a missing row
scores 0 and an exception would lose the whole submission).
"""
import time
import traceback

import numpy as np
import pandas as pd

from casmi.channels.fpmodel import candidate_scores
from casmi.rank.features import add_context, add_priors


class Predictor:
    def __init__(self, engine, pool, keys, ranker=None, fp_models=None, device="cpu"):
        self.engine, self.pool, self.keys = engine, pool, keys
        self.ranker, self.fp_models, self.device = ranker, fp_models or [], device
        self.timings = {"features": 0.0, "d3": 0.0, "rank": 0.0, "dedupe": 0.0}
        self.failures = []

    # ------------------------------------------------------------------ one molecule
    def _d3(self, f, spectra):
        from casmi.channels.fpinfer import molecule_logp
        t = time.perf_counter()
        lp_on, lp_off = molecule_logp(self.fp_models, spectra, self.device)
        rows = np.fromiter((self.pool.row_of[k] for k in f.pool_key), np.int64, len(f))
        bits = np.unpackbits(self.engine.fps[rows], axis=1).astype(np.float32)
        f["d3_ll"], f["d3_cos"] = candidate_scores(lp_on, lp_off, bits)
        for c in ("d3_ll", "d3_cos"):
            v = f[c].to_numpy()
            f[f"{c}_rank"] = pd.Series(v).rank(ascending=False, method="min").to_numpy()
            f[f"{c}_gap"] = v.max() - v
            sd = v.std()
            f[f"{c}_z"] = (v - v.mean()) / sd if sd > 0 else 0.0
        self.timings["d3"] += time.perf_counter() - t
        return f

    def rank(self, spectra, top_n=25):
        f = self.engine.features([dict(s) for s in spectra])
        if f.empty:
            return []
        f = add_context(add_priors(f, self.pool))
        if self.fp_models:
            f = self._d3(f, spectra)
        t = time.perf_counter()
        if self.ranker is not None:
            order = np.argsort(-self.ranker.predict(f), kind="stable")
            ranked = f.pool_key.to_numpy()[order][:60].tolist()
        else:
            ranked = self.engine.rank(f)
        self.timings["rank"] += time.perf_counter() - t
        t = time.perf_counter()
        out = self.keys.dedupe_ranked(ranked, self.pool.smiles_of, top_n)
        self.timings["dedupe"] += time.perf_counter() - t
        return out                                   # [(pool_key, metric_key)] best first

    def fallback_keys(self, spectra, top_n=25):
        """Mass-window candidates ordered by mass error: the answer when the model path fails."""
        from casmi.engine import molecule_mass
        win = self.pool.window(molecule_mass(spectra), 20.0)
        keys = self.pool.key[win.idx][np.argsort(np.abs(win.ppm))].tolist()
        return self.keys.dedupe_ranked(keys, self.pool.smiles_of, top_n)

    # ------------------------------------------------------------------ whole test file
    def predict_file(self, test: pd.DataFrame, log_every=50) -> pd.DataFrame:
        rows, t0 = [], time.time()
        groups = test.groupby("molecule_id", sort=False)
        for i, (mid, g) in enumerate(groups):
            spectra = [dict(mz=np.asarray(r.ms2_mzs, np.float64), it=np.asarray(r.ms2_normalized_intensities, np.float64),
                            precursor_mz=float(r.precursor_mz), adduct=str(r.adduct),
                            ce_ev=float(np.mean(r.collision_energy_ev)) if r.collision_energy_ev is not None and len(r.collision_energy_ev) else np.nan,
                            mode=str(r.ionization_mode)) for r in g.itertuples()]
            try:
                picked = self.rank(spectra)
            except Exception:
                self.failures.append((mid, traceback.format_exc(limit=2)))
                picked = []
            if len(picked) < 25:                     # refill by metric key, so no structure appears twice
                seen = {mk for _, mk in picked}
                picked += [(pk, mk) for pk, mk in self.fallback_keys(spectra) if mk not in seen][: 25 - len(picked)]
            smiles = [self.pool.smiles_of(pk) for pk, _ in picked]
            rows.append((mid, ";".join(smiles) if smiles else "C"))   # never an empty field
            if log_every and (i + 1) % log_every == 0:
                print(f"  {i + 1}/{groups.ngroups} molecules, {time.time() - t0:.0f}s", flush=True)
        self.keys.save()
        return pd.DataFrame(rows, columns=["molecule_id", "smiles"])


def validate_submission(sub: pd.DataFrame, test: pd.DataFrame, max_rank=25):
    """The format rules the scorer enforces, checked before writing."""
    assert list(sub.columns) == ["molecule_id", "smiles"], sub.columns
    assert not sub.isna().any().any(), "nulls"
    assert sub.molecule_id.is_unique, "duplicate molecule_id"
    missing = set(test.molecule_id) - set(sub.molecule_id)
    assert not missing, f"{len(missing)} molecules missing"
    n = sub.smiles.str.split(";").map(len)
    assert (n >= 1).all() and (n <= max_rank).all(), "guess count out of range"
