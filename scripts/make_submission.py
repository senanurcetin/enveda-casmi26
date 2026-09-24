"""Produce submission.csv for a test file with the full pipeline, timing every stage.

  python scripts/make_submission.py [--test PATH] [--ranker TAG] [--no-d3] [--sanity]

--sanity  the public test.parquet is 400 enveda-180 molecules copied from train, so their answers
          can be recovered; this reports the MRR on them. They are all Class 1 and drug-like, so
          the number checks plumbing and runtime only — it must never be used to pick a model.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, DATA_DIR, OUTPUT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE  # noqa: E402
from casmi.channels.analog import AnalogIndex  # noqa: E402
from casmi.channels.library import LibraryIndex  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.engine import BaselineEngine, EngineConfig  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.predict import Predictor, validate_submission  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402


class LoadedRanker:
    """The seed-bagged LightGBM models saved by scripts/train_ranker.py."""

    def __init__(self, tag):
        import lightgbm as lgb
        self.models = [lgb.Booster(model_file=str(p)) for p in sorted((ARTIFACT_DIR / "models").glob(f"{tag}_seed*.txt"))]
        assert self.models, f"no ranker models for {tag}"
        self.columns = self.models[0].feature_name()

    def predict(self, feats):
        return np.mean([m.predict(feats[self.columns]) for m in self.models], axis=0)


def sanity_answers(test):
    """Recover the public test answers by exact spectrum match against train (see --sanity)."""
    keys = {(round(r.precursor_mz, 4), len(r.ms2_mzs), round(float(r.ms2_mzs[0]), 4)): r.molecule_id for r in test.itertuples()}
    pf, found = pq.ParquetFile(DATA_DIR / "train.parquet"), {}
    precs = {k[0] for k in keys}
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=["precursor_mz", "normalized_smiles"]).to_pandas()
        cand = np.flatnonzero(t.precursor_mz.round(4).isin(precs).to_numpy())
        if len(cand):
            mz = pf.read_row_group(rg, columns=["ms2_mzs"]).column(0)
            for i in cand:
                m = mz[int(i)].as_py()
                mid = keys.get((round(t.precursor_mz.iat[i], 4), len(m), round(float(m[0]), 4)))
                if mid:
                    found[mid] = t.normalized_smiles.iat[i]
    return pd.DataFrame({"molecule_id": list(found), "smiles": list(found.values())})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", default=str(DATA_DIR / "test.parquet"))
    ap.add_argument("--ranker", default=None, help="ranker tag in artifacts/models (default: interim rule)")
    ap.add_argument("--no-d3", action="store_true")
    ap.add_argument("--sanity", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    stage = lambda name: print(f"[{time.time() - t0:6.0f}s] {name}", flush=True)
    cfg = EngineConfig()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE)
    library = LibraryIndex(store, key_of_row, cfg.library)
    analog = AnalogIndex(store, key_of_row, cfg.analog, cache=ARTIFACT_DIR / "analog_index.npz")
    engine = BaselineEngine(pool, library, cfg, analog, fps)
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    stage("assets loaded")

    ranker = LoadedRanker(args.ranker) if args.ranker else None
    fp_models, device = [], "cpu"
    if ranker is not None and not args.no_d3 and any(c.startswith("d3_") for c in ranker.columns):
        import torch
        from casmi.channels.fpmodel import load
        device = "cuda" if torch.cuda.is_available() else "cpu"
        fp_models = [load(p, device) for p in sorted((ARTIFACT_DIR / "models").glob("fp_fold*.pt"))]
    stage(f"ranker: {args.ranker or 'interim rule'} | D3 models: {len(fp_models)}")

    test = pd.read_parquet(args.test)
    predictor = Predictor(engine, pool, keys, ranker, fp_models, device)
    sub = predictor.predict_file(test)
    validate_submission(sub, test)
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    out = Path(OUTPUT_DIR) / "submission.csv"
    sub.to_csv(out, index=False)
    stage(f"wrote {out} ({len(sub)} molecules, {len(predictor.failures)} fell back)")
    print("stage timings (s):", {k: round(v, 1) for k, v in predictor.timings.items()})
    for mid, tb in predictor.failures[:3]:
        print(f"  fallback {mid}: {tb.splitlines()[-1]}")

    if args.sanity:
        from casmi.validate.metric import score_submission
        answers = sanity_answers(test)
        rr = score_submission(answers, sub)
        print(f"sanity MRR on the public (Class-1, drug-like) file: {rr.mean():.4f} over {len(rr)} molecules")


if __name__ == "__main__":
    main()
