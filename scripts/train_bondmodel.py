"""Fit D5 bond-event propensities on library spectra of training structures, then screen D5 alone on
V-in / V-real Class-2 candidate lists (same protocol as probe_d4.py; D4 int_max there: 0.404 / 0.558).

  python scripts/train_bondmodel.py [--n-spectra 30000]
Writes artifacts/models/bondmodel.pkl.
"""
import argparse
import importlib.util
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, SEED  # noqa: E402
from casmi.channels.bondmodel import BondModel, cleavage_events, observed  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402

spec = importlib.util.spec_from_file_location("brd", ROOT / "scripts" / "build_rank_data.py")
brd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(brd)


def mrr(v, lab):
    r = (v > v[lab].max()).sum() + 1
    return 1.0 / r if r <= 25 else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-spectra", type=int, default=30000)
    ap.add_argument("--n-probe", type=int, default=250)
    args = ap.parse_args()
    t0 = time.time()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    meta = store.meta
    key_of_row = meta.train_ik14.to_numpy(dtype=object)
    smiles = brd.structure_smiles()
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet"))
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    held = vin_skel | brd.vreal_skeletons()

    ok = np.flatnonzero(meta.ingest_lib.isin(brd.TRAIN_LIBS).to_numpy() & meta.test_adduct.to_numpy()
                        & (meta.n_peaks.to_numpy() >= 5))
    rng = np.random.default_rng(SEED)
    rows = rng.choice(ok, min(args.n_spectra * 2, len(ok)), replace=False)
    model, used, ev_cache = BondModel(), 0, {}
    for r in rows:
        k = key_of_row[r]
        if k not in smiles:
            continue
        if k not in ev_cache:
            ev_cache[k] = None if skeleton_key(smiles[k]) in held else cleavage_events(smiles[k])
        ev = ev_cache[k]
        if not ev:
            continue
        mz, it = store.peaks(r)
        m = meta.iloc[r]
        mode = "neg" if str(m.ionization_mode).lower().startswith("n") else "pos"
        rel = it / it.max() if len(it) and it.max() > 0 else it
        model.add(ev, observed(ev, np.asarray(mz, float), np.asarray(rel, float), mode, str(m.adduct), float(m.precursor_mz)), mode)
        used += 1
        if used >= args.n_spectra:
            break
    model.finalize()
    (ARTIFACT_DIR / "models").mkdir(exist_ok=True)
    (ARTIFACT_DIR / "models" / "bondmodel.pkl").write_bytes(pickle.dumps(model.counts))
    top = sorted(((c[0], k) for k, c in model.counts.items()), reverse=True)[:5]
    print(f"{used:,} spectra, {len(model.counts):,} event types, base rate {model.base:.3f} ({time.time() - t0:.0f}s)")
    print("most frequent:", [(k, round(model.p(*k), 3)) for _, k in top], flush=True)

    spectra = {q.qid: q.spectra for q in vin}
    spectra.update({q.qid: q.spectra for q, _ in brd.vreal_queries(store, key_of_row, vin_skel)[0]})
    for name in ("vin", "vreal"):
        t = pd.read_parquet(ARTIFACT_DIR / f"rank_{name}_gen2_0.parquet",
                            columns=["qid", "cls", "pool_key", "label", "gen"] + (["stratum"] if name == "vreal" else []))
        t = t[(t.cls == "C2") & (t.gen == 0)]
        if name == "vreal":
            t = t[t.stratum == "coco"]
        t = t[t.groupby("qid").label.transform("max") > 0]
        t = t[t.qid.isin(t.qid.unique()[: args.n_probe])]
        r_ll, r_w = [], []
        for qid, g in t.groupby("qid", sort=False):
            sc = np.array([model.score(cleavage_events(pool.smiles_of(k)), spectra[qid]) for k in g.pool_key])
            lab = g.label.to_numpy() == 1
            r_ll.append(mrr(sc[:, 0], lab))
            r_w.append(mrr(sc[:, 1], lab))
        print(f"{name}: D5 ll {np.mean(r_ll):.3f} | D5 weighted-observed {np.mean(r_w):.3f}  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
