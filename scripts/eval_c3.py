"""Class-3 simulation for the biotransformation hedge.

For every V-in and V-real query: remove the structure (and its tautomer twins) from the library,
search spectral analogs exactly as the submission does, propose transformed analogs at the query
mass, and record whether the truth is among the proposals and where each scoring rule puts it.
Scores for a proposal (all available at test time):
  sim     spectral similarity of the analog it came from
  d3      spectrum → fingerprint log-likelihood (out-of-fold model)
  d4      fragment-explained intensity share
Output: runs/c3_proposals.parquet (one row per proposal, with truth flag) + a summary.
"""
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, PROJECT_ROOT  # noqa: E402
from casmi.channels.analog import AnalogIndex  # noqa: E402
from casmi.channels.fpinfer import molecule_logp  # noqa: E402
from casmi.channels.fpmodel import candidate_scores, load  # noqa: E402
from casmi.channels.fragments import FragmentScorer  # noqa: E402
from casmi.channels.library import LibraryIndex  # noqa: E402
from casmi.chem.canon import canonical_ik14, skeleton_key  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.engine import EngineConfig, molecule_mass  # noqa: E402
from casmi.generate.transforms import propose  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.clean import clean_peaks_fast  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries, equivalent_keys, knockout_for  # noqa: E402

spec = importlib.util.spec_from_file_location("brd", ROOT / "scripts" / "build_rank_data.py")
brd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(brd)


def morgan_bits(smiles_list):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out = np.zeros((len(smiles_list), 2048), np.float32)
    for i, s in enumerate(smiles_list):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out[i] = gen.GetFingerprintAsNumPy(m)
    return out


def main():
    t0 = time.time()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    cfg = EngineConfig()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet"))
    library = LibraryIndex(store, key_of_row, cfg.library)
    analog = AnalogIndex(store, key_of_row, cfg.analog, cache=ARTIFACT_DIR / "analog_index.npz")
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    smiles = brd.structure_smiles()
    models = {f: load(ARTIFACT_DIR / "models" / f"fp_fold{f}.pt", device) for f in (0, 1)}
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    vreal, stratum_of = brd.vreal_queries(store, key_of_row, vin_skel)
    items = [(q, -1, "vin") for q in vin] + [(q, f, stratum_of[q.pool_key]) for q, f in vreal]
    frag = FragmentScorer(lambda k: gen_smiles[k])
    gen_smiles = {}
    print(f"{len(items)} queries loaded ({time.time() - t0:.0f}s)", flush=True)

    rows, summary = [], []
    for qi, (q, fold, stratum) in enumerate(items):
        mass = molecule_mass(q.spectra, cfg.mass_offset_ppm)
        truth, same = equivalent_keys(q.pool_key, mass, pool, keys, cfg.tol_ppm * 3, truth_smiles=smiles.get(q.pool_key))
        ko = knockout_for(q, "C3", library.rows_of, same | {q.pool_key})
        spectra = [dict(s) for s in q.spectra]
        for s in spectra:
            s["mz_c"], s["it_c"] = clean_peaks_fast(s["mz"], s["it"], s["precursor_mz"], library.cfg.clean)
        found = analog.search(spectra, mass, ko.excluded_rows)
        found = [(k, v) for k, v in found.items() if k not in same and k != q.pool_key]
        analogs = [(smiles[k], sim, am) for k, (sim, am) in found if k in smiles]
        props = propose(analogs, mass)
        truth_skel = skeleton_key(smiles[q.pool_key])
        hit = np.zeros(len(props), bool)
        for i, (s, _, _) in enumerate(props):
            if skeleton_key(s) == truth_skel and canonical_ik14(s) == truth:
                hit[i] = True
        summary.append((q.qid, stratum, len(analogs), len(props), bool(hit.any())))
        if props:
            use = [models[0], models[1]] if fold == -1 else [models[1 - fold]]
            lp_on, lp_off = molecule_logp(use, q.spectra, device)
            ll, cos = candidate_scores(lp_on, lp_off, morgan_bits([p[0] for p in props]))
            gkeys = []
            for i, (s, _, _) in enumerate(props):
                gk = f"gen{qi}_{i}"
                gen_smiles[gk] = s
                gkeys.append(gk)
            d4 = frag.score(q.spectra, gkeys)
            for i, (s, sim, path) in enumerate(props):
                rows.append((q.qid, stratum, s, sim, path, ll[i], cos[i], d4[i, 0], d4[i, 1], hit[i]))
            for gk in gkeys:
                frag.cache.pop(gk, None)
                gen_smiles.pop(gk, None)
        if (qi + 1) % 100 == 0:
            s_ = pd.DataFrame(summary, columns=["qid", "stratum", "n_analogs", "n_props", "hit"])
            print(f"  {qi + 1}/{len(items)}  truth proposed {s_.hit.mean():.3f}  median proposals "
                  f"{s_.n_props.median():.0f} ({time.time() - t0:.0f}s)", flush=True)
    keys.save()
    props = pd.DataFrame(rows, columns=["qid", "stratum", "smiles", "sim", "path", "d3_ll", "d3_cos", "d4_int", "d4_int_max", "hit"])
    props.to_parquet(PROJECT_ROOT / "runs" / "c3_proposals.parquet", index=False)
    s_ = pd.DataFrame(summary, columns=["qid", "stratum", "n_analogs", "n_props", "hit"])
    s_.to_parquet(PROJECT_ROOT / "runs" / "c3_summary.parquet", index=False)
    print(s_.groupby("stratum")[["n_analogs", "n_props", "hit"]].mean().round(3))
    print(f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
