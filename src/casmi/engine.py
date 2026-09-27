"""Per-molecule engine: neutral mass -> candidates -> channel features -> ranking.

The baseline ranks by D1 (direct library similarity), breaking ties by mass error. Later stages
replace `rank` with the learned ranker; the candidate/feature interface stays the same, so every
channel is validated through exactly the code path the submission uses.
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from casmi.channels.analog import AnalogConfig, propagate
from casmi.channels.library import LibraryConfig, LibraryIndex
from casmi.chem.masses import neutral_mass
from casmi.spectra.clean import clean_peaks_fast


@dataclass(frozen=True)
class EngineConfig:
    tol_ppm: float = 10.0
    mass_offset_ppm: float = 0.0      # calibrated instrument offset (measured +1.5 ppm on np-examples)
    top_n: int = 25
    library: LibraryConfig = field(default_factory=LibraryConfig)
    analog: AnalogConfig = field(default_factory=AnalogConfig)
    generate: bool = False            # Class-3 hedge: add biotransformed analogs as candidates
    gen_max_analogs: int = 15
    gen_cap: int = 300
    gen_consensus: bool = False       # v17 features: seed / route consensus and transform type
    np_score: bool = False            # RDKit NP-likeness for every candidate (chemistry, not membership)


def molecule_mass(spectra, offset_ppm=0.0):
    masses = [float(neutral_mass(s["precursor_mz"], s["adduct"])) for s in spectra]
    masses = [m for m in masses if np.isfinite(m)]
    m = float(np.median(masses)) if masses else np.nan
    return m / (1 + offset_ppm * 1e-6)


class BaselineEngine:
    def __init__(self, pool, library: LibraryIndex, cfg: EngineConfig = EngineConfig(), analog=None, fps=None):
        """`analog` (AnalogIndex) and `fps` (packed pool fingerprints) enable the D2 channel."""
        self.pool, self.library, self.cfg = pool, library, cfg
        self.analog, self.fps = analog, fps

    def features(self, spectra, excluded_rows=None, excluded_keys=None) -> pd.DataFrame:
        for s in spectra:   # clean once per query, with the library's cleaning settings
            s["mz_c"], s["it_c"] = clean_peaks_fast(s["mz"], s["it"], s["precursor_mz"], self.library.cfg.clean)
        mass = molecule_mass(spectra, self.cfg.mass_offset_ppm)
        win = self.pool.window(mass, self.cfg.tol_ppm, excluded_keys)
        keys = self.pool.key[win.idx]
        d1 = np.array([self.library.score(spectra, k, excluded_rows) for k in keys]).reshape(-1, 3)
        f = pd.DataFrame({"pool_key": keys, "ppm": win.ppm, "d1_max": d1[:, 0], "d1_mean": d1[:, 1],
                          "d1_nref": d1[:, 2].astype(int), "n_cand": len(keys)})
        analogs = None
        if self.analog is not None and len(keys):
            analogs = self.analog.search(spectra, mass, excluded_rows)
            for name, v in propagate(analogs, win.idx, self.pool.mass[win.idx], self.fps,
                                     self.pool.row_of, self.cfg.analog.powers).items():
                f[name] = v
        if self.cfg.generate:
            if analogs is None and self.analog is not None:
                analogs = self.analog.search(spectra, mass, excluded_rows)
            f = f.assign(gen=0, gen_sim=0.0, gen_depth=0)
            if self.cfg.gen_consensus:
                f = f.assign(gen_sim_sum=0.0, gen_n_seeds=0, gen_n_paths=0, gen_t1=-1, gen_t2=-1)
            g = self._generated(analogs or {}, mass, excluded_keys or set())
            if len(g):
                f = pd.concat([f, g], ignore_index=True)
                f["n_cand"] = len(f)
        if self.cfg.np_score and len(f):
            from casmi.chem.npscore import np_likeness
            f["np_score"] = [np_likeness(k, self.pool.smiles_of(k)) for k in f.pool_key]
        return f

    def _generated(self, analogs, mass, excluded_keys):
        """Biotransformed analogs at the query mass that are not already pool candidates."""
        from rdkit import Chem
        from rdkit.Chem.Descriptors import ExactMolWt
        from casmi.generate.transforms import propose, propose_agg
        order = sorted(analogs.items(), key=lambda kv: -kv[1][0])[: self.cfg.gen_max_analogs]
        seeds = [(self.pool.smiles_of(k), sim, am) for k, (sim, am) in order if k in self.pool.row_of]
        rows, seen = [], {}
        if self.cfg.gen_consensus:
            props = propose_agg(seeds, mass, self.cfg.tol_ppm, self.cfg.gen_max_analogs, self.cfg.gen_cap)
            items = sorted(props.items(), key=lambda kv: -kv[1]["sim_max"])
        else:
            items = [(smi, dict(sim_max=sim, depth=path.count("+") + 1))
                     for smi, sim, path in propose(seeds, mass, self.cfg.tol_ppm, self.cfg.gen_max_analogs, self.cfg.gen_cap)]
        for smi, r in items:
            mol = Chem.MolFromSmiles(smi)
            ik = Chem.MolToInchiKey(mol)[:14] if mol is not None else ""
            if not ik or (ik in self.pool.row_of and ik not in excluded_keys):
                continue
            if ik in seen:                                     # same structure reached as another SMILES
                if self.cfg.gen_consensus:
                    row = rows[seen[ik]]
                    row[5] += r["sim_sum"]
                    row[6] += r["n_seeds"]
                    row[7] += r["n_paths"]
                continue
            seen[ik] = len(rows)
            row = ["G" + ik, smi, (mass - ExactMolWt(mol)) / ExactMolWt(mol) * 1e6, r["sim_max"], r["depth"]]
            if self.cfg.gen_consensus:
                row += [r["sim_sum"], r["n_seeds"], r["n_paths"], r["t1"], r["t2"]]
            rows.append(row)
        if not rows:
            return pd.DataFrame()
        cols = ["pool_key", "gen_smiles", "ppm", "gen_sim", "gen_depth"]
        if self.cfg.gen_consensus:
            cols += ["gen_sim_sum", "gen_n_seeds", "gen_n_paths", "gen_t1", "gen_t2"]
        g = pd.DataFrame(rows, columns=cols)
        self.pool.register_generated(g.pool_key, g.gen_smiles)
        g = g.drop(columns="gen_smiles").assign(d1_max=0.0, d1_mean=0.0, d1_nref=0, gen=1)
        if analogs:
            keys = [k for k in analogs if k in self.pool.row_of]
            local = np.vstack([self.pool.packed_fps(g.pool_key.to_numpy(), self.fps),
                               self.fps[[self.pool.row_of[k] for k in keys]]])
            row_of = {k: len(g) + i for i, k in enumerate(keys)}
            gm = mass / (1 + g.ppm.to_numpy() * 1e-6)
            for name, v in propagate({k: analogs[k] for k in keys}, np.arange(len(g)), gm, local, row_of,
                                     self.cfg.analog.powers).items():
                g[name] = v
        return g

    def rank(self, feats: pd.DataFrame, depth=60):
        """Ranked pool keys (deeper than 25, so output deduplication can refill the list)."""
        if feats.empty:
            return []
        primary = feats.d1_max.to_numpy()
        secondary = feats.d1_mean.to_numpy()
        if "an_p4" in feats:   # interim rule until the learned ranker: best of direct and propagated evidence
            primary = np.maximum(primary, feats.an_p4.to_numpy())
            secondary = feats.an_p2.to_numpy()
        order = np.lexsort((np.abs(feats.ppm.to_numpy()), -secondary, -primary))
        return feats.pool_key.to_numpy()[order][:depth].tolist()
