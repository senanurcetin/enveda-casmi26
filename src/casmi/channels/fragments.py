"""D4 — fragment explanation: how much of the query spectrum a candidate structure can explain.

Candidates in one mass window are mostly isomers, so formula-level evidence cannot separate them;
bond topology can. Every candidate is broken at one or two bonds (combinatorial fragmentation in
the spirit of MetFrag, depth 2), every resulting connected piece becomes a neutral fragment mass,
and each fragment is ionised with the precursor's charge carrier and ±2 hydrogen rearrangements.
A query peak is explained when some fragment ion lies within tolerance.

Features per candidate (averaged and maxed over the molecule's spectra):
  d4_int   share of fragment-peak intensity explained
  d4_pk    share of fragment peaks explained
Peaks within 1.5 Da of the precursor are ignored (every candidate explains the precursor).
"""
import numpy as np
from numba import njit

from casmi.chem.masses import ELECTRON, PROTON
from casmi.spectra.clean import CleanConfig, clean_peaks_fast

H = 1.00782503207
NA = 22.9897692820
FRAG_CLEAN = CleanConfig(max_above_precursor=-1.5, min_rel_intensity=0.01, max_peaks=60)
TOL_DA, TOL_PPM = 0.01, 10.0
MAX_ATOMS = 160


@njit(cache=True)
def _find(parent, a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]
        a = parent[a]
    return a


@njit(cache=True)
def _pieces(n, atom_mass, bu, bv, skip1, skip2, out, n_out, total):
    parent = np.arange(n)
    for b in range(len(bu)):
        if b == skip1 or b == skip2:
            continue
        ra, rb = _find(parent, bu[b]), _find(parent, bv[b])
        if ra != rb:
            parent[ra] = rb
    acc = np.zeros(n)
    for a in range(n):
        acc[_find(parent, a)] += atom_mass[a]
    for a in range(n):
        m = acc[a]
        if m > 0.0 and m < total - 1e-6 and n_out < len(out):
            out[n_out] = m
            n_out += 1
    return n_out


@njit(cache=True)
def _fragments(atom_mass, bu, bv):
    n, m = len(atom_mass), len(bu)
    total = atom_mass.sum()
    out = np.empty(3 * (m + m * (m - 1) // 2) + 1)
    k = 0
    for i in range(m):
        k = _pieces(n, atom_mass, bu, bv, i, -1, out, k, total)
        for j in range(i + 1, m):
            k = _pieces(n, atom_mass, bu, bv, i, j, out, k, total)
    return out[:k]


def fragment_masses(mol):
    """Sorted unique neutral fragment masses of an RDKit molecule (heavy atoms carry their Hs)."""
    n = mol.GetNumAtoms()
    if n < 2 or n > MAX_ATOMS:
        return np.empty(0)
    atom_mass = np.array([_mono(a) for a in mol.GetAtoms()])
    bu = np.array([b.GetBeginAtomIdx() for b in mol.GetBonds()], np.int64)
    bv = np.array([b.GetEndAtomIdx() for b in mol.GetBonds()], np.int64)
    if len(bu) == 0:
        return np.empty(0)
    f = _fragments(atom_mass, bu, bv)
    return np.unique(np.round(f, 4))


_MONO = {}


def _mono(atom):
    from rdkit.Chem import GetPeriodicTable
    s = atom.GetSymbol()
    if s not in _MONO:
        pt = GetPeriodicTable()
        _MONO[s] = pt.GetMassForIsotope(s, pt.GetMostCommonIsotope(s))
    return _MONO[s] + atom.GetTotalNumHs() * H


def ion_masses(frags, mode, adduct):
    """Fragment ions for the precursor's charge carrier, with -2..+2 hydrogen rearrangements."""
    if mode.lower().startswith("n") or adduct.endswith("-"):
        carriers = [-PROTON]
    else:
        carriers = [PROTON] + ([NA - ELECTRON] if "Na" in adduct else [])
    shifts = np.array([c + k * H for c in carriers for k in (-2, -1, 0, 1, 2)])
    return np.sort((frags[:, None] + shifts[None, :]).ravel())


def explained(peaks_mz, peaks_it, ions):
    """(intensity share, peak share) of the peaks within tolerance of some ion."""
    if len(peaks_mz) == 0 or len(ions) == 0:
        return 0.0, 0.0
    pos = np.searchsorted(ions, peaks_mz)
    lo = ions[np.clip(pos - 1, 0, len(ions) - 1)]
    hi = ions[np.clip(pos, 0, len(ions) - 1)]
    d = np.minimum(np.abs(peaks_mz - lo), np.abs(peaks_mz - hi))
    hit = d <= np.maximum(TOL_DA, peaks_mz * TOL_PPM * 1e-6)
    return float(peaks_it[hit].sum() / peaks_it.sum()), float(hit.mean())


class FragmentScorer:
    """Caches fragment masses per structure; scores (query spectra × candidates)."""

    def __init__(self, smiles_of):
        self.smiles_of = smiles_of
        self.cache = {}

    def frags(self, key):
        f = self.cache.get(key)
        if f is None:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(self.smiles_of(key))
            f = fragment_masses(mol) if mol is not None else np.empty(0)
            self.cache[key] = f
        return f

    def prepare(self, spectra):
        """Cleaned fragment peaks per spectrum (sqrt intensity, precursor region removed)."""
        out = []
        for s in spectra:
            mz, it = clean_peaks_fast(s["mz"], s["it"], s["precursor_mz"], FRAG_CLEAN)
            mz, it = np.asarray(mz, np.float64), np.sqrt(np.asarray(it, np.float64))
            if len(mz):
                out.append((mz, it, str(s.get("mode", "")), str(s.get("adduct", ""))))
        return out

    def score(self, spectra, keys):
        peaks = self.prepare(spectra)
        d4 = np.zeros((len(keys), 4))
        if not peaks:
            return d4
        for i, k in enumerate(keys):
            f = self.frags(k)
            v = np.array([explained(mz, it, ion_masses(f, mode, adduct)) for mz, it, mode, adduct in peaks])
            d4[i] = (v[:, 0].mean(), v[:, 0].max(), v[:, 1].mean(), v[:, 1].max())
        return d4


D4_COLUMNS = ["d4_int", "d4_int_max", "d4_pk", "d4_pk_max"]


def add_d4(f, scorer, spectra):
    """D4 columns plus within-list rank / gap / z for the molecule's candidate table `f`."""
    d4 = scorer.score(spectra, f.pool_key.to_numpy())
    f = f.assign(**{c: d4[:, j] for j, c in enumerate(D4_COLUMNS)})
    import pandas as pd
    for c in ("d4_int", "d4_pk"):
        v = f[c].to_numpy()
        sd = v.std()
        f[f"{c}_rank"] = pd.Series(v).rank(ascending=False, method="min").to_numpy()
        f[f"{c}_gap"] = v.max() - v
        f[f"{c}_z"] = (v - v.mean()) / sd if sd > 0 else 0.0
    return f


def add_d4_table(tbl, spectra_of_qid, smiles_of, cache_path=None):
    """D4 for every (qid, cls) group of a ranker table; fragment masses cached on disk."""
    import pickle
    import pandas as pd
    scorer = FragmentScorer(smiles_of)
    if cache_path is not None and cache_path.exists():
        scorer.cache.update(pickle.loads(cache_path.read_bytes()))
    parts = [add_d4(g.reset_index(drop=True), scorer, spectra_of_qid[qid])
             for (qid, _), g in tbl.groupby(["qid", "cls"], sort=False)]
    if cache_path is not None:
        cache_path.write_bytes(pickle.dumps(scorer.cache, protocol=pickle.HIGHEST_PROTOCOL))
    return pd.concat(parts, ignore_index=True)
