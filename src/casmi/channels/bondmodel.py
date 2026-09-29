"""D5 — learned bond-cleavage likelihood.

D4 counts how much of the spectrum a candidate *can* explain, which favours candidates that
fragment into many pieces. D5 asks the reverse question as well: of the cleavages this candidate
*should* show, which do we see? Each cleavage event (one acyclic bond, or two bonds of the same
ring) has a type; its probability of producing an observed fragment ion is learned from library
spectra of training structures. A candidate's score is the Bernoulli log-likelihood of the
observed/unobserved pattern of its events, averaged over events and spectra.
"""
import numpy as np
from rdkit import Chem

from casmi.chem.masses import ELECTRON, PROTON

H = 1.00782503207
NA = 22.9897692820
TOL_DA, TOL_PPM = 0.01, 10.0
MIN_REL = 0.01


def _atom_masses(mol):
    pt = Chem.GetPeriodicTable()
    return np.array([pt.GetMassForIsotope(a.GetSymbol(), pt.GetMostCommonIsotope(a.GetSymbol()))
                     + a.GetTotalNumHs() * H for a in mol.GetAtoms()])


def _bond_key(b):
    a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
    ends = sorted([(a1.GetSymbol(), a1.GetDegree(), a1.GetIsAromatic()), (a2.GetSymbol(), a2.GetDegree(), a2.GetIsAromatic())])
    return f"{ends[0][0]}{ends[0][1]}{'a' if ends[0][2] else ''}-{b.GetBondTypeAsDouble():g}-{ends[1][0]}{ends[1][1]}{'a' if ends[1][2] else ''}"


def _split(mol, bonds, masses):
    """Masses of the two pieces after removing `bonds`, or None if the molecule stays connected."""
    em = Chem.RWMol(mol)
    for b in bonds:
        bd = mol.GetBondWithIdx(b)
        em.RemoveBond(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx())
    frags = Chem.GetMolFrags(em, asMols=False, sanitizeFrags=False)
    if len(frags) != 2:
        return None
    return tuple(float(masses[list(f)].sum()) for f in frags)


def cleavage_events(smiles):
    """[(event key, (piece mass 1, piece mass 2))]: acyclic single bonds + bond pairs within a ring."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() > 120:
        return []
    masses = _atom_masses(mol)
    out = []
    ring_bonds = set()
    for ring in mol.GetRingInfo().BondRings():
        ring_bonds.update(ring)
    for b in mol.GetBonds():
        if b.GetIdx() in ring_bonds or b.GetBondTypeAsDouble() > 1.5:
            continue
        pieces = _split(mol, [b.GetIdx()], masses)
        if pieces:
            out.append(("A:" + _bond_key(b), pieces))
    for ring in mol.GetRingInfo().BondRings():
        if len(ring) > 7 or any(mol.GetBondWithIdx(i).GetIsAromatic() for i in ring):
            continue
        for i in range(len(ring)):
            for j in range(i + 2, len(ring)):
                if i == 0 and j == len(ring) - 1:
                    continue
                pieces = _split(mol, [ring[i], ring[j]], masses)
                if pieces:
                    k = "|".join(sorted([_bond_key(mol.GetBondWithIdx(ring[i])), _bond_key(mol.GetBondWithIdx(ring[j]))]))
                    out.append(("R:" + k, pieces))
    return out


def _carriers(mode, adduct):
    if mode.lower().startswith("n") or adduct.endswith("-"):
        return [-PROTON]
    return [PROTON] + ([NA - ELECTRON] if "Na" in adduct else [])


def observed(events, mz, rel, mode, adduct, precursor_mz):
    """Per event: 1 if one of its pieces (±2 H, precursor carrier) matches a peak >= MIN_REL."""
    keep = (rel >= MIN_REL) & (mz < precursor_mz - 1.5)
    mz = np.sort(mz[keep])
    obs = np.zeros(len(events), np.int8)
    if len(mz) == 0:
        return obs
    shifts = np.array([c + k * H for c in _carriers(mode, adduct) for k in (-2, -1, 0, 1, 2)])
    for e, (_, pieces) in enumerate(events):
        ions = (np.array(pieces)[:, None] + shifts[None, :]).ravel()
        pos = np.searchsorted(mz, ions)
        lo = mz[np.clip(pos - 1, 0, len(mz) - 1)]
        hi = mz[np.clip(pos, 0, len(mz) - 1)]
        d = np.minimum(np.abs(ions - lo), np.abs(ions - hi))
        obs[e] = int((d <= np.maximum(TOL_DA, ions * TOL_PPM * 1e-6)).any())
    return obs


class BondModel:
    """Event-type propensities: P(observed | event key, ion mode), Beta-smoothed toward the global rate."""

    def __init__(self, counts=None, prior=4.0):
        self.counts = counts or {}          # (key, mode) -> [n, observed]
        self.prior = prior
        self._p = None

    def add(self, events, obs, mode):
        for (k, _), o in zip(events, obs):
            c = self.counts.setdefault((k, mode), [0, 0])
            c[0] += 1
            c[1] += int(o)

    def finalize(self):
        n = sum(c[0] for c in self.counts.values())
        o = sum(c[1] for c in self.counts.values())
        base = o / max(n, 1)
        self.base = base
        self._p = {k: (c[1] + self.prior * base) / (c[0] + self.prior) for k, c in self.counts.items()}
        return self

    def p(self, key, mode):
        return self._p.get((key, mode), self.base)

    def score(self, events, spectra):
        """(mean log-likelihood per event, observed share weighted by propensity) averaged over spectra."""
        if not events or not spectra:
            return 0.0, 0.0
        lls, ws = [], []
        for s in spectra:
            mode = "neg" if (str(s["mode"]).lower().startswith("n") or str(s["adduct"]).endswith("-")) else "pos"
            it = np.asarray(s["it"], float)
            rel = it / it.max() if len(it) and it.max() > 0 else it
            obs = observed(events, np.asarray(s["mz"], float), rel, mode, str(s["adduct"]), float(s["precursor_mz"]))
            p = np.clip(np.array([self.p(k, mode) for k, _ in events]), 1e-3, 1 - 1e-3)
            lls.append(float(np.mean(np.where(obs == 1, np.log(p), np.log(1 - p)))))
            ws.append(float((p * obs).sum() / p.sum()))
        return float(np.mean(lls)), float(np.mean(ws))
