"""Class-3 hedge: propose unseen structures by applying one or two common natural-product
biotransformations to spectral analogs of the query.

A Class-3 molecule is in no database, so retrieval scores 0 on it by construction. Natural products
come in families, though: the analog search (D2) often finds a library structure whose spectrum
matches the query up to a precursor shift Δ. When Δ equals a common biotransformation (methylation,
hydroxylation, glycosylation, acylation, ...) the query is plausibly that analog with the
modification at one of its sites. Every site gives one candidate; the metric ignores stereo
(InChIKey block 1), so connectivity is all that has to be right.

`propose(analogs, query_mass)` returns generated SMILES with the analog they came from; scoring
and placement in the ranked list happen elsewhere (they need D3/D4 and a confidence gate).
"""
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt

RDLogger.DisableLog("rdApp.*")

HEX = "C1OC(CO)C(O)C(O)C1O"          # hexose (glucose/galactose connectivity)
DHEX = "C1OC(C)C(O)C(O)C1O"          # 6-deoxyhexose (rhamnose/fucose)
PENT = "C1OCC(O)C(O)C1O"             # pentose (pyranose)
GLCA = "C1OC(C(=O)O)C(O)C(O)C1O"     # glucuronic acid


@dataclass(frozen=True)
class Transform:
    name: str
    delta: float
    smarts: tuple


def _t(name, formula_delta, *smarts):
    return Transform(name, formula_delta, tuple(smarts))


# forward transforms; each is also applied in reverse where a removal SMARTS is given
TRANSFORMS = [
    _t("methyl", 14.01565, "[OX2H1:1]>>[O:1]C", "[NX3;H1,H2;!$(N-C=O):1]>>[N:1]C", "[cH1:1]>>[c:1]C"),
    _t("demethyl", -14.01565, "[O:1]-[CH3]>>[O:1]", "[N:1]-[CH3]>>[N:1]"),
    _t("hydroxyl", 15.99491, "[cH1:1]>>[c:1]O", "[CX4;H1,H2:1]>>[C:1]O"),
    _t("dehydroxyl", -15.99491, "[c:1]-[OH]>>[c:1]", "[CX4:1]-[OH]>>[C:1]"),
    _t("methoxyl", 30.01056, "[cH1:1]>>[c:1]OC"),
    _t("demethoxyl", -30.01056, "[c:1]-[OX2]-[CH3]>>[c:1]"),
    _t("hydrogenation", 2.01565, "[C:1]=[C:2]>>[C:1]-[C:2]"),
    _t("dehydrogenation", -2.01565, "[CX4;!H0:1]-[CX4;!H0:2]>>[C:1]=[C:2]", "[CX4;H1,H2:1]-[OH:2]>>[C:1]=[O:2]"),
    _t("hydration", 18.01056, "[C:1]=[C:2]>>[C:1](O)-[C:2]"),
    _t("dehydration", -18.01056, "[CX4;!H0:1]-[CX4:2]-[OH]>>[C:1]=[C:2]"),
    _t("hexose", 162.05282, f"[OX2H1:1]>>[O:1]{HEX}"),
    _t("dehexose", -162.05282, "[O:1]-[CH1;R1]1-[OX2;R1]-[CH1;R1](-[CH2]-[OH])-[CH1;R1](-[OH])-[CH1;R1](-[OH])-[CH1;R1]1-[OH]>>[O:1]"),
    _t("deoxyhexose", 146.05791, f"[OX2H1:1]>>[O:1]{DHEX}"),
    _t("dedeoxyhexose", -146.05791, "[O:1]-[CH1;R1]1-[OX2;R1]-[CH1;R1](-[CH3])-[CH1;R1](-[OH])-[CH1;R1](-[OH])-[CH1;R1]1-[OH]>>[O:1]"),
    _t("pentose", 132.04226, f"[OX2H1:1]>>[O:1]{PENT}"),
    _t("glucuronide", 176.03209, f"[OX2H1:1]>>[O:1]{GLCA}"),
    _t("acetyl", 42.01057, "[OX2H1:1]>>[O:1]C(C)=O"),
    _t("deacetyl", -42.01057, "[O:1]-C(=O)-[CH3]>>[O:1]"),
    _t("malonyl", 86.00039, "[OX2H1:1]>>[O:1]C(=O)CC(=O)O"),
    _t("sulfate", 79.95682, "[OX2H1:1]>>[O:1]S(=O)(=O)O"),
    _t("prenyl", 68.06260, "[cH1:1]>>[c:1]CC=C(C)C", "[OX2H1:1]>>[O:1]CC=C(C)C"),
    _t("coumaroyl", 146.03678, "[OX2H1:1]>>[O:1]C(=O)C=Cc1ccc(O)cc1"),
    _t("caffeoyl", 162.03169, "[OX2H1:1]>>[O:1]C(=O)C=Cc1ccc(O)c(O)c1"),
    _t("feruloyl", 176.04734, "[OX2H1:1]>>[O:1]C(=O)C=Cc1ccc(O)c(OC)c1"),
    _t("galloyl", 152.01096, "[OX2H1:1]>>[O:1]C(=O)c1cc(O)c(O)c(O)c1"),
]
_RXN = {t.name: [AllChem.ReactionFromSmarts(s) for s in t.smarts] for t in TRANSFORMS}
DELTAS = np.array([t.delta for t in TRANSFORMS])


def matching(delta, tol_da):
    """Single transforms and ordered pairs whose total mass change equals `delta`."""
    single = [(t,) for t in TRANSFORMS if abs(t.delta - delta) <= tol_da]
    pairs = [(a, b) for i, a in enumerate(TRANSFORMS) for b in TRANSFORMS[i:]
             if abs(a.delta + b.delta - delta) <= tol_da]
    return single + pairs


def apply(mol, transform, cap=60):
    """Every single-site product of `transform` on `mol` (sanitised, unique by SMILES)."""
    out = {}
    for rxn in _RXN[transform.name]:
        try:
            products = rxn.RunReactants((mol,), maxProducts=cap)
        except Exception:
            continue
        for (p,) in products:
            try:
                Chem.SanitizeMol(p)
                smi = Chem.MolToSmiles(p)
            except Exception:
                continue
            if smi and "." not in smi:
                out.setdefault(smi, p)
        if len(out) >= cap:
            break
    return out


def propose(analogs, query_mass, tol_ppm=10.0, max_analogs=8, cap=200):
    """analogs: [(smiles, similarity, analog_mass)], best first.
    Returns [(smiles, similarity, transform names)] of structures at the query mass."""
    tol_da = query_mass * tol_ppm * 1e-6
    seen, out = set(), []
    for smi, sim, amass in analogs[:max_analogs]:
        paths = matching(query_mass - amass, tol_da * 2)
        if not paths:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        for path in paths:
            level = {smi: mol}
            for t in path:
                nxt = {}
                for m in level.values():
                    nxt.update(apply(m, t))
                level = nxt
                if not level:
                    break
            for psmi, p in level.items():
                if psmi in seen:
                    continue
                seen.add(psmi)
                if abs(ExactMolWt(p) - query_mass) <= tol_da:
                    out.append((psmi, sim, "+".join(t.name for t in path)))
            if len(out) >= cap:
                return out
    return out
