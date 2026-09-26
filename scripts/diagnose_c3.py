"""Where does the Class-3 hedge miss? For C3-simulated queries (V-in + V-real COCONUT half):
  - hit rate of the current transforms when seeding from the top-K spectral analogs (K = 8, 15, 30)
  - for the misses: the formula difference between the truth and its structurally closest analog
    among the top 30 (Tanimoto), i.e. the transformation that would have been needed
"""
import importlib.util
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdMolDescriptors  # noqa: E402

from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.channels.analog import AnalogIndex, tanimoto_matrix  # noqa: E402
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


def formula_delta(a, b):
    ca = Counter(dict((el, int(n or 1)) for el, n in __import__("re").findall(r"([A-Z][a-z]?)(\d*)", a)))
    cb = Counter(dict((el, int(n or 1)) for el, n in __import__("re").findall(r"([A-Z][a-z]?)(\d*)", b)))
    d = {el: ca[el] - cb[el] for el in set(ca) | set(cb) if ca[el] != cb[el]}
    return " ".join(f"{el}{d[el]:+d}" for el in sorted(d)) or "isomer"


def main(n_max=400):
    cfg = EngineConfig()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet"))
    fps = np.load(ARTIFACT_DIR / "pool_fp2048_pcl.npy")
    library = LibraryIndex(store, key_of_row, cfg.library)
    analog = AnalogIndex(store, key_of_row, cfg.analog, cache=ARTIFACT_DIR / "analog_index.npz")
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")
    smiles = brd.structure_smiles()
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    vreal, stratum_of = brd.vreal_queries(store, key_of_row, vin_skel)
    items = [(q, "vin") for q in vin] + [(q, "coco") for q, _ in vreal if stratum_of[q.pool_key] == "coco"]
    items = items[:n_max]
    hits = {8: 0, 15: 0, 30: 0}
    n_props = []
    needed = Counter()
    for q, st in items:
        mass = molecule_mass(q.spectra, cfg.mass_offset_ppm)
        truth, same = equivalent_keys(q.pool_key, mass, pool, keys, 30, truth_smiles=smiles.get(q.pool_key))
        ko = knockout_for(q, "C3", library.rows_of, same | {q.pool_key})
        spectra = [dict(s) for s in q.spectra]
        for s in spectra:
            s["mz_c"], s["it_c"] = clean_peaks_fast(s["mz"], s["it"], s["precursor_mz"], library.cfg.clean)
        found = [(k, v) for k, v in analog.search(spectra, mass, ko.excluded_rows).items()
                 if k not in same and k != q.pool_key and k in pool.row_of]
        tskel = skeleton_key(smiles[q.pool_key])
        got = False
        for k_ in (8, 15, 30):
            if not got or k_ == 15:
                props = propose([(pool.smiles_of(k), s, m) for k, (s, m) in found], mass, max_analogs=k_, cap=400)
                if k_ == 15:
                    n_props.append(len(props))
                got = any(skeleton_key(p) == tskel and canonical_ik14(p) == truth for p, _, _ in props)
            hits[k_] += got
        if not got and found:
            top = [k for k, _ in found[:30]]
            T = tanimoto_matrix(fps[[pool.row_of[q.pool_key]]] if q.pool_key in pool.row_of else fps[:1],
                                fps[[pool.row_of[k] for k in top]])[0]
            best = top[int(T.argmax())]
            ft = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles[q.pool_key]))
            fa = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(pool.smiles_of(best)))
            needed[(formula_delta(ft, fa), round(float(T.max()), 1))] += 1
    n = len(items)
    print(f"median / mean proposals at top-15: {np.median(n_props):.0f} / {np.mean(n_props):.1f}")
    print(f"{n} queries | truth proposed with top-8 {hits[8] / n:.3f}, top-15 {hits[15] / n:.3f}, top-30 {hits[30] / n:.3f}")
    print("misses: formula change truth - closest analog (Tanimoto), most common:")
    for (d, t), c in needed.most_common(40):
        print(f"  {c:4d}  {d:<28} T~{t}")
    keys.save()


if __name__ == "__main__":
    main()
