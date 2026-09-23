"""A3 — candidate pool = training structures ∪ COCONUT 2.0, keyed on standard InChIKey14.

COCONUT (CC0, coconut.naturalproducts.net, 09-2026 lite CSV) contributes natural-product
structures that have no reference spectra — exactly the Class-2 search space. Per-structure
COCONUT annotations (NP-likeness, NPClassifier pathway, sugars) are kept as ranker priors.

Output: artifacts/pool.parquet
"""
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.chem.masses import formula_mass  # noqa: E402
from casmi.pool.pool import SOURCE_BITS  # noqa: E402

COCONUT_ZIP = Path("C:/Users/MONSTER/kaggle-data/external/coconut/coconut_csv_lite-09-2026.zip")
MASS_RANGE = (100.0, 1300.0)          # test monoisotopic masses span 157–1,159 Da
ENVEDA180_BIT = 1                     # bit 0 in build_structures.LIBS


def largest_fragment(smiles):
    from rdkit import Chem
    frags = smiles.split(".")
    return max(frags, key=lambda f: (Chem.MolFromSmiles(f) or Chem.MolFromSmiles("C")).GetNumHeavyAtoms())


def load_coconut():
    cols = ["identifier", "canonical_smiles", "standard_inchi_key", "molecular_formula",
            "exact_molecular_weight", "formal_charge", "np_likeness", "np_classifier_pathway",
            "np_classifier_superclass", "contains_sugar"]
    with zipfile.ZipFile(COCONUT_ZIP).open("coconut_csv_lite-09-2026.csv") as f:
        c = pd.read_csv(f, usecols=cols, low_memory=False)
    n0 = len(c)
    c = c[(c.formal_charge == 0) & c.canonical_smiles.notna() & c.standard_inchi_key.notna()].copy()
    multi = c.canonical_smiles.str.contains(".", regex=False)
    # salts / counter-ions: keep the parent. Its formula no longer matches, so mass is recomputed below.
    c.loc[multi, "canonical_smiles"] = c.loc[multi, "canonical_smiles"].map(largest_fragment)
    c["mass"] = c.molecular_formula.fillna("").map(formula_mass)
    if multi.any():
        from rdkit.Chem.Descriptors import ExactMolWt
        from rdkit import Chem
        c.loc[multi, "mass"] = c.loc[multi, "canonical_smiles"].map(lambda s: ExactMolWt(Chem.MolFromSmiles(s)))
    c["pool_key"] = c.standard_inchi_key.str[:14]
    c = c[c.mass.between(*MASS_RANGE)].drop_duplicates("pool_key")
    print(f"COCONUT: {n0:,} rows -> {len(c):,} neutral structures in {MASS_RANGE} Da ({multi.sum():,} de-salted)")
    return c.rename(columns={"canonical_smiles": "smiles"})


def main():
    t0 = time.time()
    s = pd.read_parquet(ARTIFACT_DIR / "structures.parquet")
    train = pd.DataFrame({"pool_key": s.pool_key, "smiles": s.smiles, "mass": s.mass,
                          "sources": SOURCE_BITS["train"] | np.where(s.libs == ENVEDA180_BIT, SOURCE_BITS["enveda180_only"], 0),
                          "n_spectra": s.n_spectra, "libs": s.libs})
    coco = load_coconut()
    coco_part = pd.DataFrame({"pool_key": coco.pool_key, "smiles": coco.smiles, "mass": coco.mass,
                              "sources": SOURCE_BITS["coconut"], "n_spectra": 0, "libs": 0,
                              "np_likeness": coco.np_likeness, "np_pathway": coco.np_classifier_pathway,
                              "np_superclass": coco.np_classifier_superclass, "has_sugar": coco.contains_sugar})
    both = pd.concat([train, coco_part], ignore_index=True)
    g = both.groupby("pool_key", sort=False)
    pool = pd.DataFrame({
        "smiles": g.smiles.first(),                       # train SMILES win: they are the labels' form
        "mass": g.mass.first(),
        "sources": g.sources.agg(np.bitwise_or.reduce),
        "n_spectra": g.n_spectra.sum(), "libs": g.libs.agg(np.bitwise_or.reduce),
        "np_likeness": g.np_likeness.first(), "np_pathway": g.np_pathway.first(),
        "np_superclass": g.np_superclass.first(), "has_sugar": g.has_sugar.first(),
    }).reset_index()
    pool = pool.dropna(subset=["mass"]).sort_values("mass", kind="stable").reset_index(drop=True)
    pool.to_parquet(ARTIFACT_DIR / "pool.parquet", index=False)

    in_train = (pool.sources & SOURCE_BITS["train"]) > 0
    in_coco = (pool.sources & SOURCE_BITS["coconut"]) > 0
    print(f"pool: {len(pool):,} structures | train {in_train.sum():,} | COCONUT {in_coco.sum():,} | "
          f"both {(in_train & in_coco).sum():,} | COCONUT-only {(in_coco & ~in_train).sum():,}")
    np_ex = (pool.libs & (1 << 9)) > 0
    print(f"enveda-np-examples structures covered by COCONUT: {(np_ex & in_coco).sum()}/{np_ex.sum()}")
    print(f"built in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
