"""A2 — structure table: one row per training structure.

Output: artifacts/structures.parquet
  pool_key    inchikey14 as shipped in train.parquet (standard InChIKey14; also the pool identity)
  smiles, formula, mass (monoisotopic, from formula)
  libs        bitmask of ingest_lib values the structure appears in
  n_spectra

The metric key (tautomer-canonical InChIKey14) is computed lazily by casmi.chem.keys.MetricKeys.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, DATA_DIR  # noqa: E402
from casmi.chem.masses import formula_mass  # noqa: E402

LIBS = ["enveda-180", "pluskal_ms2", "riken", "gnps", "massbank", "mona", "spectraverse",
        "msdial", "drug_plus", "enveda-np-examples", "masaryk"]
LIB_BIT = {lib: 1 << i for i, lib in enumerate(LIBS)}


def main():
    t0 = time.time()
    pf = pq.ParquetFile(DATA_DIR / "train.parquet")
    frames = []
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=["inchikey14", "normalized_smiles", "molecular_formula", "ingest_lib"]).to_pandas()
        t["bit"] = t.ingest_lib.map(LIB_BIT).astype("int64")
        g = t.groupby("inchikey14")
        frames.append(pd.DataFrame({"smiles": g.normalized_smiles.first(), "formula": g.molecular_formula.first(),
                                    "libs": g.bit.agg(np.bitwise_or.reduce), "n_spectra": g.size()}))
    s = pd.concat(frames)
    g = s.groupby(level=0)
    s = pd.DataFrame({"smiles": g.smiles.first(), "formula": g.formula.first(),
                      "libs": g.libs.agg(np.bitwise_or.reduce), "n_spectra": g.n_spectra.sum()})
    s.index.name = "pool_key"
    s = s.reset_index()
    s["mass"] = s.formula.map(formula_mass)
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    s.to_parquet(ARTIFACT_DIR / "structures.parquet", index=False)
    print(f"{len(s):,} structures, mass missing {s.mass.isna().sum()}, built in {time.time() - t0:.0f}s")
    for lib, bit in LIB_BIT.items():
        print(f"  {lib:<20} {((s.libs & bit) > 0).sum():>7,}")


if __name__ == "__main__":
    main()
