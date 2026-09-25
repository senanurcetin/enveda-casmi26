"""Extra structure sources for the candidate pool, normalised to pool rows (pool_key, smiles, mass).

  LIPID MAPS LMSD ids/SMILES (lipidmaps_ids_cc0.tsv, CC0)       -> artifacts/extra_lipidmaps.parquet
  ChEBI structure registry (structure_registry.tsv.gz, CC BY 4.0) -> artifacts/extra_chebi.parquet

Keys are recomputed with RDKit (standard InChIKey, first block) so every source is keyed exactly
like the rest of the pool. Filters: single component, *net* charge 0 (zwitterionic phospholipids
are kept), all elements tabulated, monoisotopic mass within the test range.
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402

EXT = Path("C:/Users/MONSTER/kaggle-data/external")
MASS_RANGE = (100.0, 1300.0)


def _rows(smiles_list):
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Descriptors import ExactMolWt
    RDLogger.DisableLog("rdApp.*")
    out = []
    for s in smiles_list:
        if not isinstance(s, str) or "." in s:
            continue
        m = Chem.MolFromSmiles(s)
        if m is None or Chem.GetFormalCharge(m) != 0:
            continue
        mass = ExactMolWt(m)
        if not (MASS_RANGE[0] <= mass <= MASS_RANGE[1]):
            continue
        key = Chem.MolToInchiKey(m)
        if key:
            out.append((key[:14], Chem.MolToSmiles(m), mass))
    return out


def normalise(smiles, name):
    t0 = time.time()
    chunks = [smiles[i:i + 4000] for i in range(0, len(smiles), 4000)]
    with ProcessPoolExecutor(max_workers=6) as ex:
        rows = [r for chunk in ex.map(_rows, chunks) for r in chunk]
    df = pd.DataFrame(rows, columns=["pool_key", "smiles", "mass"]).drop_duplicates("pool_key")
    df.to_parquet(ARTIFACT_DIR / f"extra_{name}.parquet", index=False)
    print(f"{name}: {len(smiles):,} input -> {len(df):,} structures ({time.time() - t0:.0f}s)")


def main():
    if not (ARTIFACT_DIR / "extra_lipidmaps.parquet").exists():
        lm = pd.read_csv(EXT / "lipidmaps" / "lipidmaps_ids_cc0.tsv", sep="\t", usecols=["smiles"])
        normalise(lm.smiles.tolist(), "lipidmaps")
    reg = pd.read_csv(EXT / "chebi" / "structure_registry.tsv.gz", sep="\t", usecols=["layers"])
    smiles = [json.loads(x).get("CANONICAL_SMILES") for x in reg.layers if isinstance(x, str)]   # some rows have no structure
    normalise(smiles, "chebi")


if __name__ == "__main__":
    main()
