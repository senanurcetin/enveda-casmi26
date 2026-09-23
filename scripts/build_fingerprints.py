"""Morgan fingerprints (radius 2, 2048 bits) for every pool structure, packed and row-aligned
with artifacts/pool.parquet. Used for analog propagation (Tanimoto candidate ↔ analog) and later
as the target of the spectrum → fingerprint model.

Output: artifacts/pool_fp2048.npy  (uint8, shape [n_pool, 256]), zero rows where RDKit fails
"""
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402

NBITS = 2048


def _fp_chunk(smiles_list):
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=NBITS)
    out = np.zeros((len(smiles_list), NBITS // 8), np.uint8)
    arr = np.zeros(NBITS, np.uint8)
    for i, s in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if mol is None:
            continue
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
        out[i] = np.packbits(arr)
    return out


def main():
    t0 = time.time()
    # CandidatePool sorts by mass (stable) and drops duplicate keys: reproduce that order exactly
    pool = pd.read_parquet(ARTIFACT_DIR / "pool.parquet", columns=["pool_key", "smiles", "mass"])
    pool = pool.dropna(subset=["pool_key", "mass"]).sort_values("mass", kind="stable")
    pool = pool.drop_duplicates("pool_key", keep="first").reset_index(drop=True)
    smiles = pool.smiles.tolist()
    chunks = [smiles[i:i + 5000] for i in range(0, len(smiles), 5000)]
    with ProcessPoolExecutor(max_workers=6) as ex:
        fps = np.concatenate(list(ex.map(_fp_chunk, chunks)))
    np.save(ARTIFACT_DIR / "pool_fp2048.npy", fps)
    pool[["pool_key"]].to_parquet(ARTIFACT_DIR / "pool_fp_keys.parquet", index=False)
    empty = (fps.sum(axis=1) == 0).sum()
    print(f"{len(fps):,} fingerprints ({fps.nbytes / 1e6:.0f} MB), {empty:,} empty, built in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
