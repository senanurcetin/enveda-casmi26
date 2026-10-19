"""Morgan fingerprints (radius 2, 2048 bits) for every pool structure, packed and row-aligned
with artifacts/pool.parquet. Used for analog propagation (Tanimoto candidate ↔ analog) and later
as the target of the spectrum → fingerprint model.

Output: artifacts/pool_fp2048.npy  (uint8, shape [n_pool, 256]), zero rows where RDKit fails
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE, POOL_TAG  # noqa: E402
from casmi.pool.fingerprint_provenance import make_fingerprint_provenance  # noqa: E402

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
    # Frozen scientific artifacts must never be silently replaced when
    # COCONUT normalization or the source pool has changed.
    manifest_path = ARTIFACT_DIR / f"pool_fp_provenance{POOL_TAG}.json"
    outputs = [ARTIFACT_DIR / POOL_FP_FILE, ARTIFACT_DIR / POOL_FP_KEYS_FILE,
               manifest_path]
    existing = [str(path.name) for path in outputs if os.path.lexists(path)]
    if existing:
        raise FileExistsError(
            "Refusing to replace existing fingerprints: " + ", ".join(existing)
            + "; use a new CASMI_POOL_TAG for an independent research variant"
        )
    # CandidatePool sorts by mass (stable) and drops duplicate keys: reproduce that order exactly
    pool = pd.read_parquet(ARTIFACT_DIR / POOL_FILE, columns=["pool_key", "smiles", "mass"])
    pool = pool.dropna(subset=["pool_key", "mass"]).sort_values("mass", kind="stable")
    pool = pool.drop_duplicates("pool_key", keep="first").reset_index(drop=True)
    if pool.empty:
        raise ValueError("Candidate pool is empty; cannot generate fingerprints")
    smiles = pool.smiles.tolist()
    chunks = [smiles[i:i + 5000] for i in range(0, len(smiles), 5000)]
    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 2) - 1)) as ex:
        fps = np.concatenate(list(ex.map(_fp_chunk, chunks)))
    provenance = make_fingerprint_provenance(pool, pool.pool_key.to_numpy(), fps)
    np.save(ARTIFACT_DIR / POOL_FP_FILE, fps)
    pool[["pool_key"]].to_parquet(ARTIFACT_DIR / POOL_FP_KEYS_FILE, index=False)
    with manifest_path.open("x", encoding="utf-8") as output:
        json.dump(provenance, output, indent=2)
        output.write("\n")
    empty = (fps.sum(axis=1) == 0).sum()
    print(f"{len(fps):,} fingerprints ({fps.nbytes / 1e6:.0f} MB), {empty:,} empty, built in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
