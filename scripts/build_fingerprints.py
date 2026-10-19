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


def _fp_chunk(job):
    """Produce packed Morgan bits, or fail on the first invalid candidate.

    Do not silently turn a malformed/empty SMILES into an all-zero vector.
    Such rows otherwise become apparently valid signed research features.
    Report only the row number: raw licensed molecule strings stay private.
    """
    start, smiles_list = job
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=NBITS)
    out = np.zeros((len(smiles_list), NBITS // 8), np.uint8)
    arr = np.zeros(NBITS, np.uint8)
    for i, s in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(s) if isinstance(s, str) else None
            if mol is None or mol.GetNumHeavyAtoms() == 0:
                raise ValueError("Invalid candidate SMILES")
            arr.fill(0)
            DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
            if not np.any(arr):
                raise ValueError("Zero-vector Morgan fingerprint")
        except (ValueError, TypeError, RuntimeError) as exc:
            raise ValueError(
                f"Cannot create a valid Morgan fingerprint for pool row {start + i}"
            ) from None
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
    # Use the actual inference pool constructor rather than maintaining a
    # second sorting/deduplication implementation that can drift.
    from casmi.pool.pool import CandidatePool
    source = pd.read_parquet(
        ARTIFACT_DIR / POOL_FILE, columns=["pool_key", "smiles", "mass"]
    )
    pool = CandidatePool(source).table
    if pool.empty:
        raise ValueError("Candidate pool is empty; cannot generate fingerprints")
    smiles = pool.smiles.tolist()
    chunks = [(i, smiles[i:i + 5000]) for i in range(0, len(smiles), 5000)]
    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 2) - 1)) as ex:
        fps = np.concatenate(list(ex.map(_fp_chunk, chunks)))
    if not np.all(fps.any(axis=1)):
        raise ValueError("Zero-vector candidate fingerprints are prohibited")
    provenance = make_fingerprint_provenance(pool, pool.pool_key.to_numpy(), fps)
    # Exclusive output files close the race between the initial existence
    # check and publication, protecting frozen historical arrays.
    with (ARTIFACT_DIR / POOL_FP_FILE).open("xb") as output:
        np.save(output, fps)
    with (ARTIFACT_DIR / POOL_FP_KEYS_FILE).open("xb") as output:
        pool[["pool_key"]].to_parquet(output, index=False)
    with manifest_path.open("x", encoding="utf-8") as output:
        json.dump(provenance, output, indent=2)
        output.write("\n")
    empty = (fps.sum(axis=1) == 0).sum()
    print(f"{len(fps):,} fingerprints ({fps.nbytes / 1e6:.0f} MB), {empty:,} empty, built in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
