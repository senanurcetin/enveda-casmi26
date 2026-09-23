"""Training tensors for the spectrum → fingerprint model (D3), precomputed once.

Cleaning every spectrum inside the training loop would make the GPU wait on Python, so each
eligible spectrum is cleaned once to ≤ 64 peaks and written to fixed-shape memmaps.

Output: artifacts/fpdata/
  mz.f32 [N, 64], it.f16 [N, 64], npk.i8 [N]           peaks (zero-padded)
  side.parquet   precursor_mz, adduct_id, ce_ev, mode, struct_idx, fold, is_timstof
  struct_fp.npy  packed Morgan-2048 of every referenced structure (same bits as the pool)
Folds come from the skeleton hash used for the ranker data (tautomers share a fold); fold = -1
marks V-in structures, which no fold model may train on.
"""
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.spectra.clean import clean_peaks_fast  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402

from casmi.channels.fpinfer import ADDUCT_ID, FP_CLEAN, MAX_PEAKS  # noqa: E402  (one definition, shared with inference)


def main():
    t0 = time.time()
    out = ARTIFACT_DIR / "fpdata"
    out.mkdir(exist_ok=True)
    store = SpectrumStore(ARTIFACT_DIR / "store")
    meta = store.meta
    ok = (meta.test_adduct.to_numpy() & (np.abs(meta.ppm.to_numpy()) <= 20) & (meta.n_peaks.to_numpy() >= 3))
    rows = np.flatnonzero(ok)
    keys = meta.train_ik14.to_numpy(dtype=object)[rows]
    print(f"{len(rows):,} eligible spectra", flush=True)

    # structures: fingerprint row in the pool, skeleton fold, V-in flag
    pool_keys = pd.read_parquet(ARTIFACT_DIR / "pool_fp_keys.parquet").pool_key.to_numpy()
    fp_row = dict(zip(pool_keys, range(len(pool_keys))))
    pool_fps = np.load(ARTIFACT_DIR / "pool_fp2048.npy", mmap_mode="r")
    smiles = pd.read_parquet(ARTIFACT_DIR / "structures.parquet", columns=["pool_key", "smiles"])
    smiles = dict(zip(smiles.pool_key, smiles.smiles))
    uniq = pd.unique(keys)
    skel = {k: skeleton_key(smiles[k]) for k in uniq}
    vin_keys = set(meta.train_ik14[meta.ingest_lib == "enveda-np-examples"])
    vin_skel = {skel[k] for k in vin_keys if k in skel}
    fold_of = {k: (-1 if (s is None or s in vin_skel) else int(hashlib.md5(s.encode()).hexdigest(), 16) % 2)
               for k, s in skel.items()}
    struct_idx = {k: i for i, k in enumerate(uniq)}
    np.save(out / "struct_fp.npy", np.asarray(pool_fps[[fp_row[k] for k in uniq]]))
    pd.DataFrame({"pool_key": uniq}).to_parquet(out / "structs.parquet", index=False)
    print(f"{len(uniq):,} structures, skeletons in {time.time() - t0:.0f}s", flush=True)

    n = len(rows)
    mz = np.lib.format.open_memmap(out / "mz.npy", "w+", np.float32, (n, MAX_PEAKS))
    it = np.lib.format.open_memmap(out / "it.npy", "w+", np.float16, (n, MAX_PEAKS))
    npk = np.zeros(n, np.int8)
    prec = meta.precursor_mz.to_numpy()
    for j, r in enumerate(rows):
        m, i = store.peaks(r)
        cm, ci = clean_peaks_fast(m, i, prec[r], FP_CLEAN)
        k = cm.size
        mz[j, :k], it[j, :k], npk[j] = cm, ci, k
        if (j + 1) % 500_000 == 0:
            print(f"  {j + 1:,}/{n:,} cleaned ({time.time() - t0:.0f}s)", flush=True)
    mz.flush(); it.flush()
    np.save(out / "npk.npy", npk)

    side = pd.DataFrame({
        "store_row": rows,
        "precursor_mz": prec[rows].astype(np.float32),
        "adduct_id": meta.adduct.astype(str).to_numpy()[rows],
        "ce_ev": meta.ce_ev.to_numpy()[rows],
        "positive": (meta.ionization_mode.astype(str).to_numpy()[rows] == "positive").astype(np.int8),
        "is_timstof": meta.is_timstof.to_numpy()[rows].astype(np.int8),
        "struct_idx": np.array([struct_idx[k] for k in keys], np.int32),
        "fold": np.array([fold_of[k] for k in keys], np.int8),
    })
    side["adduct_id"] = side.adduct_id.map(ADDUCT_ID).astype(np.int8)
    side.to_parquet(out / "side.parquet", index=False)
    print(side.fold.value_counts().rename({-1: "V-in (held out)", 0: "fold 0", 1: "fold 1"}).to_string())
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
