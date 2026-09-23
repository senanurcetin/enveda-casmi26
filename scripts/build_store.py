"""Build the on-disk spectrum store (A1) from train.parquet, one row group at a time.

Peaks are read through Arrow's flat value/offset buffers instead of Python lists, which keeps
peak memory at a few hundred MB per row group. Masses and ppm errors are recomputed from the
formula and adduct, because the shipped precursor_error_ppm column follows different conventions
in different libraries (see reports/).
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, DATA_DIR  # noqa: E402
from casmi.chem.masses import TEST_ADDUCTS, formula_mass, neutral_mass, parse_adduct  # noqa: E402
from casmi.spectra.clean import CleanConfig, clean_peaks_fast  # noqa: E402

STORE_CFG = CleanConfig(min_rel_intensity=0.001, max_peaks=150, top_k_per_window=0)
META_COLS = ["ingest_lib", "inchikey14", "molecular_formula", "ionization_mode", "instrument_type",
             "adduct", "precursor_mz", "base_peak_intensity", "collision_energy_orig_units"]


def _flat(col):
    """(values, offsets) of a list column as NumPy, without materialising Python lists."""
    arr = col.combine_chunks() if isinstance(col, pa.ChunkedArray) else col
    return arr.values.to_numpy(zero_copy_only=False), arr.offsets.to_numpy()


def main():
    out = ARTIFACT_DIR / "store"
    out.mkdir(parents=True, exist_ok=True)
    pf = pq.ParquetFile(DATA_DIR / "train.parquet")
    f_mz, f_it = open(out / "mz.f32", "wb"), open(out / "it.f32", "wb")
    offsets, metas, n_peaks_total = [0], [], 0
    t0 = time.time()

    for rg in range(pf.num_row_groups):
        tbl = pf.read_row_group(rg, columns=META_COLS + ["ms2_mzs", "ms2_normalized_intensities", "collision_energy_ev"])
        meta = tbl.select(META_COLS).to_pandas()
        mz_v, mz_o = _flat(tbl.column("ms2_mzs"))
        it_v, it_o = _flat(tbl.column("ms2_normalized_intensities"))
        ce_v, ce_o = _flat(tbl.column("collision_energy_ev"))
        prec = meta.precursor_mz.to_numpy()

        buf_mz, buf_it, lens = [], [], np.zeros(len(meta), np.int32)
        for i in range(len(meta)):
            cm, ci = clean_peaks_fast(mz_v[mz_o[i]:mz_o[i + 1]], it_v[it_o[i]:it_o[i + 1]], prec[i], STORE_CFG)
            buf_mz.append(cm)
            buf_it.append(ci)
            lens[i] = cm.size
        np.concatenate(buf_mz).tofile(f_mz)
        np.concatenate(buf_it).tofile(f_it)
        offsets.extend((n_peaks_total + np.cumsum(lens)).tolist())
        n_peaks_total += int(lens.sum())

        # collision energy: mean of the listed eV values, count of merged acquisitions
        # (bincount over segment ids, not np.add.reduceat, which mishandles empty trailing lists)
        ce_n = np.diff(ce_o)
        seg = np.repeat(np.arange(len(meta)), ce_n)
        ce_sum = np.bincount(seg, weights=np.nan_to_num(ce_v[ce_o[0]:ce_o[-1]]), minlength=len(meta))
        ce_sum = np.where(ce_n > 0, ce_sum, np.nan)
        meta["ce_ev"] = (ce_sum / np.maximum(ce_n, 1)).astype(np.float32)
        meta["ce_count"] = ce_n.astype(np.int8)
        meta["n_peaks"] = lens
        meta["row_group"] = np.int16(rg)
        metas.append(meta)
        print(f"row group {rg + 1}/{pf.num_row_groups}: {len(meta):,} spectra, {lens.mean():.0f} peaks avg, "
              f"{time.time() - t0:.0f}s", flush=True)
        del tbl, buf_mz, buf_it

    f_mz.close()
    f_it.close()
    np.asarray(offsets, dtype=np.int64).tofile(out / "offsets.i64")

    meta = pd.concat(metas, ignore_index=True)
    meta = meta.rename(columns={"inchikey14": "train_ik14", "molecular_formula": "formula"})
    meta["is_timstof"] = meta.ingest_lib.isin(["enveda-180", "enveda-np-examples"])
    meta["test_adduct"] = meta.adduct.isin(TEST_ADDUCTS)

    # our own masses: neutral mass implied by precursor + adduct, and signed ppm vs the formula
    meta["formula_mass"] = meta.formula.map(formula_mass).astype(np.float64)
    nm = np.full(len(meta), np.nan)
    for adduct, idx in meta.groupby("adduct").indices.items():
        if parse_adduct(adduct) is not None:
            nm[idx] = neutral_mass(meta.precursor_mz.to_numpy()[idx], adduct)
    meta["neutral_mass"] = nm
    meta["ppm"] = ((nm - meta.formula_mass) / meta.formula_mass * 1e6).astype(np.float32)

    for c in ["ingest_lib", "ionization_mode", "instrument_type", "adduct", "collision_energy_orig_units"]:
        meta[c] = meta[c].astype("category")
    meta.to_parquet(out / "meta.parquet", index=False)
    print(f"\nstore: {len(meta):,} spectra, {n_peaks_total:,} peaks "
          f"({n_peaks_total * 8 / 1e9:.2f} GB), built in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
