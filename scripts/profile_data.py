"""Memory-safe profile of the CASMI 2026 data.

Streams train.parquet one row group at a time and reads only the columns each
question needs, so it runs with ~2 GB of free RAM. Writes reports/data_profile.json.

Usage:  python scripts/profile_data.py            (data dir from $CASMI_DATA)
"""
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pyarrow.compute as pc
import pyarrow.parquet as pq

DATA = Path(os.environ.get("CASMI_DATA", "C:/Users/MONSTER/kaggle-data/casmi26"))
OUT = Path(__file__).resolve().parents[1] / "reports" / "data_profile.json"

TEST_ADDUCTS = {"[M+H]+", "[M+NH4]+", "[M-H2O+H]+", "[M-2H2O+H]+", "[M+Na]+", "[M+K]+",
                "[M-H]-", "[M-H2O-H]-", "[M+CH2O2-H]-", "[M+Cl]-"}
TIMSTOF_LIBS = {"enveda-180", "enveda-np-examples"}


def q(values, qs=(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return {f"p{int(x * 100)}": round(float(np.quantile(v, x)), 3) for x in qs} if len(v) else {}


# ------------------------------------------------------------------ test file
test = pq.read_table(DATA / "test.parquet").to_pandas()
per_mol = test.groupby("molecule_id").size()
test_profile = {
    "spectra": len(test),
    "molecules": int(test.molecule_id.nunique()),
    "spectra_per_molecule": q(per_mol.values),
    "adducts": test.adduct.value_counts().to_dict(),
    "ionization_mode": test.ionization_mode.value_counts().to_dict(),
    "instrument_type": test.instrument_type.value_counts().to_dict(),
    "precursor_mz": q(test.precursor_mz),
    "num_peaks": q(test.ms2_mzs.map(len)),
    "ce_values_per_spectrum": test.collision_energy_ev.map(lambda x: 0 if x is None else len(x)).value_counts().to_dict(),
    "base_peak_intensity": q(test.base_peak_intensity),
}
# fingerprint used to find the public test spectra inside train (it is sampled from train)
test_keys = {(round(r.precursor_mz, 4), len(r.ms2_mzs), round(float(r.ms2_mzs[0]), 4)): r.spectrum_id
             for r in test.itertuples()}
test_prec = {k[0] for k in test_keys}

# ------------------------------------------------------------------ train, streamed
pf = pq.ParquetFile(DATA / "train.parquet")
LIGHT = ["ingest_lib", "inchikey14", "ionization_mode", "instrument_type", "adduct",
         "precursor_mz", "precursor_error_ppm", "num_peaks", "base_peak_intensity",
         "collision_energy_orig_units"]

lib_spectra, lib_structs = Counter(), defaultdict(set)
adducts, adducts_timstof, modes, units = Counter(), Counter(), Counter(), Counter()
struct_spectra = Counter()
ppm_by_lib, peaks_by_lib, prec_timstof = defaultdict(list), defaultdict(list), []
n_test_adduct = 0
matches = {}

for rg in range(pf.num_row_groups):
    t = pf.read_row_group(rg, columns=LIGHT).to_pandas()
    lib_spectra.update(t.ingest_lib.value_counts().to_dict())
    for lib, keys in t.groupby("ingest_lib").inchikey14:
        lib_structs[lib].update(keys.unique())
    struct_spectra.update(t.inchikey14.value_counts().to_dict())
    adducts.update(t.adduct.value_counts().to_dict())
    adducts_timstof.update(t[t.ingest_lib.isin(TIMSTOF_LIBS)].adduct.value_counts().to_dict())
    modes.update(t.ionization_mode.fillna("null").value_counts().to_dict())
    units.update(t.collision_energy_orig_units.fillna("null").value_counts().to_dict())
    n_test_adduct += int(t.adduct.isin(TEST_ADDUCTS).sum())
    rng = np.random.default_rng(rg)
    for lib, g in t.groupby("ingest_lib"):
        take = g.sample(min(len(g), 4000), random_state=rg)
        ppm_by_lib[lib].extend(take.precursor_error_ppm.tolist())
        peaks_by_lib[lib].extend(take.num_peaks.tolist())
    prec_timstof.extend(t.loc[t.ingest_lib.isin(TIMSTOF_LIBS), "precursor_mz"].sample(frac=.05, random_state=rg).tolist())

    # locate public-test spectra: filter on precursor first, read peaks only for those rows
    cand = t.index[t.precursor_mz.round(4).isin(test_prec)]
    if len(cand):
        peaks = pf.read_row_group(rg, columns=["ms2_mzs"]).column(0)
        for i in cand:
            mzs = peaks[int(i)].as_py()
            key = (round(t.precursor_mz.iat[i], 4), len(mzs), round(float(mzs[0]), 4))
            if key in test_keys:
                matches[test_keys[key]] = (t.ingest_lib.iat[i], t.inchikey14.iat[i])
    print(f"row group {rg + 1}/{pf.num_row_groups} done", flush=True)
    del t

# ------------------------------------------------------------------ overlap questions
np_ex = lib_structs["enveda-np-examples"]
other = set().union(*(s for lib, s in lib_structs.items() if lib != "enveda-np-examples"))
non_timstof = set().union(*(s for lib, s in lib_structs.items() if lib not in TIMSTOF_LIBS))
e180 = lib_structs["enveda-180"]

matched = {sid: v for sid, v in matches.items()}
test["lib"] = test.spectrum_id.map(lambda s: matched.get(s, (None, None))[0])

profile = {
    "train": {
        "spectra": pf.metadata.num_rows,
        "structures": len(struct_spectra),
        "spectra_per_structure": q(list(struct_spectra.values())),
        "by_library": {lib: {"spectra": lib_spectra[lib], "structures": len(lib_structs[lib]),
                             "ppm_error": q(ppm_by_lib[lib]), "num_peaks": q(peaks_by_lib[lib])}
                       for lib, _ in lib_spectra.most_common()},
        "share_with_test_adducts": round(n_test_adduct / pf.metadata.num_rows, 4),
        "adducts_top20": dict(adducts.most_common(20)),
        "adducts_timstof": dict(adducts_timstof.most_common(12)),
        "ionization_mode": dict(modes),
        "ce_units": dict(units),
        "timstof_precursor_mz": q(prec_timstof),
    },
    "overlap": {
        "np_examples_structures": len(np_ex),
        "np_examples_also_in_other_libs": len(np_ex & other),
        "np_examples_also_in_non_timstof_libs": len(np_ex & non_timstof),
        "enveda180_structures_also_in_non_timstof_libs": len(e180 & non_timstof),
        "enveda180_vs_np_examples": len(e180 & np_ex),
    },
    "test_public": test_profile,
    "test_public_located_in_train": {
        "matched_spectra": len(matched),
        "of": len(test),
        "by_library": test.lib.value_counts(dropna=False).to_dict(),
        "molecules_fully_matched": int(test.groupby("molecule_id").lib.apply(lambda s: s.notna().all()).sum()),
    },
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(profile, indent=2, default=str))
print(f"\nwrote {OUT}")
