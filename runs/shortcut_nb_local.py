import matplotlib; matplotlib.use('Agg')
import glob, re, time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import lightgbm as lgb
import matplotlib.pyplot as plt

SEED = 42
rng = np.random.default_rng(SEED)
T0 = time.time()
hits = sorted(glob.glob("/kaggle/input/**/train.parquet", recursive=True))
DATA = Path(hits[0]).parent if hits else Path("C:/Users/MONSTER/kaggle-data/casmi26")
TRAIN = pq.ParquetFile(DATA / "train.parquet")
print(DATA, TRAIN.metadata.num_rows, "spectra")

TEST_ADDUCTS = {"[M+H]+", "[M+NH4]+", "[M-H2O+H]+", "[M-2H2O+H]+", "[M+Na]+", "[M+K]+",
                "[M-H]-", "[M-H2O-H]-", "[M+CH2O2-H]-", "[M+Cl]-"}
MASS = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956, "P": 30.97376163,
        "S": 31.97207100, "F": 18.99840322, "Cl": 34.96885268, "Br": 78.9183371, "I": 126.904473,
        "Si": 27.9769265325, "B": 11.0093054, "Se": 79.9165213, "Na": 22.9897692809, "K": 38.96370668}
E, PROTON = 0.000548579909, 1.007276466621
# ion m/z = M + shift for the ten test adducts (electron mass included)
SHIFT = {"[M+H]+": PROTON, "[M+NH4]+": 18.033826, "[M-H2O+H]+": PROTON - 18.010565,
         "[M-2H2O+H]+": PROTON - 36.021129, "[M+Na]+": 22.989221, "[M+K]+": 38.963158,
         "[M-H]-": -PROTON, "[M-H2O-H]-": -PROTON - 18.010565, "[M+CH2O2-H]-": 44.998203, "[M+Cl]-": 34.969401}

def formula_mass(f):
    tot = 0.0
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if el not in MASS:
            return np.nan
        tot += MASS[el] * (int(n) if n else 1)
    return tot

cols = ["ingest_lib", "inchikey14", "molecular_formula", "adduct", "precursor_mz", "num_peaks"]
meta = pd.concat([TRAIN.read_row_group(i, columns=cols).to_pandas() for i in range(TRAIN.num_row_groups)], ignore_index=True)
meta["row"] = np.arange(len(meta))
meta["fmass"] = meta.molecular_formula.map(formula_mass)
meta["nmass"] = meta.precursor_mz - meta.adduct.map(SHIFT)
meta["ppm"] = (meta.nmass - meta.fmass) / meta.fmass * 1e6
ok = meta.adduct.isin(TEST_ADDUCTS) & (meta.num_peaks >= 5) & (meta.ppm.abs() <= 20)
lib = meta[ok]
print(f"usable reference spectra: {ok.sum():,} of {len(meta):,}")

pool = lib.groupby("inchikey14").fmass.first().sort_values()
pool_keys, pool_mass = pool.index.to_numpy(), pool.to_numpy()
rows_of = lib.groupby("inchikey14").row.apply(np.asarray).to_dict()
lib_of_row = meta.ingest_lib.to_numpy()
print(f"pool: {len(pool):,} structures")

def make_queries(libs, keys=None, n=None, max_spectra=5):
    q = lib[lib.ingest_lib.isin(libs)]
    if keys is not None:
        q = q[q.inchikey14.isin(keys)]
    best = q.groupby(["inchikey14", "ingest_lib"]).size().reset_index(name="k").sort_values("k", ascending=False)
    best = best.drop_duplicates("inchikey14")
    if n is not None and len(best) > n:
        best = best.sample(n, random_state=SEED)
    out = []
    for key, l in zip(best.inchikey14, best.ingest_lib):
        r = q.row[(q.inchikey14 == key) & (q.ingest_lib == l)].to_numpy()
        if len(r) > max_spectra:
            r = np.sort(rng.choice(r, max_spectra, replace=False))
        out.append((key, l, r))
    return out

vin_q = make_queries(["enveda-np-examples"])
vin_keys = {k for k, _, _ in vin_q}
np_libs = ["riken", "gnps", "massbank", "mona", "msdial", "spectraverse"]
cand = set(lib.inchikey14[lib.ingest_lib.isin(np_libs)]) - vin_keys
train_q = make_queries(np_libs, keys=cand, n=700)
print(f"V-in queries: {len(vin_q)} | training queries: {len(train_q)}")

def window(m, tol=10.0):
    d = m * tol * 1e-6
    return np.arange(np.searchsorted(pool_mass, m - d), np.searchsorted(pool_mass, m + d, side="right"))

def query_mass(rows):
    return float(np.median(meta.nmass.to_numpy()[rows]))

need = set()
for key, _, r in vin_q + train_q:
    need.update(r.tolist())
    for i in window(query_mass(r)):
        need.update(rows_of[pool_keys[i]][:12].tolist())
need = np.array(sorted(need))

rg_start = np.cumsum([0] + [TRAIN.metadata.row_group(i).num_rows for i in range(TRAIN.num_row_groups)])
PEAKS = {}
for g in range(TRAIN.num_row_groups):
    sel = need[(need >= rg_start[g]) & (need < rg_start[g + 1])] - rg_start[g]
    if len(sel):
        t = TRAIN.read_row_group(g, columns=["ms2_mzs", "ms2_normalized_intensities"]).take(sel)
        for j, (mz, it) in enumerate(zip(t.column(0).to_pylist(), t.column(1).to_pylist())):
            b, inv = np.unique(np.round(np.asarray(mz) * 100).astype(np.int64), return_inverse=True)   # 0.01 Da bins
            v = np.bincount(inv, weights=np.sqrt(np.asarray(it)))
            PEAKS[int(sel[j] + rg_start[g])] = (b, v / np.linalg.norm(v))
print(f"loaded {len(PEAKS):,} spectra in {time.time() - T0:.0f}s")

def cosine(a, b):
    ia, va = a; ib, vb = b
    common, xa, xb = np.intersect1d(ia, ib, assume_unique=True, return_indices=True)
    return float(va[xa] @ vb[xb])

HIDDEN_B = set(rng.choice(pool_keys, int(0.6 * len(pool_keys)), replace=False))

def group(key, qlib, qrows, cls, design):
    qspec = [PEAKS[r] for r in qrows]
    idx = window(query_mass(qrows))
    rows = []
    for i in idx:
        k = pool_keys[i]
        refs = [r for r in rows_of[k][:12] if lib_of_row[r] != qlib]          # never the query's own library
        if cls == "C2" and k == key:
            refs = []
        if design == "B" and k in HIDDEN_B and not (cls == "C1" and k == key):
            refs = []
        d1 = max((cosine(q, PEAKS[r]) for q in qspec for r in refs), default=0.0)
        rows.append((k, d1, len(refs), abs(query_mass(qrows) - pool_mass[i]) / pool_mass[i] * 1e6, int(k == key)))
    g = pd.DataFrame(rows, columns=["cand", "d1", "n_refs", "abs_ppm", "label"])
    if g.label.sum() == 0 or (cls == "C1" and g.loc[g.label == 1, "n_refs"].iat[0] == 0):
        return None                                   # not a valid query of this class
    g["d1_gap"] = g.d1.max() - g.d1
    g["d1_rank"] = g.d1.rank(ascending=False, method="min")
    g["log_n_cand"] = np.log1p(len(g))
    g["qid"], g["cls"], g["design"] = key, cls, design
    return g

tables = []
for name, qs in (("train", train_q), ("vin", vin_q)):
    for design in ("A", "B"):
        for cls in ("C1", "C2"):
            for key, qlib, qrows in qs:
                g = group(key, qlib, qrows, cls, design)
                if g is not None:
                    tables.append(g.assign(split=name))
T = pd.concat(tables, ignore_index=True)
print(T.groupby(["split", "design", "cls"]).qid.nunique().unstack(), f"\n{time.time() - T0:.0f}s")

EVIDENCE = ["d1", "d1_gap", "d1_rank", "abs_ppm", "log_n_cand"]
WITH_PRIOR = EVIDENCE + ["n_refs"]

def fit(df, feats):
    df = df.sort_values(["qid", "cls"], kind="stable")
    ds = lgb.Dataset(df[feats], df.label, group=df.groupby(["qid", "cls"], sort=False).size().to_numpy())
    params = dict(objective="lambdarank", learning_rate=0.05, num_leaves=15, min_data_in_leaf=30,
                  lambdarank_truncation_level=25, verbose=-1, seed=SEED)
    return lgb.train(params, ds, 300)

def mrr(df, model, feats):
    df = df.assign(s=model.predict(df[feats]))
    df["r"] = df.groupby(["qid", "cls"]).s.rank(ascending=False, method="first")
    r = df[df.label == 1].groupby(["cls", "qid"]).r.min()
    return (1 / r).where(r <= 25, 0).groupby(level="cls").mean()

rows = []
for train_design in ("A", "B"):
    for fname, feats in (("evidence only", EVIDENCE), ("+ n_refs", WITH_PRIOR)):
        m = fit(T[(T.split == "train") & (T.design == train_design)], feats)
        for eval_design in ("A", "B"):
            s = mrr(T[(T.split == "vin") & (T.design == eval_design)], m, feats)
            rows.append({"trained on": train_design, "features": fname, "evaluated on": eval_design,
                         "V-in C1": s.get("C1", np.nan), "V-in C2": s.get("C2", np.nan)})
res = pd.DataFrame(rows)
print(res.round(3).to_string(index=False))

def within_query_auc(df, feat):
    vals = []
    for _, g in df.groupby("qid"):
        pos, neg = g.loc[g.label == 1, feat].to_numpy(), g.loc[g.label == 0, feat].to_numpy()
        if len(pos) and len(neg):
            vals.append((pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean())
    return float(np.mean(vals))

def knockout_audit(table, feats, cls="C2"):
    rows = []
    for f in feats:
        a = within_query_auc(table[(table.design == "A") & (table.cls == cls)], f)
        b = within_query_auc(table[(table.design == "B") & (table.cls == cls)], f)
        rows.append({"feature": f, "AUC, answer vs decoys (A)": a, "same, placebo knockout (B)": b,
                     "separation lost": (abs(a - 0.5) - abs(b - 0.5)) / max(abs(a - 0.5), 1e-9)})
    return pd.DataFrame(rows).sort_values("separation lost", ascending=False)

print(knockout_audit(T[T.split == "train"], WITH_PRIOR).round(3).to_string(index=False))

print(f"total runtime {time.time() - T0:.0f}s")