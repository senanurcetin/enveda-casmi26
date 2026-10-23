"""Generates notebooks/casmi26-ranker-shortcut.ipynb (public Kaggle notebook #1)."""
import sys
from pathlib import Path

import nbformat
from research_style import finalize
from nbformat.v4 import new_code_cell as C, new_markdown_cell as M, new_notebook

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("casmi26-ranker-shortcut.ipynb")
cells = []
A = cells.append

A(M(r'''# CASMI26: Your Ranker Learns How You Built the Pool

Most strong pipelines here share a shape: pull candidates from a structure pool by mass, describe each candidate with evidence channels, and train a learning-to-rank model on **simulated queries** — a molecule's spectra are knocked out of the library to fake a Class-2 molecule, left in to fake a Class-1 one.

This notebook is about a failure mode of that setup that I hit in my own pipeline and that no validation number warned me about:

> **Any feature whose value for the answer is *set by the knockout* — a direct-match score that is zero because you removed the spectra, a spectrum count, a library-membership flag — can tell the ranker which simulated class it is looking at. The ranker then scores your simulation, not the chemistry.**

How much it can exploit depends on *your simulation design*, not on the test set. So your CV can be excellent and wrong at the same time.

**TL;DR — measured below on the competition data (no external data, ~5 min on CPU):**

| Ranker trained on | Evaluated on | V-in Class 1 | V-in Class 2 |
|---|---|---|---|
| naive knockout (A) | naive knockout (A) — *what your CV shows* | 0.75 | **0.94** |
| naive knockout (A) | realistic pool (B) — *what the test does* | **0.33** | **0.48** |
| realistic pool (B) | realistic pool (B) | 0.64 | 0.44 |
| realistic pool (B) **+ a spectrum-count feature** | realistic pool (B) | **0.48** | 0.44 |

Plus a **knockout audit** you can run on your own ranker table, and the numbers from my full pipeline, where removing pool-membership features moved held-out Class-1 MRR from **0.375 to 0.839**.'''))

A(M("## 1. Setup"))
A(C(r'''import glob, re, time, os
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
DATA = Path(hits[0]).parent if hits else Path(os.environ.get("CASMI_DATA_DIR", "data"))
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
    return tot'''))

A(M(r'''## 2. A deliberately small retrieval engine

Everything uses the training data only:

- **pool** — every training structure (`inchikey14`), keyed by the neutral monoisotopic mass of its formula;
- **candidates** — pool structures within ±10 ppm of the query's neutral mass;
- **evidence** — `d1`: the best cosine between the query's spectra and the candidate's own reference spectra (0.01 Da bins, √intensity);
- **a pool-membership feature** — `n_refs`: how many reference spectra the candidate has in the (knocked-out) library. It looks like harmless metadata. It is the subject of this notebook.

Queries are molecule-like: one structure, up to 5 of its spectra from a single library, test adducts only. **V-in** = the 250 `enveda-np-examples` structures (the only instrument- and pipeline-matched natural products). **Training queries** = natural-product structures from the other public libraries, with any V-in structure excluded.'''))

A(C(r'''cols = ["ingest_lib", "inchikey14", "molecular_formula", "adduct", "precursor_mz", "num_peaks"]
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
print(f"pool: {len(pool):,} structures")'''))

A(C(r'''def make_queries(libs, keys=None, n=None, max_spectra=5):
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
print(f"V-in queries: {len(vin_q)} | training queries: {len(train_q)}")'''))

A(M(r'''Only the peaks we actually need are read: the query spectra and the reference spectra of every candidate within the mass windows.'''))
A(C(r'''def window(m, tol=10.0):
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
    return float(va[xa] @ vb[xb])'''))

A(M(r'''## 3. Two ways to simulate a Class-2 molecule

**Design A — the natural first attempt.** Hide the true structure's spectra; leave everything else. A simulated Class-2 query then has exactly one candidate with no reference spectra: the answer.

**Design B — closer to the test.** In the real candidate pool most structures have *no* spectra at all (COCONUT / PubChem entries; in my pool 440k of 715k). So additionally hide the spectra of a fixed random 60% of all structures. The answer is now one spectrum-less candidate among many.

Class 1 is simulated in both designs by removing only the query library's spectra of the structure (so matches are cross-library, as in the test).'''))
A(C(r'''HIDDEN_B = set(rng.choice(pool_keys, int(0.6 * len(pool_keys)), replace=False))

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
print(T.groupby(["split", "design", "cls"]).qid.nunique().unstack(), f"\n{time.time() - T0:.0f}s")'''))

A(M("## 4. Train on one design, evaluate on the other"))
A(C(r'''EVIDENCE = ["d1", "d1_gap", "d1_rank", "abs_ppm", "log_n_cand"]
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
res.style.format({"V-in C1": "{:.3f}", "V-in C2": "{:.3f}"}).background_gradient(subset=["V-in C1", "V-in C2"], cmap="RdYlGn", vmin=0, vmax=1)'''))

A(C(r'''fig, ax = plt.subplots(figsize=(10, 3.8))
lab = res["trained on"] + "→" + res["evaluated on"] + " | " + res["features"]
x = np.arange(len(res))
ax.bar(x - 0.2, res["V-in C1"], 0.4, label="Class 1", color="#00E0FF")
ax.bar(x + 0.2, res["V-in C2"], 0.4, label="Class 2", color="#8F9BFF")
ax.set_xticks(x, lab, fontsize=8, rotation=20, ha="right")
ax.set_ylabel("V-in MRR@25"); ax.set_ylim(0, 1.05); ax.legend(loc="upper right")
ax.set_title("Same queries, same model family — only the simulation design and one feature change")
plt.tight_layout()'''))
A(M(r'''Read the table in pairs.

**1. The naive design flatters everything — including honest evidence.** Trained and evaluated on design A, even the *evidence-only* ranker scores 0.94 on Class 2. Nothing chemical is being recognised: when every decoy has reference spectra and only the answer has none, "direct-match similarity = 0" *is* the answer. The same model evaluated on the realistic pool (B) drops to 0.48.

**2. It also breaks Class 1, in the opposite direction.** A ranker that learned "no direct match ⇒ answer" sees hundreds of such candidates in the realistic pool and ranks them above a Class-1 answer that *does* match: 0.75 → 0.33.

**3. A membership feature costs Class 1 even with a realistic design.** Trained and evaluated on B, adding `n_refs` leaves Class 2 unchanged (0.442 → 0.442) and takes Class 1 from **0.642 to 0.478**. The feature separates answers from decoys in *opposite directions* for the two classes, so the ranker has to guess the class first — and a spectrum count is a guess about the simulation, not the molecule.

Nothing here is specific to `n_refs`: library flags, source tags and spectrum counts all behave this way, and so does any property populated for only one source (NaN encodes the source).'''))

A(M(r'''## 5. An audit you can run on your own ranker table: is the feature detecting the knockout?

The question to ask of every feature is not "does it separate answers from decoys?" — in Class 2 even honest evidence does, because an unseen molecule has no direct match by definition. The question is **whether it separates them because of the knockout itself**.

So score the same queries twice: under your knockout (design A), and under a *placebo* knockout that also strips a large share of decoys the way the real pool does (design B). For each feature, compute the within-query AUC of "answer vs decoys" for Class 2 in both. A feature whose AUC collapses towards 0.5 when decoys are knocked out too was detecting "has been knocked out", not the molecule.'''))
A(C(r'''def within_query_auc(df, feat):
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

audit = knockout_audit(T[T.split == "train"], WITH_PRIOR).round(3)
print(audit.to_string(index=False))
audit'''))
A(M(r'''Two groups separate cleanly:

- **Knockout detectors** — `d1`, its gap and rank, and `n_refs` — lose about 60% of their separation once decoys are knocked out too. Under the naive design they were not measuring spectral agreement; they were measuring *whether a candidate had been knocked out*.
- **Genuinely structural features** — here the mass error — are untouched by the placebo.

Note that `d1` is a perfectly good feature; it is the naive *design* that turned it into a detector. The fix is therefore not to delete evidence but to **simulate the pool as well as the answer**, and then to drop features that remain class-identifying under the realistic design (point 3 above).'''))

A(M(r'''## 6. What this cost in my full pipeline

My pipeline (not reproduced here) ranks a 715k-structure pool (train ∪ COCONUT) with direct-match, mass-shifted analog and spectrum→fingerprint channels, plus prior features: in-train, in-COCONUT, only-in-enveda-180, spectrum count, and COCONUT's NP-likeness. Held-out MRR@25 on the 250 V-in molecules, structure-level knockout:

| Ranker training data | Priors | V-in C1 | V-in C2 |
|---|---|---|---|
| NP libraries, Class 1 = query's own spectra removed | all | 0.680 | 0.676 |
| Class 1 from `enveda-180` (timsTOF, cross-instrument) | all | **0.375** | 0.735 |
| both Class-1 sources | none of the membership priors | **0.839** | 0.662 |

The middle row is the instructive one. Its Class-1 answers were drug-like `enveda-180` structures, which are not in COCONUT — so no training answer was ever "in the library **and** in COCONUT", which is exactly what a natural-product Class-1 answer looks like. The ranker learned that combination as a decoy signature. The similarity distribution of those Class-1 answers matched V-in almost exactly (median 0.832 vs 0.848); the evidence was fine, the membership features were not. NP-likeness leaked the same information a second way: it is only defined for COCONUT rows, so NaN meant "not from COCONUT".

### Update — does it hold on the leaderboard?

I submitted the full pipeline with and without the membership priors (same code otherwise; the pipeline is deterministic — two identical submissions scored identically):

| Full pipeline | V-in C1 | V-in C2 | Public LB |
|---|---|---|---|
| simple rule, no learned ranker | 0.908 | 0.252 | 0.190 |
| ranker **with** membership priors | 0.680 | 0.676 | 0.245 |
| ranker **without** membership priors | 0.839 | 0.662 | **0.258** |

The leaderboard orders the three exactly as the knockout validation does, and removing the priors helped on the hidden test too. One comparison on the public split is not proof, and the absolute validation numbers are clearly optimistic (V-in's natural products are deliberately common, with close analogs in the libraries) — but the direction held.

## Takeaways

1. **Simulate the pool, not just the answer.** If most real candidates have no reference spectra, most simulated candidates must have none either — otherwise "no direct match" becomes the answer key, for evidence features as much as for metadata.
2. **Audit with a placebo knockout.** Score your table under your knockout and under one that also strips decoys; a feature whose separation collapses was detecting the knockout.
3. **Drop features that remain class-identifying.** Spectrum counts, library/source flags, and anything that is only defined for one source (NaN is a value) — they separate answers from decoys in opposite directions per class.
4. **Report per class.** A mixture score can hide a collapse in one class behind a gain in the other; in my pipeline Class 2 improved while Class 1 fell by half.
5. **Validate on queries whose answers look like the test's answers** — natural products, cross-instrument Class-1 matches.

*Questions and counter-examples are very welcome in the comments — especially if your pipeline shows something different.*'''))
A(C(r'''print(f"total runtime {time.time() - T0:.0f}s")'''))

nb = new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                         "language_info": {"name": "python", "version": "3.12"}})
nbformat.write(finalize(nb), OUT)
print("wrote", OUT, len(cells), "cells")
