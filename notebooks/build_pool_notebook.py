"""Build the canonical explanatory notebook under kaggle/publications/pool/."""
import sys
from pathlib import Path

import nbformat
from research_style import finalize
from nbformat.v4 import new_code_cell as C, new_markdown_cell as M, new_notebook

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "kaggle/publications/pool/notebook.ipynb"
cells = []
A = cells.append

A(M(r'''# CASMI26: A Bigger Candidate Pool Lowered My Score

If the answer is not in your candidate pool, your ranker scores 0 on that molecule. So the obvious next step is to **add more structure databases**. I tried three, with the same pipeline and the same ranker each time. The only thing that changed was the pool:

| Pool (on top of train + COCONUT) | New structures | Public LB |
|---|---|---|
| — | — | 0.258 |
| + PubChemLite for Exposomics | +361 k | **0.267** (ranker retrained on it: **0.273**) |
| + PubChemLite + ChEBI + LIPID MAPS | +415 k | 0.268 (−0.005) |
| + PubChemLite + top-250k PubChem compounds by PubMed count | +567 k | **0.252 (−0.021)** |

Two of the three additions *hurt*. My validation set got all three wrong: it said PubChemLite would **hurt**, and that the other two would change nothing.

This notebook diagnoses candidate coverage and decoy density from the competition data on CPU. The leaderboard tables are quoted historical results. It then shows the validation set that *did* track these leaderboard moves, and how to build one.

**TL;DR**
1. MRR is a ranking metric. Every structure you add inside a query's ±10 ppm window is a new competitor for rank 1. A new source helps only if the **answers** it adds outweigh the **decoys**.
2. How well a source covers the training libraries is a poor guide. ChEBI covered those structures about as well as PubChemLite per structure added, yet it lost on the leaderboard.
3. A validation set of *famous* natural products (like `enveda-np-examples`) barely sees decoy costs. Its answers are easy to find, and nearly all are in COCONUT, so it cannot see a new source's benefit either. Harder held-out queries tracked both leaderboard losses.'''))

A(M("## 1. Setup — competition data plus four public structure lists"))
A(C(r'''import glob, re, time, os
from pathlib import Path
import numpy as np, pandas as pd, pyarrow.parquet as pq
import matplotlib.pyplot as plt

T0 = time.time()
def find(pattern):
    hits = sorted(glob.glob(f"/kaggle/input/**/{pattern}", recursive=True)) or \
           sorted(glob.glob(str(Path(os.environ.get("CASMI_DATA_DIR", "data")) / "**" / pattern), recursive=True))
    assert hits, pattern
    return Path(hits[0])

TRAIN = find("train.parquet")
MASS = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956, "P": 30.97376163,
        "S": 31.97207100, "F": 18.99840322, "Cl": 34.96885268, "Br": 78.9183371, "I": 126.904473,
        "Si": 27.9769265325, "B": 11.0093054, "Se": 79.9165213, "Na": 22.9897692809, "K": 38.96370668,
        "As": 74.9215965}

def formula_mass(f):
    m = 0.0
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f or ""):
        if el not in MASS:
            return np.nan
        m += MASS[el] * (int(n) if n else 1)
    return m if m > 0 else np.nan'''))

A(M(r'''**The queries.** A test molecule is a natural product measured on a timsTOF. The closest labelled stand-in is the set of structures in the natural-product spectral libraries of `train.parquet`. We use their masses and InChIKey first blocks, restricted to the test's mass range.'''))
A(C(r'''t = pq.read_table(TRAIN, columns=["ingest_lib", "inchikey14", "molecular_formula"]).to_pandas()
NP_LIBS = ["riken", "gnps", "massbank", "mona", "msdial", "spectraverse"]
train_keys = set(t.inchikey14)
q = t[t.ingest_lib.isin(NP_LIBS)].drop_duplicates("inchikey14")
q = q.assign(mass=q.molecular_formula.map(formula_mass)).dropna(subset=["mass"])
q = q[q.mass.between(150, 1170)]
print(f"{len(t):,} spectra | {len(train_keys):,} training structures | {len(q):,} NP-library query structures")'''))

A(M(r'''**The sources.** Each pool keeps single-component, uncharged structures of 100–1300 Da, keyed by InChIKey block 1. Every structure in `train` is always in the pool, because the libraries themselves are candidates.
- COCONUT (CC0): `coconut_csv_lite-09-2026`
- PubChemLite for Exposomics (CC BY 4.0, Schymanski et al., Zenodo 22143470)
- ChEBI (CC BY 4.0) and LIPID MAPS (CC0), from the dataset `casmi-extra-structures`
- The 250,000 PubChem compounds with the most PubMed-linked articles (PubChem `CID-PMID` + `CID-SMILES`), from the dataset `casmi-pubchem-literature-250k`'''))
A(C(r'''def neutral_single(smiles):
    s = smiles.fillna("")
    charged = s.str.count(r"\+\]") != s.str.count(r"-\]")          # net charge != 0 (zwitterions pass)
    return ~s.str.contains(".", regex=False) & ~charged

coco = pd.read_csv(find("coconut_csv_lite-09-2026.*"),
                   usecols=["canonical_smiles", "standard_inchi_key", "exact_molecular_weight", "formal_charge"])
coco = coco[(coco.formal_charge == 0) & ~coco.canonical_smiles.str.contains(".", regex=False)]
coco = pd.DataFrame({"key": coco.standard_inchi_key.str[:14], "mass": coco.exact_molecular_weight})

pcl = pd.read_csv(find("PubChemLite_CCSbase_*.csv"), usecols=["FirstBlock", "SMILES", "MonoisotopicMass"])
pcl = pcl[neutral_single(pcl.SMILES)]
pcl = pd.DataFrame({"key": pcl.FirstBlock, "mass": pcl.MonoisotopicMass})

chebi = pd.read_parquet(find("extra_chebi.parquet"), columns=["pool_key", "mass"]).rename(columns={"pool_key": "key"})
lipidmaps = pd.read_parquet(find("extra_lipidmaps.parquet"), columns=["pool_key", "mass"]).rename(columns={"pool_key": "key"})
publit = pd.read_parquet(find("extra_pubchem.parquet"), columns=["pool_key", "mass"]).rename(columns={"pool_key": "key"})

train_struct = t.drop_duplicates("inchikey14")
train_struct = pd.DataFrame({"key": train_struct.inchikey14, "mass": train_struct.molecular_formula.map(formula_mass)})
SOURCES = {"COCONUT": coco, "PubChemLite": pcl, "ChEBI": chebi, "LIPID MAPS": lipidmaps, "PubChem-lit 250k": publit}
for name, d in SOURCES.items():
    d.dropna(inplace=True); d.query("100 <= mass <= 1300", inplace=True); d.drop_duplicates("key", inplace=True)
    print(f"{name:<17} {len(d):>9,} structures")'''))

A(M(r'''## 2. Answers versus decoys, per pool

For each pool we measure two things on the NP-library queries:
- **coverage**: the share of query structures that the pool contains through a non-library source. In the test, a Class-2 answer by definition has no spectra, so it can only be found this way.
- **window size**: how many pool structures fall within ±10 ppm of the query mass. Each one competes with the answer for rank 1.'''))
A(C(r'''POOLS = {
    "train + COCONUT":            ["COCONUT"],
    "+ PubChemLite":              ["COCONUT", "PubChemLite"],
    "+ PubChemLite, ChEBI, LM":   ["COCONUT", "PubChemLite", "ChEBI", "LIPID MAPS"],
    "+ PubChemLite, PubChem-lit": ["COCONUT", "PubChemLite", "PubChem-lit 250k"],
}
LB = {"train + COCONUT": 0.258, "+ PubChemLite": 0.273, "+ PubChemLite, ChEBI, LM": 0.268, "+ PubChemLite, PubChem-lit": 0.252}

qmass, qkeys = q.mass.to_numpy(), set(q.inchikey14)
rows = []
for name, srcs in POOLS.items():
    ext = pd.concat([SOURCES[s] for s in srcs]).drop_duplicates("key")
    pool = pd.concat([train_struct, ext]).drop_duplicates("key")
    m = np.sort(pool.mass.to_numpy())
    tol = qmass * 10e-6
    window = np.searchsorted(m, qmass + tol, "right") - np.searchsorted(m, qmass - tol, "left")
    rows.append({"pool": name, "structures": len(pool), "coverage": len(qkeys & set(ext.key)) / len(qkeys),
                 "median window": np.median(window), "mean window": window.mean(), "public LB": LB[name]})
res = pd.DataFrame(rows)
base = res.iloc[0]
res["Δ coverage (pts)"] = 100 * (res.coverage - base.coverage)
res["Δ window (%)"] = 100 * (res["mean window"] / base["mean window"] - 1)
print(res.round(3).to_string(index=False))'''))

A(M(r'''PubChemLite adds a lot of coverage. The other two sources add much less, while the windows keep growing. The PubChem literature subset grows the windows most of all and adds very few answers.

In the test, Class-1 answers are already in the pool (their spectra are in the libraries). So a new source can only help through Class-2 answers, and the decoys it adds are paid for by every molecule. Here is the trade-off per source, measured against the pool that already holds PubChemLite:'''))
A(C(r'''ref = pd.concat([train_struct, coco, pcl]).drop_duplicates("key")
ref_ext = set(pd.concat([coco, pcl]).key)
m_ref = np.sort(ref.mass.to_numpy())
tol = qmass * 10e-6
w_ref = (np.searchsorted(m_ref, qmass + tol, "right") - np.searchsorted(m_ref, qmass - tol, "left")).mean()
# reference: PubChemLite measured against the base pool (train + COCONUT), from the first table
pcl_row = res.iloc[1]
rows = [{"source": "PubChemLite (vs base pool)", "new structures": int(pcl_row.structures - base.structures),
         "answers gained": int(round((pcl_row.coverage - base.coverage) * len(qkeys))),
         "coverage gain (pts)": pcl_row["Δ coverage (pts)"], "decoys per window (+%)": pcl_row["Δ window (%)"],
         "answers per 1% more decoys": (pcl_row.coverage - base.coverage) * len(qkeys) / pcl_row["Δ window (%)"]}]
for name in ["ChEBI", "LIPID MAPS", "PubChem-lit 250k"]:
    d = SOURCES[name]
    new = d[~d.key.isin(ref.key)]
    m_new = np.sort(new.mass.to_numpy())
    added = (np.searchsorted(m_new, qmass + tol, "right") - np.searchsorted(m_new, qmass - tol, "left")).mean()
    gained = len((qkeys & set(d.key)) - ref_ext)
    rows.append({"source": name, "new structures": len(new), "answers gained": gained,
                 "coverage gain (pts)": 100 * gained / len(qkeys),
                 "decoys per window (+%)": 100 * added / w_ref,
                 "answers per 1% more decoys": gained / max(100 * added / w_ref, 1e-9)})
print(pd.DataFrame(rows).round(2).to_string(index=False))'''))

A(M(r'''The last column is the one that matters, and it tells only half the story:
- **PubChem-lit** buys very few answers for a lot of decoys. The leaderboard punished it hardest (−0.021), as this predicts.
- **ChEBI and LIPID MAPS** look *better* than PubChemLite on this measure. Yet adding them lost 0.005 on the leaderboard.

So coverage of *training-library* structures is itself an optimistic proxy. Those structures exist because someone had a reference standard to measure, and that is exactly what curated databases like ChEBI and LIPID MAPS contain. The test's Class-2 molecules have never been measured. They are natural-product isolates, not lipid standards, so the answers these sources add on paper mostly never appear in the test, while their decoys do.

**Rule I now follow:** a pool change is a hypothesis about the hidden test. Measure its cost on a validation set that can see decoys, and let the leaderboard judge its benefit. Change one source at a time.'''))

A(M(r'''## 3. Why my validation set was blind, and the one that was not

My first validation set (**V-in**) was the natural choice: the 250 structures of `enveda-np-examples`, the only data measured on the test instrument. For each of them I rebuilt the Class-2 situation: I removed every spectrum of the structure from every library and kept it in the pool.

It rejected PubChemLite (−0.035 on Class 2, the only effect it can see is the cost) and saw nothing for the other two (+0.001, −0.009). Two reasons, both visible in the data:
- the hosts describe these as *deliberately common* natural products. They have many close relatives in the libraries, so the answer usually sits far ahead of any new isomer;
- 249 of the 250 answers are in COCONUT, so V-in can never see the *benefit* of a new source either.

The second set (**V-real**) uses 600 structures from the training libraries that my ranker never trained on (grouped by the historical structure-identity protocol; this is not a verified scaffold split). Half are COCONUT members; half are known only through PubChemLite. They are harder and less famous, which is closer to what the leaderboard implies about the real Class 2.

Numbers from my pipeline (logged runs; the ranker was retrained once, for the PubChemLite pool, as on the leaderboard; V-real = Class-2 MRR on its COCONUT half):'''))
A(C(r'''runs = pd.DataFrame({
    "pool": ["train + COCONUT", "+ PubChemLite", "+ PubChemLite, ChEBI, LM", "+ PubChemLite, PubChem-lit"],
    "V-in Class 2":   [0.662, 0.619, 0.620, 0.610],
    "V-real Class 2": [0.594, 0.548, 0.529, 0.518],
    "public LB":      [0.258, 0.273, 0.268, 0.252],
})
print(runs.to_string(index=False))

fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
x = np.arange(len(runs))
for col, mk in [("V-in Class 2", "o"), ("V-real Class 2", "s")]:
    ax[0].plot(x, runs[col] - runs[col].iloc[1], marker=mk, label=col)
ax[0].plot(x, 2 * (runs["public LB"] - runs["public LB"].iloc[1]), marker="^", color="#DDB26F", label="public LB (×2)")
ax[0].axhline(0, color="grey", lw=0.8)
ax[0].set_xticks(x, ["base", "+PCL", "+PCL+ChEBI+LM", "+PCL+PubChem"], rotation=15)
ax[0].set_title("change vs. the PubChemLite pool"); ax[0].legend(fontsize=8)
ax[1].bar(res.pool.str.replace("+ ", "+\n"), res["mean window"], color="#4EA8FF")
ax[1].set_title("mean candidates in a ±10 ppm window"); ax[1].tick_params(axis="x", labelsize=8)
plt.tight_layout(); plt.show()'''))

A(M(r'''V-real's Class-2 score falls with every source added after PubChemLite, and it falls further for the source the leaderboard punished more. It cannot judge PubChemLite itself on this half, because every answer there is in COCONUT: like V-in, it only sees the cost. Its *level* is informative too. For the three PubChemLite-based pools, half of V-real's Class-2 MRR lands within 0.01 of the public LB. For the base pool it overshoots (0.30 vs 0.26), consistent with answers that only PubChemLite contains. Half of V-in's level would predict about 0.31 everywhere.

## 4. Takeaways
- **Pool size is a trade, not a free win.** Check how many decoys a source adds inside your query windows, not just how many answers it covers.
- **Validate on queries that look like the test, not like a showcase.** A held-out set of famous compounds tells you how well you retrieve famous compounds.
- **Group metric-equivalent structures before splitting.** Tautomer-canonical keys prevent one identity leak; scaffold independence needs a separate scaffold definition and audit.
- **One change per submission** when the leaderboard is your only judge of a benefit.

*Previous notebook in this series: [Your Ranker Learns How You Built the Pool](https://www.kaggle.com/code/senanuretin/casmi26-your-ranker-learns-how-you-built-the-pool): shortcut features in knockout-trained rankers.*'''))
A(C(r'''print(f"notebook ran in {time.time() - T0:.0f}s")'''))

nb = new_notebook(cells=cells, metadata={"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                                         "language_info": {"name": "python"}})
nbformat.write(finalize(nb), OUT)
print(f"wrote {OUT}")
