# Enveda CASMI 2026 — Solution Architecture

*Status: design v1 · 2026-09-23 · every number below is measured (`reports/data_profile.json`) or quoted from the competition pages, and labelled as such.*

---

## 1. The problem, precisely

| | |
|---|---|
| Input | 1–16 MS/MS spectra per molecule (median 3), Bruker timsTOF, one of 10 adducts |
| Output | up to 25 SMILES per `molecule_id`, best first |
| Metric | MRR@25 on the **tautomer-canonical InChIKey14** (RDKit **2026.03.3**, pinned) |
| Hidden test | ~400 molecules / ~1,500 spectra, natural-product-like, mass 157–1,159 Da (median 348) |
| Runtime | Code competition: ≤ 9 h CPU or GPU, **internet off**, public external data and models allowed |
| Budget | 5 submissions/day, 2 final selections, deadline **2026-12-14** (merger/entry 2026-12-07) |

Three hidden novelty classes decide what is even reachable:

| Class | Definition | What can find it |
|---|---|---|
| 1 | has public reference spectra | spectral library search |
| 2 | in PubChem / COCONUT, no reference spectra | database retrieval + ranking |
| 3 | not in PubChem | de novo generation only |

### Metric details that constrain the design (from the official scoring notebook)

1. Rank = **position in the list**. An unparseable SMILES still consumes its slot, and so does a repeated structure → every output must be parsed and deduplicated on InChIKey14 **with RDKit 2026.03.3** before writing.
2. Only the first correct hit counts, so after rank ~10 each slot is worth ≤ 0.1 — the tail is for hedging, the head is where score lives.
3. A missing `molecule_id` scores 0 rather than failing, so a crash for one molecule must degrade to a fallback list, never to a missing row.

---

## 2. What the data actually looks like (measured)

**Train:** 2,539,608 spectra · 275,810 structures · median 6 spectra per structure.

| Library | Spectra | Structures | ppm error p50 / p95 | Peaks p50 / p95 | Note |
|---|---|---|---|---|---|
| enveda-180 | 1,153,785 | 182,941 | 0.99 / 3.35 | 136 / 830 | timsTOF, **drug-like, not NP** |
| pluskal_ms2 | 527,581 | 46,821 | 1.66 / **2,432** | 21 / 118 | Orbitrap, NCE |
| riken | 347,171 | 15,892 | 3.67 / 296 | 9 / 113 | plant metabolites |
| gnps | 220,849 | 45,750 | 1.70 / **1,877** | 46 / 668 | heterogeneous |
| massbank · mona · spectraverse · msdial · drug_plus · masaryk | 289,038 | — | small–heavy tails | 9–108 | mixed |
| **enveda-np-examples** | **1,184** | **250** | 1.54 / 3.64 | 142 / 847 | **the only in-domain labelled data** |

Findings that shape the architecture:

1. **The public `test.parquet` is useless for validation.** All 400 molecules / 1,213 spectra come from **enveda-180** and all 1,213 were located verbatim in train. It is drug-like, 100% Class 1 by construction, and lacks the water-loss adducts the real test contains. Use it only to test format and runtime.
2. **`enveda-np-examples` leaks when held out by library.** All 250 structures also exist in non-timsTOF libraries, and none in enveda-180. Holding it out by `ingest_lib` tests Class 1 only — this is why community CV on it does not track the leaderboard.
3. **Label quality varies by library.** timsTOF libraries sit under ~3.6 ppm at p95; pluskal and gnps have p95 errors in the thousands of ppm. Library spectra need a precursor-error filter before they are trusted as evidence.
4. **Peak counts are a domain shift.** timsTOF spectra carry ~140 peaks at the median and ~850 at p95 (public test sample median 230, p99 1,610); most other libraries carry 9–54. Peak cleaning is not a detail — it decides whether test spectra and library spectra are comparable.
5. **timsTOF training adducts are dominated by dimers** (`[2M+Na]+` 188k, `[2M+H]+` 89k) that never occur in test; 85.7% of all train spectra use one of the 10 test adducts.
6. **Train has no `spectrum_id`/`molecule_id`.** Structure identity (`inchikey14`) is the only grouping key, so every split must be a structure-level split.

### The class mix is not known — design for a range

- A library-only submission scores **0.151** on the public LB (reported by prvsiyan). Since Class 2/3 score 0 there, **f₁ ≥ 0.151**; ~16% if Class-1 MRR is 0.93, higher if cross-instrument Class-1 matching is weaker.
- The widely repeated "45% Class 2 / 39% Class 3" split is **not supported** by its own derivation (the posted arithmetic yields f₂ ≈ 18%). We treat f₂/f₃ as unknown.
- Public LB top: **0.417**; public-notebook ceiling ~0.33–0.34; identical code re-run moved 0.335 → 0.328 (seed noise ±0.007). With ~400 molecules split public/private, LB sampling noise is on the order of ±0.03 — larger than most single-change effects.

---

## 3. Design principles

1. **Retrieve, then rank.** An exact-InChIKey14 metric makes free generation a long shot; ranking real molecules is tractable. Generation is only a gated tail hedge for Class 3.
2. **Recall is a first-class metric.** Every stage reports pool recall ρ (is the truth in the list at all?) separately from ranking quality. A ranker cannot fix what retrieval dropped.
3. **Validation decides, the leaderboard confirms.** Constants are never tuned on the LB. LB submissions are spent on calibration probes and paired A/B tests of changes that already won locally.
4. **Structure-level isolation everywhere.** Held-out structures are removed from *every* library (and, for Class 3, from the pool) — never by library name.
5. **One config, logged runs.** Every run writes its config hash, local scores with confidence intervals, and LB score if submitted.

---

## 4. Pipeline

```
                    OFFLINE (local / Kaggle GPU)                           ONLINE (submission notebook, ≤ 9 h, no internet)
 ┌──────────────────────────────────────────────┐        ┌──────────────────────────────────────────────────────────────┐
 │ A1 spectrum store  (cleaned peaks, memmap)   │        │ B1 clean query spectra  → B2 neutral mass per molecule        │
 │ A2 structure table (RDKit 2026.03.3, FPs)    │  ───►  │ C  candidate generation (mass window over pool)  → recall ρ   │
 │ A3 candidate pool  (train ∪ COCONUT ∪ …)     │        │ D  evidence channels → per-candidate features                  │
 │ A4 models (spec→FP, ranker)                  │        │ E  LambdaRank ranker (seed-bagged)                             │
 └──────────────────────────────────────────────┘        │ F  tail hedge (de novo, gated) → G dedupe/validate → CSV       │
                                                          └──────────────────────────────────────────────────────────────┘
```

### A. Offline assets (built once, shipped as Kaggle datasets)

| Asset | Content | Key decisions |
|---|---|---|
| **A1 Spectrum store** | train spectra after cleaning, as flat `float32` peak arrays + offsets (memmap) and a light metadata table | cleaning identical to the query path; drop spectra with \|ppm\| > 20 or non-test adducts from the *evidence* library (keep them for model pretraining) |
| **A2 Structure table** | `inchikey14 → smiles, formula, monoisotopic mass, fingerprints` | computed with **RDKit 2026.03.3** so our keys equal the metric's keys |
| **A3 Candidate pool** | train structures ∪ **COCONUT 2.0** (CC0) ∪ ChEBI/LIPID MAPS (optional) ∪ **formula-filtered PubChem** (host confirmed prize-eligible) | stored as a mass-sorted array for O(log n) window lookup; each structure carries its **source tags** so the ranker can learn source priors |
| **A4 Models** | spectrum → fingerprint network; LambdaRank ranker | trained on structure-level splits only |

Two licensing constraints: competition data may not be redistributed to non-participants, so anything derived directly from `train.parquet` is rebuilt inside the notebook or kept in a **private** dataset; COCONUT/ChEBI/PubChem-derived assets can be public.

### B. Query processing

- **B1 Spectrum cleaning** (one function shared with A1): drop peaks > precursor + 2 Da, relative-intensity floor, window-based top-k (top ~6 per 50 Da keeps low-mass fragments that a global top-N discards), optional deisotoping, √-intensity. The settings are tuned on in-domain validation, not assumed.
- **B2 Neutral mass**: per spectrum from precursor m/z and the adduct table (10 test adducts); per molecule a robust consensus across its spectra, with a calibrated timsTOF offset (np-examples: median +1.5 ppm).
- **B3 Merged spectrum**: one consensus spectrum per molecule across collision energies and adducts, alongside the individual spectra. (The largest single effect measured publicly — +0.019 — came from combining per-spectrum and merged evidence.)

### C. Candidate generation

Mass window over the pool (start at ±10 ppm around the calibrated neutral mass, tuned on recall). Pool size and window are tuned **jointly on ρ and MRR**, because adding candidates raises recall but dilutes ranking — public experiments show naive PubChem expansion dropping Class-2 MRR from 0.52 to 0.35. PubChem therefore enters as a *separately tagged* source the ranker can learn to discount, not as an undifferentiated extension.

### D. Evidence channels (features per molecule × candidate)

| # | Channel | Signal | Main target |
|---|---|---|---|
| D1 | **Direct library match** | entropy similarity / cosine of query spectra vs the candidate's own library spectra; max & mean over spectra; instrument-match and adduct-match flags | Class 1 |
| D2 | **Analog propagation** | modified (mass-shifted) similarity to library spectra of *other* structures, propagated as `max sim(a)^p · Tanimoto(c, a)` | Class 2 |
| D3 | **Spectrum → fingerprint model** | network predicts P(bit) over ~7k fingerprint bits; candidate score = log-likelihood; per-spectrum and merged-spectrum heads averaged | Class 2 |
| D4 | **Fragment explanation** | share of query intensity explained by 1–2 bond cleavages / sub-formulas of the candidate | Class 2, isomer separation |
| D5 | **Priors** | source tags (COCONUT / enveda-180 / PubChem), NP-likeness, formula plausibility, mass error | all |
| D6 | **List context** | within-molecule ranks and z-scores of every channel, candidate count | ranking calibration |

D6 matters as much as any channel: a similarity of 0.4 means something different in a list of 5 candidates than in a list of 500.

### E. Ranker

LightGBM **LambdaRank** (objective aligned with reciprocal rank), grouped by molecule, trained on simulated queries (§5) with the **class mix drawn from the plausible range**, not a single guess. Seed-bagged for stability. Monotone constraints on D1/D2 similarity features to stop the ranker from learning spurious inversions on few examples.

### F. Class-3 hedge (research track, gated)

A generator conditioned on the predicted fingerprint and formula proposes novel structures for the **tail slots (21–25)** only when all channels are weak. It ships only if it improves the mixture score on the Class-3 simulation without costing Class 1/2 — each tail slot is worth ≤ 0.048 per hit, so the bar is "does not hurt".

### G. Output safety

Parse every SMILES and canonicalize with RDKit 2026.03.3 → drop failures → deduplicate on InChIKey14 → fill to 25 from the next-best candidates → assert the format rules from the metric → per-molecule `try/except` falling back to a mass-window list, so one bad molecule can never cost a row.

---

## 5. Validation framework — the core of this project

### Query sets

| Set | Queries | Why |
|---|---|---|
| **V-in** | enveda-np-examples: 250 structures, 1,184 spectra | only instrument- and pipeline-matched NP data |
| **V-ext** | NP-like structures from riken / gnps / massbank / mona (COCONUT members), ~2,000 held out | variance reduction; cross-instrument, so reported separately |
| Public `test.parquet` | 400 enveda-180 molecules | **format and runtime only** — it is Class 1 and drug-like |

### Class simulation by knockout (structure-level, all libraries)

| Simulated class | Library | Candidate pool |
|---|---|---|
| C1 | remove only the query spectra themselves | keep structure |
| C2 | remove **every** spectrum of the structure, in every library | keep structure |
| C3 | remove every spectrum | **remove structure** (retrieval scores 0; measures the hedge) |

### Scoring and decisions

- Report MRR per class, pool recall ρ, and the **mixture** `S(f) = f₁·MRR₁ + f₂·MRR₂ + f₃·MRR₃` over a grid of plausible `f` (f₁ ∈ [0.15, 0.30], remainder split across C2/C3).
- Compare variants **paired on the same queries**, with bootstrap confidence intervals on the difference. A change is accepted only if it wins across the whole `f` grid.
- Spend LB submissions on (a) calibration probes that pin down `f` (e.g. library-only, pool-only), and (b) confirming changes that already won locally. Never on tuning.
- Experiment log: `runs/<date>_<confighash>.json` with config, local scores + CIs, LB score.

This is also our publishing angle: a leak-free validation setup with measured class mixtures is exactly what the community threads are asking for.

---

## 6. Modeling plan and compute budget

### 6.1 Hardware, measured and documented

| | Local | Kaggle GPU session | Kaggle CPU session |
|---|---|---|---|
| Accelerator | RTX 3060 Laptop, 6 GB, ~10 TFLOPS FP32 | P100 16 GB (~9.3 TFLOPS) or 2× T4 | — |
| CPU / RAM | 8 cores / 16 GB (often ~3 GB free) | 4 cores / 29 GB | 4 cores / 30 GB |
| Limits | none | 12 h per session; weekly GPU quota shown in account settings | 12 h per session |
| Submission run | — | ≤ 9 h, internet off | ≤ 9 h, internet off |

The local GPU is in the same FP32 class as the P100, so **training runs locally by default** and costs no quota. The Kaggle quota is kept for what only Kaggle can do: parallel experiments, anything needing > 6 GB VRAM or > 16 GB RAM, and the submission itself. Planning assumes a conservative **15 GPU-hours/week** until the real quota is confirmed from settings.

### 6.2 The spectrum → fingerprint model (D3)

| Decision | Choice | Reason |
|---|---|---|
| Input | set of ≤ 64 cleaned peaks (m/z, √intensity, neutral-loss m/z) + precursor mass, adduct, CE, mode tokens | peaks are a set, not a sequence; neutral losses carry most of the structural signal |
| Encoder | 4-layer set transformer, d = 256, no positional encoding | permutation-invariant, ~3 M parameters, fits 6 GB with batch 512 in FP16 |
| Heads | 4,096-bit ECFP4 ‖ MACCS (sigmoid, focal loss) + formula-count auxiliary | the fingerprint is what candidates are scored against; the formula head regularises |
| Molecule level | per-spectrum and merged-spectrum predictions averaged | the one effect publicly measured above noise (+0.019) |
| Training data | ~2.1 M spectra with test adducts and \|ppm\| ≤ 20; timsTOF up-weighted | matches the test instrument without discarding cross-instrument chemistry |
| Leakage control | **structure-level 2-fold split** (A/B halves by metric key, V-in structures always held out) | ranker features for training queries must be out-of-fold, or the ranker learns to trust an overfit channel |

Cost per model version: three trainings (fold A, fold B, full) × ~20 epochs. At ~6 k spectra/s on either GPU class that is ~6 min/epoch → **~2 h per training, ~6 h per version**, i.e. 2–3 versions per week locally, or on Kaggle within a 15 h/week quota.

### 6.3 The ranker (E)

- **Queries:** Class-2/Class-1 knockout simulations on NP-like training structures (COCONUT members from non-enveda libraries, ~5 k), never the V-in structures.
- **Features:** all channels (D1–D6), computed with the query's knockout applied, D3 from the out-of-fold model.
- **Model:** LightGBM LambdaRank, group = molecule, early stopping on a structure-disjoint validation fold; monotone constraints on similarity features; 5 seeds bagged.
- **Selection:** by the mixture score on V-in across the `f` grid, with paired CIs against the current best.

### 6.4 Submission-run budget (≤ 9 h, target < 3 h)

| Stage | Work for ~400 molecules | Estimate |
|---|---|---|
| Load store, pool, models | ~3 GB from input datasets | 5–10 min |
| D1 direct match | ~60 candidates × ≤ 40 references × ~3 spectra | < 2 min |
| D2 analog search | sparse binned prefilter → top 500 → modified entropy | 10–20 min |
| D3 fingerprint model | ~1,500 spectra forward pass | < 1 min |
| D4 fragments, ranker, output | per-candidate CPU work | 10–20 min |

Timings are logged on every run of the public test file; any stage that grows past 3× its estimate blocks release.

### 6.5 Leaderboard budget

5 submissions/day × ~80 days ≈ 400 available; planned use ≤ 2/day. Each submission has a stated purpose in the run log: *calibration probe* (pins down the class mix), *confirmation* (a change that already won locally), or *final candidate*.

Local environment: project `.venv` created with `--system-site-packages`, so the heavy global stack (PyTorch + CUDA, LightGBM) is shared while `rdkit==2026.3.3` is pinned per project. Offline wheels are shipped to Kaggle as a dataset, mirroring the official metric notebook.

---

## 7. Repository layout

```
enveda-casmi26/
├── ARCHITECTURE.md          ← this document
├── configs/                 base.yaml + experiment overrides
├── src/casmi/
│   ├── config.py  io.py     data paths (env CASMI_DATA), parquet streaming
│   ├── chem/                RDKit utils, adduct table, masses, fingerprints, ik14 canonicalizer
│   ├── spectra/             cleaning, merging, similarity (entropy, cosine, modified cosine)
│   ├── pool/                candidate DB build + mass-window lookup
│   ├── channels/            library.py analog.py fpmodel.py fragments.py priors.py
│   ├── rank/                features.py (incl. list context), lambdarank.py
│   ├── validate/            knockout.py (C1/C2/C3), metric.py (official-equivalent), mixture.py
│   └── submit/              writer.py (parse → canonicalize → dedupe → fill → assert)
├── scripts/                 profile_data.py, build_store.py, build_pool.py, train_*.py, run_validation.py
├── notebooks/               kaggle_submission.ipynb (thin wrapper around src/)
├── reports/                 data_profile.json, validation reports
└── runs/                    experiment log (config hash → scores)
```

Data lives outside OneDrive (`C:/Users/MONSTER/kaggle-data/casmi26`) and outside git.

---

## 8. Milestones

| When | Deliverable | Exit criterion |
|---|---|---|
| Week 1 (→ 09-29) | env, official-equivalent metric, spectrum store, pool v1 (train ∪ COCONUT), **knockout validation harness**, D1 baseline | local metric matches official on the doctests; V-in C1 MRR reproduced |
| **Kaggle notebook #1** | "A leak-free validation setup for CASMI26" (+ the public-test-file finding) | published while the competition is young |
| Week 2–3 | D2 analog propagation, D5/D6, LambdaRank; first real submission | mixture score beats D1-only across the `f` grid |
| Week 4–6 | D3 spectrum→FP model, D4 fragments, pool expansion (ChEBI, tagged PubChem) | paired CI-positive gains |
| Week 7–10 | Class-3 hedge experiments, robustness, seed bagging | no regression on C1/C2 |
| Final 2 weeks | freeze; choose 2 finals: one best-mixture, one most-robust | runtime < 6 h with margin |

---

## 9. Risks

| Risk | Mitigation |
|---|---|
| Local validation does not track the LB (community's main complaint) | structure-level knockout, mixture over an `f` range, paired CIs; LB only for calibration |
| RDKit version drift changes keys | pin 2026.03.3 locally and on Kaggle; assert `rdkit.__version__` at startup |
| Runtime overrun on the hidden test | budget per stage, log timings on the public file, hard fallbacks |
| Overfitting to public LB noise (±0.03) | never tune on LB; final pick favours robustness |
| Licensing of derived assets | train-derived artifacts private or rebuilt in-notebook; COCONUT/ChEBI/PubChem public |
| Local RAM (16 GB, often ~3 GB free) | row-group streaming, memmaps, float32; heavy jobs on Kaggle |

## 10. As built (2026-09-30) — where the system departs from this plan

- **Pool:** train + COCONUT + PubChemLite. PubChem, ChEBI and LIPID MAPS were built and rejected on
  the leaderboard (they add decoys faster than answers).
- **D3** stayed at 2,048 Morgan bits and the small model; a 12 M-parameter version was better alone
  but not in the ranker.
- **D4** is implemented as fragment *explanation* (depth-2 bond cleavage, share of intensity/peaks);
  it was the largest single Class-2 gain. A learned bond-propensity likelihood (D5) was weaker.
- **Class-3 hedge (F)** is no longer a gated tail: generated candidates (42 biotransformations of the
  top-15 spectral analogs) enter the candidate list with full features and the ranker places them.
  Training includes Class-3 groups; up-weighting them hurt.
- **Priors (D5 in section 4)** were removed: membership features are simulation shortcuts.
- **Validation:** a second held-out set, V-real (600 frozen library queries), was added because
  V-in's showcase compounds could not see pool decoys. Decisions: V-real paired bootstrap, a
  significant V-in Class-1 drop as a veto, the leaderboard for anything V-real cannot see.
