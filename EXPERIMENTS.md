# Experiment log

Every decision below is backed by a logged run in `runs/` and judged with `scripts/compare_runs.py`:
paired bootstrap over the same 250 V-in molecules (enveda-np-examples, the only instrument- and
pipeline-matched labelled data), mixture score over the class-share grid
f₁ ∈ {0.15 … 0.30}, f₃ ∈ {0.2, 0.3, 0.4}.
**ACCEPT** = the 95% CI lower bound of the expected-LB difference is > 0 at every grid point.

## Results

| # | Change | C1 MRR | C2 MRR | Expected-LB Δ (CI lower bounds) | Verdict |
|---|---|---|---|---|---|
| 0 | D1 direct library match, tie-break by mass error | 0.910 | 0.170 | — | baseline |
| 1 | + D2 mass-shifted analog propagation (interim max rule) | 0.908 | 0.252 | +0.015 … +0.034 | **ACCEPT** |
| 2 | LambdaRank ranker over D1 + D2 + priors + list context | 0.667 | 0.663 | — | superseded by 3 |
| 3 | ranker + D3 spectrum→fingerprint (out-of-fold) | 0.680 | 0.676 | +0.038 … +0.209 vs #1 | **ACCEPT** |
| 3b | D3 alone (#2 → #3) | +0.014 | +0.013 | −0.005 … −0.001 | NOT PROVEN |
| 4 | v2: Class-1 training groups with cross-library knockout | 0.590 | 0.719 | −0.025 … +0.001 vs #3 | NOT PROVEN |
| 4b | v2, Class-2 weight ×2 | 0.516 | 0.724 | < 0 vs #4 | NOT PROVEN (worse point) |
| 5 | v3: Class-1 groups from enveda-180 (timsTOF, cross-instrument) | 0.375 | 0.735 | −0.074 … −0.037 vs #4 | **REJECT** |
| 6 | #3 data without membership priors | 0.718 | 0.682 | — | diagnostic |
| 7 | **v13-clean**: C1 from both sources, no membership priors (incl. NaN-coded NP-likeness) | **0.839** | **0.662** | +0.085 … +0.225 vs #1; −0.006 … +0.028 vs #3 | **ACCEPT vs #1; chosen over #3 on robustness** |
| 9 | #7 on a pool + PubChemLite (1.08 M structures, median 96 candidates vs 58) | 0.832 | 0.627 | −0.038 … −0.022 | **REJECT on V-in** — but V-in cannot see the benefit (all its answers are in COCONUT); LB probe decides |
| 8 | fusion / gate of #7 with the interim rule, tuned on cross-fitted training OOF | 0.839 | 0.662 | ±0 (tuning chose "ranker only") | no change |

## Findings that changed the design

1. **The public `test.parquet` cannot validate anything.** All 400 molecules / 1,213 spectra are
   enveda-180 spectra copied verbatim from train: drug-like, 100% Class 1, no water-loss adducts.
   Our end-to-end run scores 0.9975 on it — a plumbing check, not a model score.
2. **`enveda-np-examples` leaks when held out by library name.** All 250 structures also occur in
   other libraries. Knockouts must be by structure identity.
3. **Structure identity must be the metric's key.** The standard InChIKey14 equals the metric's
   tautomer-canonical key for only 88.5% of COCONUT structures; 33 of 250 V-in molecules (13%) also
   sit in the pool under a *different* standard key. Knocking out only the query's own key would
   leak those into the Class-2 simulation. A heavy-atom *skeleton key* (exact necessary condition,
   ~0.2 ms) makes the full equivalence check affordable.
4. **`precursor_error_ppm` follows different conventions per library.** It equals our signed ppm
   exactly on enveda-np-examples (100%), equals |ppm| computed *without the electron mass* on gnps,
   masaryk and drug_plus (94–99%), and matches neither on enveda-180 or pluskal. We recompute it.
5. **Class-1 simulation must be cross-library (run #2/#3 diagnosis).** Removing only the query's
   own spectra left ≥ 25% of training "Class-1" groups with no reference at all (Class 2 in
   disguise) and let the rest match near-duplicate spectra from the same library and protocol
   (truth similarity median 0.912 vs 0.848 on V-in). The ranker learned to distrust direct
   matches and V-in C1 fell from 0.91 to 0.68. v2 removes the structure's whole query library and
   drops groups with no cross-library reference.

6. **Membership priors are shortcuts (runs #4–#7).** Features such as "in the training library",
   "in COCONUT", "only in enveda-180" or a spectrum count describe how the candidate pool was
   *assembled*, and the ranker used them to recognise the simulated class instead of weighing
   evidence. v3 made it visible: its Class-1 truths were drug-like enveda-180 structures (not in
   COCONUT), so no training truth was ever "in train and in COCONUT" — exactly what V-in's Class-1
   truths are — and V-in C1 fell to 0.375 although the similarity distribution now matched V-in
   (truth median 0.832 vs 0.848). NP-likeness leaked the same information through NaN (it is only
   defined for COCONUT rows). Dropping all of them lifted C1 to 0.839 with C2 at 0.662.
   **Decision on #7 vs #3:** not a proven improvement under the rule (the lower CI bound touches
   zero at the Class-2-heavy end), but #3 relies on features whose class-conditional distribution
   is set by our simulation design and cannot be checked against the test. #7 is chosen as a
   robustness decision, recorded as such.

7. **V-in's Class 1 is unusually easy for the simple rule (run #8).** On cross-fitted training
   queries the ranker beats the direct-match rule on Class 1 (0.873 vs 0.790), while on V-in the
   rule wins (0.908 vs 0.839). The hosts describe enveda-np-examples as *deliberately common*
   natural products with references in many libraries. Tuning a fusion weight or gate on V-in
   would have bought V-in-specific score; tuned on training OOF, every setting lost to the ranker
   alone, so none was adopted. Real Class-1 test molecules are likelier to resemble the training
   queries than V-in's showcase compounds.

## Leaderboard

| Date | Submission | Public LB | Note |
|---|---|---|---|
| 2026-09-24 | #7 v13-clean (D1+D2+D3, no membership priors) | **0.258** | two identical submissions scored identically — the pipeline is deterministic on the hidden test |
| 2026-09-24 | #1 interim rule max(D1, D2) — calibration | 0.190 | the learned ranker transfers: +0.068 |
| 2026-09-24 | #3 v1 ranker *with* membership priors — calibration | 0.245 | dropping the priors helped on the real test: +0.013 |

**Calibration.** Validation predicted 0.42–0.50 over the plausible class mixes; the leaderboard gave
0.258. With Class 1 near its V-in level (~0.84) at ~16%, Class 1 contributes ~0.14, leaving ~0.12
for Class 2 — i.e. real Class-2 MRR ≈ 0.25–0.30 against 0.66 in simulation. Leading hypotheses:
(1) V-in's Class 2 is too easy — deliberately common natural products with many close analogs in
the libraries, whereas real Class-2 molecules have never had a spectrum published; (2) dropping all
membership priors also dropped a legitimate one — enveda-180-only structures are drug-like and
cannot be answers for a natural-product test set.

**After the calibration submissions.** The leaderboard *orders* the three models exactly as
validation does (rule < v1 < v13-clean), so validation is a reliable guide for choices; only its
level is optimistic. Solving LB = a + w₁·MRR₁ + w₂·MRR₂ over the three submissions gives
w₁ ≈ 0.10, w₂ ≈ 0.18, a ≈ 0.055 — a calibration (3 points, 3 unknowns), not yet a test. It implies
a V-in Class-2 gain is worth ~1.8× the same Class-1 gain, and that even a perfect simulated Class 2
would add only ~0.06: most of the gap to the top (0.417) lies where V-in cannot see — pool
coverage for PubChem-only answers and Class 3. Hypothesis (2) above was refuted: the priors did not
carry legitimate signal on the real test. The nearest spectral neighbour of a Class-2 answer is
closer for V-in (median Tanimoto 0.783) than for training queries (0.731), a modest part of V-in's
optimism.

## Engineering notes

- Spectrum store: 2.54 M spectra → 0.87 GB memmaps, built in 138 s after replacing the NumPy
  cleaner with a numba version (verified identical on 900 random spectra; ~65× faster).
- Two memory-heavy processes (rank-data builder ~3.5 GB, GPU training ~8.7 GB commit incl. the CUDA
  driver) exceeded the 16 GB laptop's commit headroom; stages now run strictly in sequence.
- D3: 3.9 M-parameter set transformer, 2.0 M training spectra, ~2,800 spectra/s on the RTX 3060
  laptop GPU, 8 min/epoch. Held-out-fold BCE 0.0672 (prior-only ≈ 0.109); the true fingerprint
  out-scores a random decoy 83% of the time. Trained locally: no Kaggle GPU quota used.
- End-to-end submission path on the public file: 643 s for 400 molecules, 0 fallbacks.
