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
| 10 | #9 with the ranker retrained on _pcl candidate lists | 0.842 | 0.619 | −0.019 … −0.008 vs #9 (point ≈ 0) | NOT PROVEN on V-in; +0.006 on LB |
| 11 | #10 on pool + ChEBI + LIPID MAPS (`_pclx`, 1.13 M structures, +53.6 k new; median 98 candidates) | 0.840 | 0.620 | −0.005 … −0.004 vs #10 (point ≈ 0) | NOT PROVEN on V-in (no cost); **REJECT on LB** (0.268 vs 0.273) |
| 12 | #10 on pool + top-250k PubChem literature structures (`_pclp`, CID-PMID ranked by article count; +205.6 k new; median 111 candidates) | 0.841 | 0.610 | −0.015 … −0.008 vs #10 (point −0.006 … −0.003) | NOT PROVEN on V-in; **REJECT on LB** (0.252 vs 0.273) |
| 13 | D4 fragment explanation (depth-2 bond cleavage, explained intensity) in the ranker (`ranker-v14-d4`) | V-in 0.807 / V-real C2 0.631 | V-in 0.636 | V-real vs #10: C2 +0.068 | V-in NOT PROVEN (C1 −0.035); V-real C2 clearly better |
| 14 | **Class-3 hedge**: biotransformed spectral analogs as candidates with full features, C3 groups in training (`ranker-v15-gen`) | V-in 0.813 / V-real 0.431 | V-in 0.625 / V-real 0.628 | vs #10: V-real +0.027 … +0.036, V-in +0.024 … +0.067 (CI lower bounds); C3 MRR V-in 0.202, V-real 0.040 | **ACCEPT on both**; LB 0.298 |
| 15 | Generator v2 (42 transforms incl. homologs, carbonyl reduction, acyl/sugar losses; 15 analog seeds) + ranker on 2× training queries (`ranker-v16-gen2`) | V-real 0.433 | V-real 0.636 | vs #14: V-real +0.003 … +0.009 (ACCEPT); V-in lower bounds −0.004 … +0.002; C3 V-real 0.040 → 0.075 | **ACCEPT**; LB 0.309 |
| 16 | v17: + consensus / transform-type / NP-likeness features, depth-3 D4 block, tuned LightGBM (`ranker-v17`) | V-in 0.767 (−0.038) | V-in 0.658 | vs #15: V-real +0.007 … +0.010 (ACCEPT; C3 +0.015); V-in C1 −0.038 (significant) | **REJECT on LB** (0.302 vs 0.309) |
| 17 | Large D3 (d_model 384, 6 layers, 6 epochs, batch 384; truth-vs-decoy 0.868/0.875 vs 0.832) + v16 recipe (`ranker-v18-d3L`) | V-real 0.431 / V-in 0.791 | V-real 0.645 / V-in 0.642 | vs #15: V-real −0.003 … −0.002 (NOT PROVEN); V-in C1 −0.015 | not submitted (guard: V-real must ACCEPT) |
| 18 | D5 learned bond-cleavage likelihood (793 event types from 30 k library spectra), screened alone on Class-2 lists | — | alone: V-in 0.058 (LL) / 0.315 (weighted observed); V-real 0.099 / 0.415 — D4 alone 0.404 / 0.558 | — | dropped at screen (weaker than D4 alone; mean log-likelihood is dominated by event count) |
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

8. **V-in cannot see pool decoys; V-real can (runs `vreal_*`).** V-real = 600 queries from the
   training libraries that no ranker ever trained on (skeleton-level exclusion): 300 COCONUT
   members, 300 structures known only through PubChemLite. Scored under the seven configurations
   the leaderboard has judged, the COCONUT Class-2 stratum falls monotonically with pool size
   (default 0.594 → +PubChemLite 0.548 → +ChEBI/LIPID MAPS 0.529 → +PubChem-lit 0.518), matching the
   leaderboard's verdicts on #11 and #12 that V-in missed (V-in: −0.001). Its level is also close to
   real: LB ≈ 0.5 × coco-C2 reproduces the leaderboard within ±0.01. First version had a design bug —
   a PubChemLite-only truth stayed in every pool because it also carries the training-library
   source bit; fixed (a Class-2 truth is a candidate only through a non-library source).

## Leaderboard

| Date | Submission | Public LB | Note |
|---|---|---|---|
| 2026-09-24 | #7 v13-clean (D1+D2+D3, no membership priors) | **0.258** | two identical submissions scored identically — the pipeline is deterministic on the hidden test |
| 2026-09-24 | #1 interim rule max(D1, D2) — calibration | 0.190 | the learned ranker transfers: +0.068 |
| 2026-09-24 | #3 v1 ranker *with* membership priors — calibration | 0.245 | dropping the priors helped on the real test: +0.013 |
| 2026-09-24 | #9 v13-clean on pool + PubChemLite | 0.267 | V-in predicted ≈ −0.007 (cost only), so real PubChem-only answers are worth ≈ +0.016 |
| 2026-09-25 | #10 ranker retrained on PubChemLite-pool candidate lists (notebook `casmi26-probe-pcl2`) | **0.273** | new best; V-in neutral (C1 +0.010, C2 −0.008, NOT PROVEN) — consistent with better ranking of PubChem-only answers, which V-in cannot see; small, treated as a consistent signal rather than proof |
| 2026-09-25 | #11 = #10 + ChEBI + LIPID MAPS in the pool (`casmi26-probe-pclx`) | 0.268 | −0.005 vs #10 although V-in saw no cost: the 53.6 k extra structures add more decoys than hidden-test answers. Coverage of *training-library* structures (+2.7 pts) does not predict LB gain — PubChemLite had a similar per-structure yield and helped. Pool additions must be judged on the LB, one source at a time |
| 2026-09-25 | #12 = #10 + top-250k PubChem literature structures (`casmi26-probe-pclp`) | 0.252 | −0.021 vs #10. Yield on training-library structures was 37 per 10 k new structures (ChEBI 224, PubChemLite 207). **Conclusion of the pool-coverage work: PubChemLite is the only addition that pays; generic chemistry databases add decoys faster than answers.** Further gains have to come from ranking, not from a bigger pool |
| 2026-09-26 | #14 D4 + Class-3 hedge (`ranker-v15-gen`, notebook `casmi26-probe-gen`) | **0.298** | new best, +0.025 vs #10; V-real predicted +0.027 … +0.051 (point), V-in's C1 loss did not show up |
| 2026-09-26 | #14b = #14 with Class-3 groups weighted ×2 (`ranker-v15-gen-c3w2`, `casmi26-probe-gen2`) | 0.281 | −0.017 vs #14. V-real saw the direction (C2 −0.009, C3 +0.004) but not the size: pushing generated candidates up costs real Class-1/2 molecules more than it gains on Class 3. Keep weight 1 |
| 2026-09-26 | #15 generator v2 + 2× training queries (`ranker-v16-gen2`, `casmi26-probe-gen3`) | **0.309** | new best, +0.011 vs #14; V-real predicted +0.012 … +0.018 (point) |
| 2026-09-27 | #16 v17 (`casmi26-probe-v17`) | 0.302 | −0.007 vs #15. Second loss for a change whose V-real gain came mostly from simulated Class 3 (#14b too). Lessons: (1) V-real's Class-3 simulation is optimistic — its queries have close library analogs, real novel molecules have fewer; discount C3 gains. (2) a significant V-in Class-1 drop is a real warning: V-in is the only instrument-matched set, and it flagged this change |
| 2026-09-28 | #16a v17 without NP-likeness (`casmi26-probe-v17nonp`) | 0.304 | NP-likeness explained the V-in Class-1 loss but not the LB loss |
| 2026-09-28 | #16b v16 features + tuned LightGBM (`casmi26-probe-v17abl`) | 0.305 | even the tuning alone lost 0.004. **All of v17 rejected; v16 (0.309) stays.** V-real gains of +0.005 … +0.010 are below what it can resolve against the leaderboard |

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

**PubChemLite (#9).** V-in rejected the larger pool (C2 −0.035: more decoys, median 96 candidates vs
58), yet the leaderboard rose 0.258 → 0.267. V-in cannot see the benefit by construction — every V-in
answer is a COCONUT structure — so this is the first measured evidence that the hidden test contains
answers reachable only through PubChem. Pool coverage is a real lever; V-in only measures its cost.
Next: retrain the ranker on the larger pool's candidate lists (its list-context features shifted).

## Engineering notes

- Spectrum store: 2.54 M spectra → 0.87 GB memmaps, built in 138 s after replacing the NumPy
  cleaner with a numba version (verified identical on 900 random spectra; ~65× faster).
- Two memory-heavy processes (rank-data builder ~3.5 GB, GPU training ~8.7 GB commit incl. the CUDA
  driver) exceeded the 16 GB laptop's commit headroom; stages now run strictly in sequence.
- D3: 3.9 M-parameter set transformer, 2.0 M training spectra, ~2,800 spectra/s on the RTX 3060
  laptop GPU, 8 min/epoch. Held-out-fold BCE 0.0672 (prior-only ≈ 0.109); the true fingerprint
  out-scores a random decoy 83% of the time. Trained locally: no Kaggle GPU quota used.
- End-to-end submission path on the public file: 643 s for 400 molecules, 0 fallbacks.

## State at pause (2026-09-28)

- **Best submission:** #15 `ranker-v16-gen2`, notebook `casmi26-probe-gen3`, public LB **0.309**
  (train + COCONUT + PubChemLite pool; D1/D2/D3/D4 channels; biotransformation hedge v2; LambdaRank on
  3,000 simulated queries). Second-best distinct submission: #16b `casmi26-probe-v17abl`, 0.305.
- **Reproduce:** `kaggle/submission_driver.py` with `CASMI_RANKER=ranker-v16-gen2`, `CASMI_USE_PUBCHEMLITE=1`
  (Kaggle datasets casmi-code, casmi-models, casmi-coconut, casmi-pubchemlite). Engine options added after v16
  (consensus, NP-likeness, depth-3 D4, large D3) switch on only for rankers trained with them.
- **What did not transfer:** pool growth beyond PubChemLite, Class-3 up-weighting, v17 features, tuned
  LightGBM, large D3. V-real resolves effects of about ±0.01 and overrates simulated Class-3 gains; a
  significant V-in Class-1 drop has been a reliable warning.
