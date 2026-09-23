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

## Engineering notes

- Spectrum store: 2.54 M spectra → 0.87 GB memmaps, built in 138 s after replacing the NumPy
  cleaner with a numba version (verified identical on 900 random spectra; ~65× faster).
- Two memory-heavy processes (rank-data builder ~3.5 GB, GPU training ~8.7 GB commit incl. the CUDA
  driver) exceeded the 16 GB laptop's commit headroom; stages now run strictly in sequence.
- D3: 3.9 M-parameter set transformer, 2.0 M training spectra, ~2,800 spectra/s on the RTX 3060
  laptop GPU, 8 min/epoch. Held-out-fold BCE 0.0672 (prior-only ≈ 0.109); the true fingerprint
  out-scores a random decoy 83% of the time. Trained locally: no Kaggle GPU quota used.
- End-to-end submission path on the public file: 643 s for 400 molecules, 0 fallbacks.
