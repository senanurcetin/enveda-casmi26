# Technical reviewer brief

**Question:** How can a retrieval-and-ranking system identify natural products from tandem mass spectra when reference-library coverage is incomplete?

**Contribution:** An original mass-constrained candidate pipeline combining library spectra, mass-shifted analog propagation, learned fingerprint likelihood, bond-cleavage features and LambdaRank. Candidate generation extends database retrieval with biotransformations. A separate experiment adds this pipeline to an attributed third-party fusion base.

## Quantitative evidence

| Measurement | Recorded value | Interpretation |
|---|---:|---|
| Training spectra / structure keys | 2,539,608 / 275,810 | Historical profile; source libraries overlap in structure identity |
| Own pipeline public MRR@25 | 0.309 | Archived competition result |
| External v4n / fused public MRR@25 | 0.420 / 0.422 | Fusion difference is not independently validated generalization |
| PubChemLite pool / larger literature pool | 0.273 / 0.252 | More candidates increased decoy costs in the recorded experiments |

Sources: [data profile](../reports/data_profile.json), [experiment log](../EXPERIMENTS.md), [handoff](../HANDOFF.md). No current leaderboard rank or final private score is claimed.

## Three substantive findings

1. **Simulated missingness can become the answer key.** Reference counts and source-membership priors distinguish artificially withheld answers from decoys. Removing those priors improved a recorded public comparison from 0.245 to 0.258. That single comparison is descriptive, not a universal causal effect.
2. **Retrieval coverage and ranking compete.** New databases can contain more answers while adding many mass-compatible isomers. Pool evaluation therefore needs both answer coverage and candidate-window density.
3. **Validation integrity is a separate deliverable.** Later audits found sampled-query C1 reference leakage and prefilter eligibility problems. Corrections have regression coverage, but their real-data effect remains unmeasured. Historic validation values must not be reused as corrected-protocol estimates.

## What a reviewer can inspect

- Molecular identity, ppm filtering and immutable stores in `src/casmi/chem`, `pool` and `spectra`.
- Channel composition and feature ranking in `src/casmi/channels`, `rank` and `engine.py`.
- Query pairing, provenance, negative controls and uncertainty in `src/casmi/validate`.
- Negative experiments and abandoned hypotheses in the original logs.
- A static NEXUS research dashboard that runs without credentials or cloud queries.

This is research software, not a chemical identification service, clinical assay or commercial deployment. The portfolio demonstrates experimental judgment and engineering controls alongside the limitations of the measured results.
