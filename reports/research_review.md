# CASMI26 research review

Reviewed: 2026-10-23. Score snapshot: 2026-10-08.

Noncommercial historical research; no new competition submission or accuracy claim.

These are source-file integrity hashes, not external verification of the claims.

| Experiment | Public MRR@25 | Attribution |
|---|---:|---|
| Library + analog rule | 0.190 | own |
| Ranker without membership priors | 0.258 | own |
| + PubChemLite, retrained | 0.273 | own |
| + Fragments and generation | 0.298 | own |
| Extended generation | 0.309 | own |
| External v4n base | 0.420 | external |
| v4n + own rank-list fusion | 0.422 | fusion |

## Historical public scores — historical

Own pipeline 0.309; fusion 0.422 includes the external v4n solution. These are archived author logs, not a fresh leaderboard readback or final private scores.

Source: [HANDOFF.md](../HANDOFF.md). SHA-256 (UTF-8 / LF): `b527d210d530c427926f0683b4ef6b66ca088b7bec0de6f7c94051680c347a18`.

## Source-library composition — historical

2,539,608 spectra and 275,810 structure keys in the committed profile. Counts are not an independent audit of the current source dataset.

Source: [reports/data_profile.json](../reports/data_profile.json). SHA-256 (UTF-8 / LF): `3c6f3ee71855a1b3045b1e821e3182114df30dfbbd1d29916eaf705fd51a5ffc`.

## Same-library C1 exclusion — implemented

Code excludes complete same-library reference groups, including unsampled equivalent-structure spectra. Synthetic regressions exist; historical scores were not recomputed.

Source: [docs/CASMI26_RESEARCH_ASSESSMENT.md](../docs/CASMI26_RESEARCH_ASSESSMENT.md). SHA-256 (UTF-8 / LF): `498aa860f51af1915fc9956fdd3ddb00132e09555d189d21a2ada3e72424c2dd`.

## Fixed-grid simultaneous uncertainty — implemented

Paired bootstrap bands account for a predeclared novelty-mixture grid. Repeated model selection and domain shift remain outside that guarantee.

Source: [EXPERIMENTS.md](../EXPERIMENTS.md). SHA-256 (UTF-8 / LF): `ab3d140294c7d7065d7c07b5014981301b605e9bd578a36bfea3db8573e15ac9`.

## Full corrected-protocol revalidation — pending

Requires licensed data, immutable input hashes, verified fingerprints, fresh analog indices and paired per-molecule outputs. No completed full-data rerun is asserted.

Source: [docs/RESEARCH_WORK_PROGRAM.md](../docs/RESEARCH_WORK_PROGRAM.md). SHA-256 (UTF-8 / LF): `3bd5ac3c439b0b10383be6a42b5a459dd1f67fcc53147118a4ca9a54753f5136`.

## Exact train / public-example overlap — pending

The archived 1,213 coarse-key matches are insufficient to establish exact peak-and-intensity identity. Recompute content-level overlap before claiming complete leakage.

Source: [src/casmi/validate/spectrum_overlap.py](../src/casmi/validate/spectrum_overlap.py). SHA-256 (UTF-8 / LF): `ed06679f55192f44b1f8db3596ec4be8ff279fce2e813544e5464e015b984aec`.
