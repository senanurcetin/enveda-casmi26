# CASMI26 research review

Reviewed: 2026-10-24. Score snapshot: 2026-10-08.

Noncommercial research: historical model scores and executed data-quality evidence; no new competition submission or accuracy claim.

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

Archived profile: 2,539,608 spectra and 275,810 structure keys. The new content audit rescanned the spectrum-row count; distinct-structure totals and library partitions were not recomputed.

Source: [reports/data_profile.json](../reports/data_profile.json). SHA-256 (UTF-8 / LF): `3c6f3ee71855a1b3045b1e821e3182114df30dfbbd1d29916eaf705fd51a5ffc`.

## Same-library C1 exclusion — implemented

Code excludes complete same-library reference groups, including unsampled equivalent-structure spectra. Synthetic regressions exist; historical scores were not recomputed.

Source: [docs/CASMI26_RESEARCH_ASSESSMENT.md](../docs/CASMI26_RESEARCH_ASSESSMENT.md). SHA-256 (UTF-8 / LF): `498aa860f51af1915fc9956fdd3ddb00132e09555d189d21a2ada3e72424c2dd`.

## Fixed-grid simultaneous uncertainty — implemented

Paired bootstrap bands account for a predeclared novelty-mixture grid. Repeated model selection and domain shift remain outside that guarantee.

Source: [EXPERIMENTS.md](../EXPERIMENTS.md). SHA-256 (UTF-8 / LF): `ab3d140294c7d7065d7c07b5014981301b605e9bd578a36bfea3db8573e15ac9`.

## Full corrected-protocol revalidation — pending

Requires licensed data, immutable input hashes, verified fingerprints, fresh analog indices and paired per-molecule outputs. No completed full-data rerun is asserted.

Source: [docs/RESEARCH_WORK_PROGRAM.md](../docs/RESEARCH_WORK_PROGRAM.md). SHA-256 (UTF-8 / LF): `e21f53c21404c68b03903305a82883d429ab33086d122f8cc708f0a5e53af8dc`.

## Exact train / public-example overlap — measured

Executed on Kaggle: all 2,539,608 training rows scanned; 1,213/1,213 public spectra and 400/400 molecules have exact precursor, full m/z and intensity matches, also consistent in adduct/mode. Two coarse candidate pairs were rejected. This concerns public examples, not hidden evaluation.

Source: [reports/public_overlap_audit.json](../reports/public_overlap_audit.json). SHA-256 (UTF-8 / LF): `1884f12c133da90612ab9312ad267a8d533fed9057e08a972775fe89ab2be895`.
