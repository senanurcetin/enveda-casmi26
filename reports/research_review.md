# CASMI26 research review

Reviewed: 2026-10-25. Score snapshot: 2026-10-08.

Noncommercial research: historical competition scores, executed data audits and retrospective frozen-model/control evidence; no new competition submission or untouched-holdout accuracy claim.

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

Source: [EXPERIMENTS.md](../EXPERIMENTS.md). SHA-256 (UTF-8 / LF): `3f3d5bc255dd1b00a57971f288b269d72327c70c44a258b950204daca61976c2`.

## Source-library composition — historical

Archived profile: 2,539,608 spectra and 275,810 structure keys. The new content audit rescanned the spectrum-row count; distinct-structure totals and library partitions were not recomputed.

Source: [reports/data_profile.json](../reports/data_profile.json). SHA-256 (UTF-8 / LF): `3c6f3ee71855a1b3045b1e821e3182114df30dfbbd1d29916eaf705fd51a5ffc`.

## Same-library C1 exclusion — implemented

Complete same-library reference groups, including unsampled equivalent spectra, are excluded. The corrected path has now run on all 250 V-in molecules; historical leaderboard scores are unchanged and the isolated causal effect of this fix was not estimated.

Source: [docs/CASMI26_RESEARCH_ASSESSMENT.md](../docs/CASMI26_RESEARCH_ASSESSMENT.md). SHA-256 (UTF-8 / LF): `53e1e17ddc5e81b8460db48ee9f5c38283923aaa75b7b23c7e03395200519ce2`.

## Fixed-grid simultaneous uncertainty — implemented

Paired bootstrap bands account for a predeclared novelty-mixture grid. Repeated model selection and domain shift remain outside that guarantee.

Source: [EXPERIMENTS.md](../EXPERIMENTS.md). SHA-256 (UTF-8 / LF): `3f3d5bc255dd1b00a57971f288b269d72327c70c44a258b950204daca61976c2`.

## Corrected-protocol retrospective revalidation — measured

Completed Kaggle run with matching code/artifact hashes: all 250 V-in molecules, rebuilt full-source artifacts, verified fingerprints and fresh D2 index. Frozen v16 C1/C2/C3 MRR: 0.8324/0.6241/0.2374. Historical selection and incomplete original weight lineage limit interpretation.

Source: [reports/revalidation_study.json](../reports/revalidation_study.json). SHA-256 (UTF-8 / LF): `32a2211c52f1f0c5f85d7556e84e12df0a917a140a50597be1ef055d8fe7055e`.

## Exact train / public-example overlap — measured

Executed on Kaggle: all 2,539,608 training rows scanned; 1,213/1,213 public spectra and 400/400 molecules have exact precursor, full m/z and intensity matches, also consistent in adduct/mode. Two coarse candidate pairs were rejected. This concerns public examples, not hidden evaluation.

Source: [reports/public_overlap_audit.json](../reports/public_overlap_audit.json). SHA-256 (UTF-8 / LF): `1884f12c133da90612ab9312ad267a8d533fed9057e08a972775fe89ab2be895`.

## Molecule-paired spectral controls — measured

Blank and intensity-shuffle treatments completed for all three models and novelty classes. Blank fragments reduce v16 C2 MRR by 0.4972, but improve simple-model C2 ranking. Controls expose dependence and failure modes; they do not certify generalization.

Source: [reports/revalidation_review.md](../reports/revalidation_review.md). SHA-256 (UTF-8 / LF): `ff669ac2152f2df7c6cc4fb1b84ab9560cf2cc6fda1ffe70e6d07e5b0fac92db`.

## Executed candidate-source quarantine — measured

1,710 mass-conflicting identities (3,891 provider rows) and 13 unparseable candidates quarantined. 1,071,863 structures fingerprint-verified; all 250 query molecules retained, including one with an absent standard truth key.

Source: [reports/revalidation_study.json](../reports/revalidation_study.json). SHA-256 (UTF-8 / LF): `32a2211c52f1f0c5f85d7556e84e12df0a917a140a50597be1ef055d8fe7055e`.

## Untouched external holdout — pending

The executed study is retrospective. Fresh training with independently attested lineage and an untouched external/scaffold holdout remains a separate research requirement; no such generalization claim is made.

Source: [MODEL_CARD.md](../MODEL_CARD.md). SHA-256 (UTF-8 / LF): `5e4039b4ba2868693f116129e36ed7182c96dfa764e1b4f51de64cff4aee44ba`.
