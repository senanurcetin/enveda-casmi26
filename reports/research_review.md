# CASMI26 research review

Reviewed: 2026-10-30. Score snapshot: 2026-10-08.

Noncommercial research: historical scores, executed audits, frozen retrospective controls, fresh internal-family benchmarks and a separately registered purged-reference external evaluation; no competition submission or independent laboratory replication.

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

## Same-library C1 exclusion — measured

Complete same-library reference groups, including unsampled equivalent spectra, are excluded. The corrected path has now run on all 250 V-in molecules; historical leaderboard scores are unchanged and the isolated causal effect of this fix was not estimated.

Source: [reports/revalidation_study.json](../reports/revalidation_study.json). SHA-256 (UTF-8 / LF): `32a2211c52f1f0c5f85d7556e84e12df0a917a140a50597be1ef055d8fe7055e`.

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

Source: [reports/revalidation_review.md](../reports/revalidation_review.md). SHA-256 (UTF-8 / LF): `492b09c3dceeeff55b52c06cf15dc6cbb78ef51d648fe26a743dd03dd22104b9`.

## Executed candidate-source quarantine — measured

1,710 mass-conflicting identities (3,891 provider rows) and 13 unparseable candidates quarantined. 1,071,863 structures fingerprint-verified; all 250 query molecules retained, including one with an absent standard truth key.

Source: [reports/revalidation_study.json](../reports/revalidation_study.json). SHA-256 (UTF-8 / LF): `32a2211c52f1f0c5f85d7556e84e12df0a917a140a50597be1ef055d8fe7055e`.

## Independent laboratory replication remains unavailable — pending

Completed public-source external-family holdout relative to a rebuilt purged training corpus is reported separately. New laboratory acquisition and independently authenticated identity/measurement replication remain unavailable; shared candidate knowledge limits unseen-chemistry claims.

Source: [MODEL_CARD.md](../MODEL_CARD.md). SHA-256 (UTF-8 / LF): `aa2df39a86ee3f6df1d2a8f926e13d2fd1ecfc64f5ff9242df295b9a10a01803`.

## Fresh training and sealed internal query families — measured

Completed source-matched Kaggle run: 600/150/250 training/development/test queries, two fresh fingerprint models and eight three-seed rankers, all 26 learned artifacts sealed before test features. Full C1/C2/C3 MRR 0.8665/0.8388/0.0373; primary C2 delta +0.6322 [0.5816, 0.6829]. Internal supervised query-family separation, not external validation.

Source: [reports/advanced_study.json](../reports/advanced_study.json). SHA-256 (UTF-8 / LF): `9f6ebf28dccf433dc4a8077f907535d11d775ac5e4700ca0760fc9b18f2c8ee9`.

## Refitted ablations and joint secondary controls — measured

27 prespecified contrasts share paired joint intervals. C2 D1/D2/D4/generator-feature removals reduce performance; D3 and NP-score effects remain unresolved. Blank/shuffle C2 deltas -0.5119/-0.0936. C3 truth coverage only 26/250.

Source: [reports/advanced_review.md](../reports/advanced_review.md). SHA-256 (UTF-8 / LF): `d229b13eb1b274b4222bdf067e56ca48a565aeece693b20803dafaea8e005fa4`.

## Hash-locked and installed scientific environment — implemented

Python 3.12 exact direct/transitive SHA-256 lock with 69 conditional entries; 50 installed Windows packages attested, compatible, and full-suite JUnit evidence recorded. Kaggle exports its separate exact package inventory.

Source: [reports/research_environment_validation.json](../reports/research_environment_validation.json). SHA-256 (UTF-8 / LF): `6e5eaafa445e70b1ba647de8dcdbb7a3007f0fbf51863f23304a7475ccd06ccf`.

## Candidate coverage and ranking loss accounting — measured

Post hoc algebraic decomposition of the authenticated frozen internal test. C3: 224/250 truths absent, 24/26 covered truths in top 25, 6 rank first. Existing-list oracle ceiling 0.1040 MRR versus observed 0.0373. Conditional metrics retain explicit denominators; no new validation, causal gain or model promotion.

Source: [reports/candidate_loss_budget.json](../reports/candidate_loss_budget.json). SHA-256 (UTF-8 / LF): `8df9ff35b89ba55de2617ee385c90bf237ed838f8762e28600891dac093d06d6`.

## Fail-closed numerical rank contracts — implemented

Summary and paired bootstrap reject noninteger, negative, out-of-top-25, nonfinite and malformed rank arrays before scoring. Fixed-seed valid results remain unchanged. Completed experiment source bundles remain frozen.

Source: [src/casmi/validate/rank_inputs.py](../src/casmi/validate/rank_inputs.py). SHA-256 (UTF-8 / LF): `c6faf51be68b80ba2f50bc5324ab802ffd1ebd9134eb71005d8dcb9123f8de62`.

## Pinned MassBank external-source feasibility — measured

Full 2,539,608-row training scan excludes all 491 eligible CC0 molecules by standard identity; zero novel molecules. Authenticated completed source. This does not establish external model accuracy.

Source: [reports/external_audit.json](../reports/external_audit.json). SHA-256 (UTF-8 / LF): `22ea8a8b213684e3868cef6548dc18c894ae0dd480f400b7c9136e27acde499e`.

## Two-library GNPS external-source feasibility — measured

Full training scan leaves 7 novel molecules among 85 eligible identities; below the predeclared minimum of 30. No model performance was inspected. The isotope-safe second execution is the authenticated source.

Source: [reports/gnps_external_audit.json](../reports/gnps_external_audit.json). SHA-256 (UTF-8 / LF): `0bd70bc61f91351fc4f18286ebc28d14b130a33aecc87318d145179bb3c34f8a`.

## Legacy-preserving C3 candidate expansion — measured

Preserved legacy candidates and zero observed losses. Internal C3 coverage 13/150 to 17/150; paired change +0.0267 [0.0067, 0.0533]. Joint external ranking improvements remain unresolved; external C3 coverage only 23/150.

Source: [reports/purged_external_study.json](../reports/purged_external_study.json). SHA-256 (UTF-8 / LF): `4a42c640699cb120f60f74f27de6301052a949091716fedff3e3bf99382bc59c`.

## Expanded GNPS original-snapshot novelty feasibility — measured

Authenticated full-training scan leaves the same seven novel molecules among 936 eligible identities, below the declared minimum of 30. Expanded GNPS includes the earlier two-library source, so these stages are not independent replications.

Source: [reports/gnps_expanded_external_audit.json](../reports/gnps_expanded_external_audit.json). SHA-256 (UTF-8 / LF): `2e5664b597a0662a0f55204b7d4b0a5558e32efdfd9637d2edcd7a62c79d124f`.

## Separately registered family-purged external evaluation — measured

The authenticated family-purged study removed 261,826 of 2,539,608 training spectra and fitted eight fresh learned artifacts. On 150 internal test molecules, expanded C3 truth coverage rose from 13 to 17 (paired change +0.0267, 95% CI [+0.0067, +0.0533]). On 150 external molecules, expanded C2/C3 MRR@25 was 0.4454/0.0466; expanded-minus-legacy ranking effects remain unresolved in the joint 11-contrast intervals. External C3 truth coverage remains only 23/150. Original-snapshot novelty gates remain unsuccessful; this is generalization relative to a rebuilt purged corpus, not independent laboratory replication.

Source: [reports/purged_external_study.json](../reports/purged_external_study.json). SHA-256 (UTF-8 / LF): `4a42c640699cb120f60f74f27de6301052a949091716fedff3e3bf99382bc59c`.

## Optional exact-assignment cosine reference baseline — implemented

Bounded one-to-one maximum intensity-product assignment; default D1 entropy remains unchanged. Exhaustive synthetic assignment, normalization and library-exclusion tests validate implementation only. No data-backed performance advantage, new C3 coverage or state-of-the-art claim is established.

Source: [src/casmi/spectra/optimal_cosine.py](../src/casmi/spectra/optimal_cosine.py). SHA-256 (UTF-8 / LF): `f2657a6902c9c10c57b593c4712df7c0470e36de82d67dac6455be0dedbe25af`.
