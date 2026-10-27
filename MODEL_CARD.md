# Model card — CASMI26 retrieval and ranking

**Owner:** Senanur Cetin / NEXUS. **Scope:** noncommercial molecular identification research. **Reviewed:** 26 October 2026.

## Intended use

Rank up to 25 candidate molecular structures per unknown molecule from MS/MS spectra. Evaluation uses reciprocal rank of the first tautomer-canonical InChIKey14 match, truncated at 25. Invalid or duplicate candidates consume their original rank slots; evaluation never compresses them to improve the score.

## Model boundaries

| System | Components | Recorded public MRR@25 |
|---|---|---:|
| Own pipeline `ranker-v16-gen2` | ±10 ppm candidate retrieval, D1 library, D2 analog, D3 fingerprint, D4 fragments, generated candidates, LambdaRank | 0.309 |
| External v4n base | Third-party fusion solution attributed in the dated experiment log | 0.420 |
| `casmi26-fusion-ours-d` | External base plus own list as a gated reciprocal-rank channel | 0.422 |

These are historical author-reported public results recorded by 8 October 2026. The registry links them to [source logs](reports/research_evidence.json). No final private score, present rank, independent server attestation or statistical significance for +0.002 is asserted. No trained model is served by the dashboard.

## Evaluation and limitations

V-in contains 250 instrument-matched natural-product structures; V-real contains 600 harder held-out library structures. They underrepresent unknown novelty mixtures and instrument/library shift. Historical fusion-weight ordering did not transfer consistently to the public leaderboard.

Later audits corrected sampled-query C1 leakage, D2 eligibility ordering, identity/provenance and numerical edge cases. Historical leaderboard scores were **not** replaced. A separate [completed retrospective study](reports/revalidation_review.md) now evaluates all 250 eligible V-in molecules through the shared ranking path using rebuilt full-source artifacts and frozen original weights. D1 / D1+D2 / v16 MRR@25 is respectively **0.9033 / 0.9040 / 0.8324 (C1)**, **0.1402 / 0.2272 / 0.6241 (C2)** and **0 / 0 / 0.2374 (C3)**. These are selected-population diagnostics, not new independent generalization estimates. Standard/tautomer keys do not establish scaffold independence.

Blank fragments lower v16 C2 MRR by 0.4972 and C3 by 0.2374, but improve simple-model C2 ranking. That counter-result shows that a negative control is not a universal pass/fail certification: residual similarity can rank decoys above mass-error ordering. Full paired intervals, intensity-shuffle controls and model contrasts are retained in the report. That retrospective study fitted no new weights; complete original training lineage is not independently attested.

Current paired comparisons require matching query identities, spectra signatures and source manifests; D2 decision-grade comparisons also require verified fingerprint and fresh-index evidence. Simultaneous bands cover only the fixed class-mixture grid, not adaptive experiment selection. [Research assessment](docs/research-assessment.md).

## Fresh supervised weights and sealed internal families

The [completed fresh-training study](reports/advanced_review.md) uses 600/150/250 train/development/test molecules, one per conservative generic chemical family, with no selected query/metric/family crossing supervised partitions. Two new D3 networks supply cross-fitted training features; eight prespecified three-seed rankers and all feature lists are sealed before test features. The exact executed source, 26 learned-file hashes and public completion are [authenticated](reports/advanced_publication.json).

Fresh full C1/C2/C3 MRR@25 is **0.8665/0.8388/0.0373**. The predeclared C2 full-minus-D1+D2 contrast is **+0.6322 [+0.5816, +0.6829]**. Among 27 jointly assessed secondary contrasts, C2 D1/D2/D4/generator-feature removals reduce performance; D3 and NP-score effects remain unresolved. Blank/shuffle C2 effects are -0.5119/-0.0936. C3 candidate coverage is only 26/250. No automatic promotion follows from testing several variants.

This is an internal query-family benchmark. Evaluation structures may occur as negatives in training candidate lists and task-defined retrieval references remain shared. The fixed NP-score prior is not freshly fitted. New population/training/architecture prevent a causal before/after comparison with historical v16. The retrospective source snapshot also retains its old connectivity-only generated cache; the representation-sensitive correction and outcome-blind amendment apply to the fresh study. Independent external spectra/laboratory replication remain absent.

## Operational and scientific failure modes

The [external/C3 follow-on](docs/studies/external-source-validation.md) adds an opt-in legacy-preserving generator budget and requires a fresh sealed experiment before claiming improvement. Completed MassBank and two-library GNPS source scans leave zero and seven novel molecules respectively; neither meets the declared minimum of 30. No external model accuracy is established. Current isotope-aware connectivity filtering is corrected; completed historical source snapshots retain their original implementation and are not silently rerun.

The [aggregate-only loss decomposition](reports/candidate_loss_review.md) shows 224/250 C3 queries lack the correct candidate. Of the 26 covered molecules, 24 reach top 25 and 6 rank first. The existing-pool oracle MRR ceiling is 0.1040, versus observed 0.0373; this hypothetical ceiling cannot be presented as achieved accuracy. Conditional MRR on the 26 covered molecules is descriptive and does not replace population MRR. No conditional uncertainty or causal generation improvement can be recovered from aggregate-only evidence.

Missing true structures in the candidate pool impose a hard recall ceiling. Isomers share masses and may have similar spectra. Adduct errors, spectral contamination, mode mismatches, limited cross-instrument references, generated-structure identity errors and unknown chemistry can all produce confident-looking but incorrect rankings. Rank scores are not calibrated identification probabilities.

No clinical, toxicological, regulatory, synthesis or other consequential chemical decision should rely on this project. Source code and source datasets have separate licensing; third-party solutions retain their original attribution and terms. [Data card](DATA_CARD.md).
