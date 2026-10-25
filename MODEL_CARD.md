# Model card — CASMI26 retrieval and ranking

**Owner:** Senanur Cetin / NEXUS. **Scope:** noncommercial molecular identification research. **Reviewed:** 25 October 2026.

## Intended use

Rank up to 25 candidate molecular structures per unknown molecule from MS/MS spectra. Evaluation uses reciprocal rank of the first tautomer-canonical InChIKey14 match, truncated at 25. Invalid or duplicate candidates consume their original rank slots; evaluation never compresses them to improve the score.

## Model boundaries

| System | Components | Recorded public MRR@25 |
|---|---|---:|
| Own pipeline `ranker-v16-gen2` | ±10 ppm candidate retrieval, D1 library, D2 analog, D3 fingerprint, D4 fragments, generated candidates, LambdaRank | 0.309 |
| External v4n base | Third-party fusion solution attributed in the handoff | 0.420 |
| `casmi26-fusion-ours-d` | External base plus own list as a gated reciprocal-rank channel | 0.422 |

These are historical author-reported public results recorded by 8 October 2026. The registry links them to [source logs](reports/research_evidence.json). No final private score, present rank, independent server attestation or statistical significance for +0.002 is asserted. No trained model is served by the dashboard.

## Evaluation and limitations

V-in contains 250 instrument-matched natural-product structures; V-real contains 600 harder held-out library structures. They underrepresent unknown novelty mixtures and instrument/library shift. Historical fusion-weight ordering did not transfer consistently to the public leaderboard.

Later audits corrected sampled-query C1 leakage, D2 eligibility ordering, identity/provenance and numerical edge cases. Historical leaderboard scores were **not** replaced. A separate [completed retrospective study](reports/revalidation_review.md) now evaluates all 250 eligible V-in molecules through the shared ranking path using rebuilt full-source artifacts and frozen original weights. D1 / D1+D2 / v16 MRR@25 is respectively **0.9033 / 0.9040 / 0.8324 (C1)**, **0.1402 / 0.2272 / 0.6241 (C2)** and **0 / 0 / 0.2374 (C3)**. These are selected-population diagnostics, not new independent generalization estimates. Standard/tautomer keys do not establish scaffold independence.

Blank fragments lower v16 C2 MRR by 0.4972 and C3 by 0.2374, but improve simple-model C2 ranking. That counter-result shows that a negative control is not a universal pass/fail certification: residual similarity can rank decoys above mass-error ordering. Full paired intervals, intensity-shuffle controls and model contrasts are retained in the report. No new weights were trained; complete original training lineage is not independently attested.

Current paired comparisons require matching query identities, spectra signatures and source manifests; D2 decision-grade comparisons also require verified fingerprint and fresh-index evidence. Simultaneous bands cover only the fixed class-mixture grid, not adaptive experiment selection. [Research assessment](docs/CASMI26_RESEARCH_ASSESSMENT.md).

## Failure modes

Missing true structures in the candidate pool impose a hard recall ceiling. Isomers share masses and may have similar spectra. Adduct errors, spectral contamination, mode mismatches, limited cross-instrument references, generated-structure identity errors and unknown chemistry can all produce confident-looking but incorrect rankings. Rank scores are not calibrated identification probabilities.

No clinical, toxicological, regulatory, synthesis or other consequential chemical decision should rely on this project. Source code and source datasets have separate licensing; third-party solutions retain their original attribution and terms. [Data card](DATA_CARD.md).
