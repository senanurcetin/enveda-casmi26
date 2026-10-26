# Fresh-weight internal chemical-family benchmark

Every fingerprint network and LambdaRank model was fitted from a fresh seeded initialization. The final test was opened only after all rankers, feature lists and learned artifact hashes were sealed. This is an internal query-family benchmark, not newly acquired external/laboratory validation or a competition score.

![NEXUS internal benchmark and paired diagnostics](../docs/assets/advanced-internal-study.png)

[Vector SVG](../docs/assets/advanced-internal-study.svg) · [Standalone PDF](../docs/assets/advanced-internal-study.pdf)

## Predeclared population and fitting lineage

| Partition | Molecules / families |
|---|---:|
| train | 600 |
| dev | 150 |
| test | 250 |

Cross-partition selected metric/query/family overlap: **0**. One representative per generic Murcko/acyclic connectivity family defines a family-balanced population, not prevalence-weighted accuracy. The shared candidate universe and task-defined references remain available; held-out query separation is not reference-free discovery.

| Fresh D3 fold | Training structures | Spectra | First / last training BCE |
|---|---:|---:|---:|
| 0 | 279 | 2996 | 0.6555 / 0.1021 |
| 1 | 321 | 3087 | 0.6545 / 0.0992 |

D3 has six fixed epochs, with no development/test checkpoint selection. Ranker-training queries use the opposite D3 fold. Training BCE describes fitted training data, not held-out fingerprint accuracy. Development data alone select ranker stopping iterations. Freshness concerns fingerprint and ranking weights; the RDKit NP-score is a fixed external chemical prior whose source-training overlap is not independently attested.

## Test MRR@25

| Model | C1, 95% CI | C2, 95% CI | C3, 95% CI |
|---|---:|---:|---:|
| mass_only | 0.3863 [0.3391, 0.4355] | 0.3863 [0.3391, 0.4355] | 0.0000 [0.0000, 0.0000] |
| d1 | 0.4305 [0.3752, 0.4861] | 0.1741 [0.1389, 0.2103] | 0.0000 [0.0000, 0.0000] |
| d1_d2 | 0.4608 [0.4068, 0.5144] | 0.2066 [0.1692, 0.2447] | 0.0000 [0.0000, 0.0000] |
| full | 0.8665 [0.8289, 0.9025] | 0.8388 [0.7971, 0.8782] | 0.0373 [0.0188, 0.0588] |
| no_d1 | 0.8195 [0.7796, 0.8587] | 0.7659 [0.7212, 0.8085] | 0.0315 [0.0158, 0.0508] |
| no_d2 | 0.8099 [0.7677, 0.8517] | 0.7512 [0.7058, 0.7966] | 0.0344 [0.0171, 0.0554] |
| no_d3 | 0.8617 [0.8236, 0.8979] | 0.8262 [0.7840, 0.8661] | 0.0372 [0.0189, 0.0585] |
| no_d4 | 0.7952 [0.7515, 0.8367] | 0.7399 [0.6933, 0.7850] | 0.0306 [0.0150, 0.0483] |
| no_gen_features | 0.7798 [0.7354, 0.8249] | 0.7073 [0.6571, 0.7558] | 0.0565 [0.0321, 0.0834] |
| no_np_score | 0.8713 [0.8342, 0.9059] | 0.8440 [0.8031, 0.8830] | 0.0381 [0.0196, 0.0596] |
| membership_only | 0.4204 [0.3773, 0.4668] | 0.1061 [0.0808, 0.1337] | 0.0054 [0.0023, 0.0094] |

The mass-only C2 baseline scores 0.3863, compared with 0.2066 for D1+D2. The primary comparison was fixed as D1+D2 before execution; its large effect should be interpreted alongside this simple-baseline weakness. No test-selected baseline replaces the declared comparison.

## Single predeclared primary contrast

Unperturbed test C2 MRR@25: fresh full model minus D1+D2 baseline.

**ΔMRR +0.6322, paired 95% interval [+0.5816, +0.6829].** This estimates the difference for the selected internal population; no external-generalization or automatic model-promotion claim follows.

## Joint secondary uncertainty

All 27 prespecified feature/control × novelty contrasts share one approximate max-deviation bootstrap family. Negative effects mean the listed treatment scored below fresh full. All feature ablations are refitted; they share generated candidates. `no_gen_features` removes features, not the candidate generator.

| Secondary comparison | C1 ΔMRR, joint CI | C2 ΔMRR, joint CI | C3 ΔMRR, joint CI |
|---|---:|---:|---:|
| no_d1_minus_full | -0.0471 [-0.1135, +0.0193] | -0.0729 [-0.1393, -0.0065] | -0.0059 [-0.0723, +0.0605] |
| no_d2_minus_full | -0.0566 [-0.1230, +0.0098] | -0.0876 [-0.1540, -0.0212] | -0.0029 [-0.0693, +0.0635] |
| no_d3_minus_full | -0.0048 [-0.0712, +0.0616] | -0.0126 [-0.0790, +0.0538] | -0.0001 [-0.0665, +0.0663] |
| no_d4_minus_full | -0.0713 [-0.1377, -0.0049] | -0.0989 [-0.1653, -0.0326] | -0.0067 [-0.0731, +0.0597] |
| no_gen_features_minus_full | -0.0867 [-0.1531, -0.0203] | -0.1315 [-0.1979, -0.0651] | +0.0192 [-0.0472, +0.0856] |
| no_np_score_minus_full | +0.0047 [-0.0616, +0.0711] | +0.0051 [-0.0613, +0.0715] | +0.0008 [-0.0656, +0.0672] |
| membership_only_minus_full | -0.4461 [-0.5125, -0.3797] | -0.7327 [-0.7991, -0.6663] | -0.0319 [-0.0983, +0.0345] |
| full_blank_minus_none | -0.6056 [-0.6720, -0.5392] | -0.5119 [-0.5783, -0.4455] | -0.0373 [-0.1037, +0.0290] |
| full_intensity_shuffle_minus_none | -0.1414 [-0.2078, -0.0750] | -0.0936 [-0.1600, -0.0272] | -0.0165 [-0.0829, +0.0499] |

## Candidate coverage

| Class | Original-pool truth present | Pool + generation truth present |
|---|---:|---:|
| C1 | 246 / 250 | 246 / 250 |
| C2 | 246 / 250 | 246 / 250 |
| C3 | 0 / 250 | 26 / 250 |

Every selected test query remains in the denominator, including empty candidate lists and recall misses. C3 original-pool truth coverage must be zero. Coverage is an upper bound on attainable ranking recall, not achieved MRR.

## Descriptive C2 error slices

Pointwise subgroup intervals are descriptive; the joint secondary guarantee above does not apply to subgroup selection. Cells below the predeclared minimum of 20 are suppressed. The test population contains 141/250 pluskal_ms2 queries: family balancing does not make sources equally represented. Small timsTOF or library cells cannot establish cross-instrument or cross-library transfer.

| Dimension | Slice | n | D1+D2 MRR | Full MRR | ΔMRR, descriptive CI |
|---|---|---:|---:|---:|---:|
| library | gnps | 56 | 0.1985 | 0.6858 | +0.4872 [+0.3721, +0.5935] |
| library | massbank | 11 | — | — | Suppressed |
| library | mona | 6 | — | — | Suppressed |
| library | msdial | 2 | — | — | Suppressed |
| library | pluskal_ms2 | 141 | 0.1363 | 0.9504 | +0.8141 [+0.7645, +0.8597] |
| library | riken | 7 | — | — | Suppressed |
| library | spectraverse | 27 | 0.4975 | 0.7903 | +0.2928 [+0.1608, +0.4294] |
| adduct | [M+CH2O2-H]- | 3 | — | — | Suppressed |
| adduct | [M+H]+ | 129 | 0.1438 | 0.8360 | +0.6922 [+0.6208, +0.7598] |
| adduct | [M+NH4]+ | 3 | — | — | Suppressed |
| adduct | [M+Na]+ | 2 | — | — | Suppressed |
| adduct | [M-H2O+H]+ | 1 | — | — | Suppressed |
| adduct | [M-H]- | 23 | 0.5438 | 0.8696 | +0.3258 [+0.1840, +0.4775] |
| adduct | mixed | 89 | 0.2063 | 0.8632 | +0.6568 [+0.5803, +0.7325] |
| mode | mixed | 67 | 0.2071 | 0.9024 | +0.6953 [+0.6118, +0.7769] |
| mode | negative | 29 | 0.4611 | 0.8477 | +0.3866 [+0.2584, +0.5219] |
| mode | positive | 154 | 0.1584 | 0.8095 | +0.6511 [+0.5846, +0.7135] |
| mass | 150-300 | 30 | 0.1682 | 0.6674 | +0.4992 [+0.3567, +0.6479] |
| mass | 300-600 | 194 | 0.1837 | 0.8859 | +0.7022 [+0.6476, +0.7554] |
| mass | 600-1170 | 26 | 0.4218 | 0.6855 | +0.2637 [+0.1459, +0.3954] |
| spectra | 1 | 43 | 0.3654 | 0.6740 | +0.3086 [+0.1827, +0.4339] |
| spectra | 2-3 | 28 | 0.2449 | 0.7700 | +0.5251 [+0.3860, +0.6639] |
| spectra | 4-5 | 179 | 0.1624 | 0.8892 | +0.7267 [+0.6718, +0.7777] |
| instrument | other_or_mixed | 250 | 0.2066 | 0.8388 | +0.6322 [+0.5816, +0.6829] |

## Reproducibility and limits

[Protocol](../protocols/advanced-internal-v1.json) · [Model/input/code seals and all aggregate conditions](advanced_study.json) · [Authenticated publication readback](advanced_publication.json)

This population, architecture and training procedure differ from historical v16. Comparing their raw scores does not isolate a model improvement. Fixed reference knowledge and the NP-score prior limit claims of wholly unseen chemistry.

- Internal source data, not newly acquired external spectra or an independent laboratory replication.
- Family-disjoint supervised query fitting does not remove the task-defined retrieval reference corpus.
- Generic Murcko topology and acyclic connectivity are computational grouping rules, not a guarantee of distant chemistry.
- Fixed training sample and small fresh fingerprint architecture trade capacity for auditable zero-cost execution; no state-of-the-art claim.
- Feature ablations reuse a fixed generated candidate universe; no_gen_features does not remove candidate generation.
- Slices are descriptive; subgroup intervals do not establish multiplicity-adjusted subgroup effects.
- No raw spectra, structures, query identifiers, per-query ranks, model weights or submitted predictions are exported.
