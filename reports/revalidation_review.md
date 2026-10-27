# Frozen-model retrospective revalidation

Executed on the designated V-in population with rebuilt full-source artifacts. These are retrospective diagnostics, not new public/private competition scores or an untouched generalization benchmark.

- Evaluation denominator: **250 molecules**, retained across all models, novelty classes and controls.
- Mass-conflicting identities quarantined: **1,710**, covering 3,891 provider rows.
- Unusable selected candidate identities quarantined: **13**.
- Retained fingerprint-verified candidate structures: **1,071,863**.
- Fresh D2 reference spectra: **1,062,995**; no index cache used.
- Complete own-library exclusions add **287 source rows** beyond the selected inputs for **74 queries**.
- Queries whose standard truth key is absent from the candidate pool: **1**; none removed from scoring.

## Unperturbed MRR@25

| Model | C1, 95% CI | C2, 95% CI | C3, 95% CI |
|---|---:|---:|---:|
| d1 | 0.9033 [0.8731, 0.9314] | 0.1402 [0.1142, 0.1678] | 0.0000 [0.0000, 0.0000] |
| d1_d2 | 0.9040 [0.8738, 0.9321] | 0.2272 [0.1913, 0.2658] | 0.0000 [0.0000, 0.0000] |
| frozen_v16 | 0.8324 [0.7891, 0.8705] | 0.6241 [0.5753, 0.6694] | 0.2374 [0.1903, 0.2823] |

Absolute intervals use 2,000 molecule bootstrap draws. C1 excludes all equivalent same-library references; C2 excludes equivalent references across libraries; C3 additionally removes candidate identities before generation.

## Paired spectral controls

Negative values mean the perturbed condition scored lower. Intervals are descriptive molecule-paired bootstrap intervals (4,000 draws), not model-promotion tests. Blank spectra retain precursor and acquisition information; intensity shuffling retains fragment positions.

| Model | Control | C1 delta, 95% CI | C2 delta, 95% CI | C3 delta, 95% CI |
|---|---|---:|---:|---:|
| d1 | blank | -0.5981 [-0.6449, -0.5495] | +0.1650 [+0.1351, +0.1966] | +0.0000 [+0.0000, +0.0000] |
| d1 | intensity_shuffle | -0.1235 [-0.1647, -0.0856] | +0.0017 [-0.0054, +0.0095] | +0.0000 [+0.0000, +0.0000] |
| d1_d2 | blank | -0.5988 [-0.6457, -0.5501] | +0.0780 [+0.0349, +0.1236] | +0.0000 [+0.0000, +0.0000] |
| d1_d2 | intensity_shuffle | -0.1240 [-0.1653, -0.0860] | -0.1258 [-0.1573, -0.0967] | +0.0000 [+0.0000, +0.0000] |
| frozen_v16 | blank | -0.7055 [-0.7499, -0.6555] | -0.4972 [-0.5463, -0.4451] | -0.2374 [-0.2814, -0.1919] |
| frozen_v16 | intensity_shuffle | -0.5732 [-0.6246, -0.5194] | -0.2094 [-0.2620, -0.1570] | -0.1886 [-0.2340, -0.1422] |

## Paired unperturbed model contrasts

| Comparison | C1 delta, 95% CI | C2 delta, 95% CI | C3 delta, 95% CI |
|---|---:|---:|---:|
| d1_d2_minus_d1 | +0.0007 [+0.0000, +0.0020] | +0.0870 [+0.0566, +0.1183] | +0.0000 [+0.0000, +0.0000] |
| frozen_v16_minus_d1_d2 | -0.0716 [-0.1075, -0.0386] | +0.3969 [+0.3495, +0.4427] | +0.2374 [+0.1919, +0.2814] |

## Interpretation

In this selected population the frozen learned ranker improves C2/C3 while scoring below D1+D2 in C1. Removing fragments substantially reduces its performance. This tradeoff is class-dependent and does not establish universal model superiority.

The simple models improve in C2 with blank spectra. With zero similarity scores their ranker falls back to mass-error tie-breaking; unperturbed residual similarities can instead promote decoys after truth references are excluded. This is an adverse diagnostic result, not a recommendation to discard spectra or a universal control pass.

The machine-readable report also retains the prespecified fixed-grid mixture bands. Neither pointwise intervals nor grid bands correct historical model selection or domain shift.

## Provenance and limits

[Executed aggregate and input/model/code hashes](revalidation_study.json) · [Authenticated Kaggle readback](revalidation_publication.json) · [Protocol](../docs/studies/frozen-model-revalidation.md)

- Retrospective designated V-in evaluation, not a new untouched external holdout or competition score.
- Frozen learned weights were not retrained; their complete training-data lineage is not independently attested.
- Cohort was used in historical research; intervals do not remove adaptive model-selection or domain-shift bias.
- Controls diagnose fragment reliance; intensity shuffling retains peak positions and precursor information.
- Mixture intervals cover only the predeclared grid; no ACCEPT/REJECT model-promotion claim.
- No raw spectra, molecular structures, query identifiers, ranks or model files are exported.
