# Portfolio maintenance verification — 26 October 2026

## Aggregate loss diagnosis and numerical scoring contracts

The [candidate/ranking loss report](../reports/candidate_loss_review.md) derives exact accounting from the authenticated, unchanged fresh-study artifact. All 11 models and three novelty classes are retained. C3 has 224 absent truths, 24 top-25 hits among 26 covered truths and 6 top-1 hits; fixed-list oracle MRR is only 0.1040. Conditional means are labeled descriptive, with no reconstructed bootstrap interval, causal improvement claim, retraining or test-selected model. Derived JSON binds the source artifact and notebook hashes. NEXUS PNG/SVG/PDF figures were inspected as generated artifacts.

Code review reproduced two numerical input failures: `summarize([26])` scored a rank outside the top-25 contract, and paired bootstrap silently treated a NaN rank as a miss. The current summary and paired-comparison paths now require nonempty one-dimensional integer arrays in 0–25 and valid integer bootstrap settings. Fixed-seed valid scores/intervals are unchanged. Completed experiment notebooks preserve their original source snapshots; the recorded runner generates bounded integer ranks, so this correction does not revise historical metrics.

The updated locked suite passed **335 tests and 132 subtests, with 2 Windows symlink skips** (469 JUnit cases, zero failures/errors). The standalone stdlib suite passed 103 tests with the same 2 skips. Coverage tests reject impossible molecule counts, inconsistent MRR/cutoff summaries, altered source provenance and zero-coverage conditional division. CI regenerates the derived reports and runs the numerical rank-contract regressions. [Current environment/JUnit attestation](../reports/research_environment_validation.json).

## Completed advanced internal research extension

The fifth owned [Kaggle study](https://www.kaggle.com/code/senanuretin/casmi26-fresh-training-family-holdout) is public and **COMPLETE**, with exact notebook source, input attachments and aggregate-artifact SHA-256 [authenticated on 26 October](../reports/advanced_publication.json). It fits two new fingerprint networks and eight three-seed rankers on 600 training / 150 development queries. All 26 learned artifacts and feature lists are sealed before final-test features for 250 query families. Selected supervised query/metric/family overlap is zero. No archived model dataset, raw record export, competition submission or production deployment was used.

The [full result report](../reports/advanced_review.md) retains every model, all three novelty classes, spectral controls, 27 joint secondary contrasts, coverage ceilings and descriptive source/adduct/mode/mass/spectrum-count/instrument slices. C2 primary delta is +0.6322 [paired 95% CI +0.5816, +0.6829]. C3 remains limited (MRR 0.0373, truth coverage 26/250); D3 and NP-score incremental benefits remain unresolved. Source concentration, shared candidate/reference knowledge and the fixed NP-score prior limit interpretation. This is internal validation, without external acquired spectra or laboratory replication.

The full hash-locked Python 3.12 suite passed **314 tests and 123 subtests, with 2 Windows symlink skips**; JUnit contains 439 cases, zero failures/errors. All 50 applicable installed package versions match the 69-entry conditional SHA-256 lock and pass dependency compatibility checks. The [environment/JUnit attestation](../reports/research_environment_validation.json) records exact versions and hashes. A separate dependency-free run passed 99 tests with the same 2 skips. Tamper regressions reject changed source, split leakage, altered seals, incomplete fitting, hidden pretrained weights, private rows, impossible coverage and contradictory primary/joint effects.

NEXUS quantitative figures were exported as standalone PNG/SVG/PDF and visually inspected as generated files. Tables are regenerated in CI; chemistry, official metric and fresh CPU neural/ranker contracts are exercised separately. Interface changes are checked through source, links and data integrity. All work used CLI, filesystem and authenticated APIs, without PC-screen interaction or Playwright. Historical screenshots below retain their original dates and scope.

## Earlier completed extension — 25 October 2026

## Completed model/control extension

The [retrospective study](../reports/revalidation_review.md) completed on Kaggle on 24 October; source, inputs, public visibility, completion status and aggregate artifact identity were authenticated on 25 October. All 250 eligible molecules were retained across three models, three novelty classes and three spectral conditions. The build quarantined 1,710 mass-conflicting identities and 13 unparseable candidates, verified 1,071,863 candidate fingerprints and rebuilt D2 without cache. No raw source rows or model files were exported.

All four owned Kaggle publications were rechecked on 25 October: public, `COMPLETE`, and matching their tracked source hashes. Dashboard JavaScript passed Node syntax checking and every local HTML asset/document link resolved without browser or screen interaction.

The complete updated local suite passed **297 tests and 101 subtests, with 2 skipped**. Source-bound evidence regressions check the executed code inventory, model/input manifest, notebook source, aggregate artifact hash, unchanged training input relative to the earlier overlap audit, paired-control summaries and retained denominators. CI regenerates both evidence reviews and notebook bundles to detect drift.

This extension used **CLI, filesystem and authenticated APIs only**. No computer-screen interaction, browser automation or Playwright was used. Dashboard changes were verified through source/data integrity checks; no new visual-rendering claim or screenshot is made. Earlier screenshots and their dated verification boundaries are retained below.

The study is retrospective, with historical selection and incomplete original weight-lineage limitations. It is not an untouched external generalization result or a new competition submission.

## Earlier verification record — 23–24 October 2026

## Delivered scope

Research positioning, reviewer documentation, model/data cards, source-linked evidence ledger, downloadable review, NEXUS dashboard evidence filters, owned Kaggle notebook editorial corrections, chart styling and reproducible publication bundles. No new full-pipeline score, paid service or competition submission.

## Executed checks

| Check | Result | Boundary |
|---|---|---|
| Full local Pytest suite | 278 passed, 2 skipped; 84 subtests passed | Synthetic/regression coverage; no licensed-data evaluation |
| New evidence-registry regressions | 3 passed | Reject duplicate IDs, unsafe/missing sources, unsupported states and invalid MRR values |
| Both notebook code-cell syntax checks | Passed | Syntax alone does not establish runtime success |
| Kaggle publication readback | Both COMPLETE, public, source hashes match | [Machine-readable readback](../reports/kaggle_publications.json); no new leaderboard scores |
| Desktop dashboard | Profile and evidence loaded; navigation and rendered status checked | Real browser accessibility and screenshot capture; no mobile/other-engine attestation |

One local skip is the Windows symlink privilege restriction; that containment case remains active on Linux CI. The other is an existing dependency/environment-conditional test. Two existing pandas dtype warnings do not affect test success. RDKit 2026.3.3 was installed into an ignored local test-dependency directory to avoid changing the shared Python installation.

The final browser review and screenshots used accessibility controls without Playwright, as requested. Prior filter/search/empty-result interaction checks were completed before that constraint was added; they were not repeated with Playwright afterward.

## Visual evidence

![NEXUS CASMI desktop research overview](assets/nexus-casmi-research-desktop.png)

![Source-linked evidence explorer](assets/nexus-casmi-evidence-explorer.png)

These are local source-rendered screenshots, not evidence of a publicly deployed web service. The historical public MRR values remain author-log records. Source-file hashes in the exported review use normalized UTF-8/LF text for consistent Windows/Linux regeneration.

CI maintenance also declares pandas/Pytest for integration modules and uses Pytest to collect their function-based regressions. Push checks run on main; pull requests run their own checks to avoid duplicate branch/PR executions.

## Executed data-evidence extension — 24 October 2026

The separate [public-example content audit](PUBLIC_EXAMPLE_OVERLAP.md) completed on Kaggle, with matching notebook source and a verified aggregate artifact hash. All 2,539,608 training rows were scanned; all 1,213 public spectra / 400 molecules matched full numerical content and acquisition mode/adduct. Two coarse candidate pairs were rejected. No hidden evaluation or model-performance conclusion follows.

Twelve initial overlap/matcher/Arrow regressions passed, followed by sixteen focused evidence/dashboard tests. The complete updated suite passed: **287 tests, 2 skipped, 96 subtests** (the same two existing pandas dtype warnings). The new measured-evidence filter was selected with browser accessibility controls and rendered exactly one of six evidence cards. No Playwright was used for this extension.

![Executed overlap evidence in the NEXUS dashboard](assets/nexus-casmi-measured-overlap.png)
