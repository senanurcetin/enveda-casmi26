# Portfolio maintenance verification — 25 October 2026

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
