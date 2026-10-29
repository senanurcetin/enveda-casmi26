# Research evidence and generated reports

This directory contains public aggregates and provenance records. Licensed source records, private molecule identifiers, per-query rankings and trained weights remain outside Git. Historical observations, executed data-backed results and local synthetic verification have distinct scopes.

| Files | Role | Maintenance |
|---|---|---|
| `data_profile.json` | Archived training/source profile | Preserve its historical scope |
| `research_evidence.json` | Claim registry linking evidence states and source paths | Update deliberately when claims or document locations change |
| `research_review.md` | Source-hashed view of that registry | Generate with `python scripts/research_evidence.py` |
| `kaggle_publications.json` | Readback for two explanatory historical publications | Authenticated `check_kaggle_publications.py`; execution does not regenerate historical score tables |
| `public_overlap_audit.json`, `public_overlap_publication.json` | Exact-content scan and publication provenance | Preserve executed artifact; authenticate with `check_overlap_publication.py` |
| `revalidation_study.json`, `revalidation_publication.json`, `revalidation_review.md` | Frozen-model retrospective evidence and review | Authenticate with `check_revalidation_publication.py`; render with `build_revalidation_report.py` |
| `advanced_study.json`, `advanced_publication.json`, `advanced_execution_history.json` | Fresh internal study, source readback and outcome-blind attempt history | Preserve source/weight/protocol bindings; authenticate with `check_advanced_publication.py` |
| `advanced_review.md` | All fresh-study outcomes, controls and slices | `python scripts/build_advanced_report.py --tables-only` |
| `candidate_loss_budget.json`, `candidate_loss_review.md` | Post hoc aggregate-only loss accounting | `python scripts/build_coverage_report.py --tables-only` |
| `external_source_preparation.json`, `external_audit.json`, `external_audit_publication.json` | Pinned MassBank CC0 input preparation and authenticated negative source-independence result | Preserve completed source; `check_external_audit_publication.py` authenticates remote evidence |
| `gnps_external_source_preparation.json`, `gnps_external_audit.json`, `gnps_external_audit_publication.json` | Two-library GNPS source preparation and authenticated insufficient external population | `check_gnps_audit_publication.py`; this is not model-performance evidence |
| `gnps_expanded_source_preparation.json`, `gnps_expanded_external_audit.json`, `gnps_expanded_external_audit_publication.json` | Expanded GNPS original-snapshot feasibility; 7/936 novel molecules | Preserve negative conclusion and completed source; `check_gnps_expanded_audit_publication.py` |
| `purged_external_execution_history.json` | Separately registered family-purged model-study attempts | Preserved pre-fit failure and authenticated completed execution; performance remains in the separately authenticated aggregate |
| `research_environment_validation.json` | Dated installed-environment/JUnit attestation | Regenerate only from an actually checked environment and completed test run |

Without `--tables-only`, the advanced/coverage builders export PNG/SVG/PDF figures to `docs/assets/` in the locked plotting environment. CI checks deterministic table regeneration; font/platform rendering is not a bitwise image guarantee.

Completed experimental artifacts bind the actual executed notebook source, not the current main branch. Do not edit result JSON to fit new code or update dates to imply rerunning an experiment. Future experiments need distinct declared protocols and evidence records. [Study methods](../docs/README.md#completed-studies) · [Scientific assessment](../docs/research-assessment.md).

## Completed purged-reference experiment

- [Authenticated aggregate](purged_external_study.json) and [publication readback](purged_external_publication.json).
- [Complete result tables, all 11 contrasts and NEXUS figures](purged_external_review.md).
- [Both execution attempts](purged_external_execution_history.json), including the preserved pre-fit failure.

The authenticated family-purged study removed 261,826 of 2,539,608 training spectra and fitted eight fresh learned artifacts. On 150 internal test molecules, expanded C3 truth coverage rose from 13 to 17 (paired change +0.0267, 95% CI [+0.0067, +0.0533]). On 150 external molecules, expanded C2/C3 MRR@25 was 0.4454/0.0466; expanded-minus-legacy ranking effects remain unresolved in the joint 11-contrast intervals. External C3 truth coverage remains only 23/150. Original-snapshot novelty gates remain unsuccessful; this is generalization relative to a rebuilt purged corpus, not independent laboratory replication.
