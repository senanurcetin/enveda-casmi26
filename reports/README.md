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
| `research_environment_validation.json` | Dated installed-environment/JUnit attestation | Regenerate only from an actually checked environment and completed test run |

Without `--tables-only`, the advanced/coverage builders export PNG/SVG/PDF figures to `docs/assets/` in the locked plotting environment. CI checks deterministic table regeneration; font/platform rendering is not a bitwise image guarantee.

Completed experimental artifacts bind the actual executed notebook source, not the current main branch. Do not edit result JSON to fit new code or update dates to imply rerunning an experiment. Future experiments need distinct declared protocols and evidence records. [Study methods](../docs/README.md#completed-studies) · [Scientific assessment](../docs/research-assessment.md).
