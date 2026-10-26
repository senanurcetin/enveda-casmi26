# Documentation

CASMI26 is a noncommercial molecular-identification research project by Senanur Cetin under the NEXUS identity. Documentation describes the current implementation and distinguishes historical competition results, completed internal experiments and unresolved external-validation questions.

| Purpose | Document |
|---|---|
| Research question, principal results and repository layout | [Project overview](../README.md) |
| Components, data flow and implementation boundaries | [Architecture](../ARCHITECTURE.md) |
| Setup, checks and report regeneration | [Reproducibility](guides/reproducibility.md) |
| Source reconciliation, spectrum stores and fingerprint provenance | [Data integrity](guides/data-integrity.md) |
| Knockouts, controls, pairing, uncertainty and coverage | [Evaluation methodology](guides/evaluation.md) |
| Approved logo, palette and attribution | [NEXUS identity](guides/brand-identity.md) |
| Model attribution and intended use | [Model card](../MODEL_CARD.md) |
| Dataset access, information boundaries and licensing | [Data card](../DATA_CARD.md) |
| Present scientific strengths and remaining limitations | [Research assessment](research-assessment.md) |
| Executed checks and their scope | [Verification](verification.md) |

## Completed studies

| Study | Method | Results |
|---|---|---|
| Exact public-example overlap | [Full-content audit](studies/public-example-overlap.md) | [Aggregate evidence](../reports/public_overlap_audit.json) |
| Frozen-model retrospective controls | [Design and source policy](studies/frozen-model-revalidation.md) | [Results and intervals](../reports/revalidation_review.md) |
| Fresh training on separated query families | [Predeclared design](studies/fresh-training-family-holdout.md) | [Results, ablations and slices](../reports/advanced_review.md) |
| Candidate versus ranking loss | Post hoc accounting on the frozen internal study | [Coverage diagnosis](../reports/candidate_loss_review.md) |

## Source and maintenance policy

[EXPERIMENTS.md](../EXPERIMENTS.md) is the original dated experiment log. Historical entries are evidence records, not current execution instructions. The [evidence registry](../reports/research_evidence.json) links claims to source files; its [generated review](../reports/research_review.md) records their integrity hashes. [Kaggle publication operations](../kaggle/README.md), [script entry points](../scripts/README.md) and [report ownership](../reports/README.md) describe the executable artifacts.

Earlier overlapping audit notes and maintenance diaries were consolidated into these guides. Their original text remains recoverable in Git history. Completed notebook source snapshots, protocol files and experimental JSON remain intact.
