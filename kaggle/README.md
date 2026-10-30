# Kaggle research publications

Nine owned publications have completed authenticated source/status readback. They use NEXUS identity and disabled-internet execution; their scientific scopes differ. None is a new competition submission.

| Publication | Bundle | Evidence |
|---|---|---|
| [Your Ranker Learns How You Built the Pool](https://www.kaggle.com/code/senanuretin/casmi26-your-ranker-learns-how-you-built-the-pool) | `publications/shortcut/` | Toy shortcut demonstration and historical author-log tables |
| [A Bigger Candidate Pool Lowered My Score](https://www.kaggle.com/code/senanuretin/casmi26-a-bigger-candidate-pool-lowered-my-score) | `publications/pool/` | Candidate coverage/decoy analysis; no newly fitted ranker |
| [Exact Public-Example Overlap Audit](https://www.kaggle.com/code/senanuretin/casmi26-exact-public-example-overlap-audit) | `audits/public-overlap/` | Full-content scan; [method](../docs/studies/public-example-overlap.md) |
| [Corrected Models and Negative Controls](https://www.kaggle.com/code/senanuretin/casmi26-corrected-models-and-negative-controls) | `audits/corrected-revalidation/` | Frozen original weights; [results](../reports/revalidation_review.md) |
| [Fresh Training Family Holdout](https://www.kaggle.com/code/senanuretin/casmi26-fresh-training-family-holdout) | `audits/advanced-internal/` | Fresh weights and sealed internal query-family test; [results](../reports/advanced_review.md) |
| [MassBank External Independence](https://www.kaggle.com/code/senanuretin/casmi26-massbank-external-independence) | `audits/external-independence/` | Authenticated negative feasibility: zero novel molecules; no model scoring |
| [GNPS2 External Independence](https://www.kaggle.com/code/senanuretin/casmi26-gnps2-external-independence) | `audits/gnps-independence/` | Authenticated insufficient population: seven novel molecules; no model scoring |
| [GNPS Community External Independence](https://www.kaggle.com/code/senanuretin/casmi26-gnps-community-external-independence) | `audits/gnps-expanded-independence/` | Authenticated original-snapshot insufficiency: 7/936 novel molecules; no model scoring |
| [Purged Family External Validation](https://www.kaggle.com/code/senanuretin/casmi26-purged-family-external-validation) | `audits/purged-external/` | Eight fresh artifacts after full-reference family purge; [all outcomes](../reports/purged_external_review.md) |

## Completed separately registered study

[Purged Family External Validation](https://www.kaggle.com/code/senanuretin/casmi26-purged-family-external-validation) completed on free CPU with disabled internet and no archived-model attachment. Its [registered design](../protocols/purged-external-v1.json) purges selected external families from all training/reference rows before fresh fitting. [Execution history](../reports/purged_external_execution_history.json) preserves both attempts; the [publication authentication](../reports/purged_external_publication.json) verifies completed model evidence. The original novelty gates remain unsuccessful.

`build_purged_external_notebook.py` preserves an existing registered attempt by default. An intentional repaired attempt requires `--new-attempt` and preserved earlier attempt history; a completed publication stays frozen. `check_purged_external_publication.py` authenticates the bounded result and `build_purged_external_report.py` renders verified outcomes.

## Build and authenticate

Run from the repository root using existing local Kaggle authentication:

```bash
python scripts/build_research_publications.py
python scripts/check_kaggle_publications.py --output reports/kaggle_publications.json
python scripts/build_overlap_notebook.py
python scripts/check_overlap_publication.py
python scripts/build_revalidation_notebook.py
python scripts/check_revalidation_publication.py
python scripts/build_advanced_notebook.py
python scripts/check_advanced_publication.py
python scripts/build_external_audit_notebook.py
python scripts/check_external_audit_publication.py
python scripts/build_gnps_audit_notebook.py
python scripts/check_gnps_audit_publication.py
python scripts/build_gnps_expanded_audit_notebook.py
python scripts/check_gnps_expanded_audit_publication.py
python scripts/build_purged_external_notebook.py
python scripts/check_purged_external_publication.py
```

Builders prepare or verify tracked source bundles. Readback commands inspect remote state and download bounded aggregate artifacts; they do not start execution or submit predictions. A notebook URL alone does not attest completion. Source, public visibility, declared attachments and relevant artifact hashes must match.

The two explanatory publications and source/content audits use CPU. The retrospective and fresh internal model studies used free Kaggle GPU execution; the separate purged-reference model study used CPU. Fresh-fitting studies attach no pretrained-model dataset. Exact runtime versions and execution-specific inputs are recorded in each artifact; the local portable environment is separate.

## Completed-source policy

Frozen model-study builders verify their completed embedded source instead of regenerating it from changing main-branch code. Version 1 of the fresh study was superseded before score inspection after a synthetic generated-cache counterexample; version 2 supplies the authenticated results. [Protocol, amendment and limits](../docs/studies/fresh-training-family-holdout.md).

New execution requires an intentionally prepared, separately declared experiment and available free compute/input permissions. Do not push a modified bundle over a completed evidence record or tune against the finished test. Credentials, licensed inputs, private row outputs and weights stay outside Git. Historical auto-submit helpers are not publication-maintenance commands. [Report ownership](../reports/README.md) · [Reproducibility](../docs/guides/reproducibility.md).

## Separately registered follow-on awaiting results

`audits/unopened-coverage/` contains the [new internal coverage design](../protocols/unopened-coverage-v1.json). Its source and attempt history are separate from all completed notebooks. No new performance result is established by preparing this bundle. The study uses free CPU, disabled internet and no pretrained weights, competition submission or deployment.
