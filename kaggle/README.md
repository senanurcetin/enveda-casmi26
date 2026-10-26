# CASMI26 Kaggle research publications

Two owned, public notebooks explain the original experiments:

- [Your Ranker Learns How You Built the Pool](https://www.kaggle.com/code/senanuretin/casmi26-your-ranker-learns-how-you-built-the-pool): simplified simulation/feature shortcut experiment; competition training input required.
- [A Bigger Candidate Pool Lowered My Score](https://www.kaggle.com/code/senanuretin/casmi26-a-bigger-candidate-pool-lowered-my-score): candidate coverage versus decoy density; competition data and the existing owned structure-source dataset required.

Both use NEXUS chart styling and an explicit research-status notice. Hard-coded tables are historical author logs; executing them does not regenerate those historical leaderboard results. The shortcut demonstration is a standard-InChIKey toy experiment, not the full corrected tautomer-aware pipeline. The pool analysis is a coverage study, not a newly fitted ranker. No notebook in this maintenance workflow submits predictions.

## Build and verify

```bash
python -m pip install nbformat
python scripts/build_research_publications.py
python scripts/check_kaggle_publications.py --output reports/kaggle_publications.json
```

`build_research_publications.py` generates deterministic publication bundles under `kaggle/publications/` from the tracked builders. The metadata preserves existing competition/dataset attachments and runs CPU-only with internet disabled. Publication uses `kaggle kernels push -p kaggle/publications/shortcut` and the equivalent `pool` directory. These are notebook runs, not leaderboard submissions. Do not run the historical auto-submit shell scripts for publication maintenance.

The readback compares ordered cell-source hashes and checks remote run status and public visibility. A matching source plus `COMPLETE` establishes that Kaggle completed that notebook version; it does not prove the old score tables or full-pipeline correctness. Failed and pending runs remain visible in the report.

Authentication uses the already configured Kaggle CLI. Tokens and local credential files are never copied into the repository. Raw data, notebook output datasets and generated structure records are not committed. Source-release availability and account CPU quotas can still block execution.

## Executed full-content data audit

[Exact Public-Example Overlap Audit](https://www.kaggle.com/code/senanuretin/casmi26-exact-public-example-overlap-audit) is a separate third research notebook. It scanned all training rows and confirmed 1,213/1,213 public spectra and 400/400 molecules as exact numeric-content matches, including adduct/mode consistency, while rejecting two coarse candidate pairs. It does not train or score a model. [Protocol, input hashes and limits](../docs/PUBLIC_EXAMPLE_OVERLAP.md).

Its source bundle lives under `kaggle/audits/public-overlap`, outside the two historical-publication bundles. Generate it with `python scripts/build_overlap_notebook.py` and verify its completed source, input attachments and aggregate JSON with `python scripts/check_overlap_publication.py`. The existing two-notebook readback remains scoped to those two publications.

## Corrected-model retrospective study

[Corrected Models and Negative Controls](https://www.kaggle.com/code/senanuretin/casmi26-corrected-models-and-negative-controls) rebuilds full-source artifacts and compares D1, D1+D2 and the frozen five-seed v16 ranker under C1/C2/C3 and paired spectral perturbations. Its first build exposed 1,710 mass-conflicting identities. The revised protocol quarantines those identities explicitly and keeps every eligible query in the denominator using independent training labels.

The bundle lives under `kaggle/audits/corrected-revalidation`. Build with `python scripts/build_revalidation_notebook.py`, push that directory and read back only a completed matching run with `python scripts/check_revalidation_publication.py`. A notebook URL alone does not attest completion. [Full protocol and interpretation limits](../docs/CORRECTED_MODEL_REVALIDATION.md). This study uses attached original owned model files, not an external solution, and submits no predictions.

**Completed and authenticated:** execution 24 October 2026, readback 25 October 2026. [Result tables and adverse controls](../reports/revalidation_review.md) · [Source/artifact readback](../reports/revalidation_publication.json). All four owned publications are complete; their distinct evidence scopes remain explicit. The model study uses free GPU execution with disabled internet and exports aggregates only.

## Fresh training and sealed internal-family benchmark

The fifth owned publication, [Fresh Training Family Holdout](https://www.kaggle.com/code/senanuretin/casmi26-fresh-training-family-holdout), completed and passed exact-source/artifact readback on 26 October 2026. It attaches structure sources and the official metric wheel, with **no pretrained models**, GPU enabled and internet disabled. Two new fingerprint models and eight three-seed rankers are sealed before test features. The 600/150/250 supervised query partitions share no selected family/metric identities.

[All outcomes, controls, slices and figures](../reports/advanced_review.md) · [Readback](../reports/advanced_publication.json) · [Protocol and amendment](../docs/ADVANCED_INTERNAL_STUDY.md). Version 1 was superseded outcome-blind after a synthetic generated-cache counterexample; version 2 supplies the authenticated results. This is internal validation with shared candidate/reference knowledge, not external spectra or competition scoring.

The executed bundle is frozen under `kaggle/audits/advanced-internal`. `python scripts/build_advanced_notebook.py` checks its identity rather than silently regenerating it from evolving source. Use `python scripts/check_advanced_publication.py` for live authenticated readback and `python scripts/build_advanced_report.py` for the report and standalone PNG/SVG/PDF figures. Future experiments require a separate declared protocol and publication; do not tune against this completed test.
