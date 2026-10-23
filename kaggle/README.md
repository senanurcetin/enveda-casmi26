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
