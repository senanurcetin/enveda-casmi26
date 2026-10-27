# Script entry points

Run commands from the repository root in the [documented scientific environment](../docs/guides/reproducibility.md). The Python tools share paths/modules with frozen notebook source inventories; keep their stable entry points when extending the project. Use each CLI's `--help` where provided before an operator-run data operation.

| Task | Entry points | Output or contract |
|---|---|---|
| Validate environment | `check_research_lock.py` | Exact pins/hashes and installed-version attestation |
| Build approved source artifacts | `build_structures.py`, `build_store.py`, `build_pool.py`, `build_fingerprints.py` | New external-to-Git artifact destinations |
| Audit stores and fingerprints | `audit_store_integrity.py`, `audit_pool_fingerprints.py`, `manifest_artifacts.py` | Bounded aggregates and explicit SHA-256 inventories |
| Generic baseline validation | `run_validation.py`, `compare_runs.py`, `compare_control_runs.py` | Signed query/source pairing and scoped comparisons |
| Complete scientific studies | `run_revalidation_study.py`, `run_advanced_study.py` | Fixed model/control or fresh-training contracts |
| External-source preparation and feasibility | `prepare_massbank_external.py`, `prepare_gnps_external.py`, `run_external_audit.py`; `build_external_audit_notebook.py`, `build_gnps_audit_notebook.py` and corresponding readback checkers | Pinned CC0 inputs, full-training molecular/content exclusions and explicit insufficient-population outcomes |
| Gated C3/external follow-on | `build_external_c3_notebook.py`, `run_external_c3_study.py`, `external_study_training.py` | Experimental implementation; requires authenticated sufficient external novelty, fresh fitting and a pre-test seal; no performance result established yet |
| Build/verify publication bundles | `build_research_publications.py`, `build_overlap_notebook.py`, `build_revalidation_notebook.py`, `build_advanced_notebook.py`; corresponding `check_*publication*.py` | [Kaggle operations](../kaggle/README.md) |
| Regenerate reviews and figures | `research_evidence.py`, `build_revalidation_report.py`, `build_advanced_report.py`, `build_coverage_report.py` | [Generated reports](../reports/README.md) |
| Optional source diagnostics | `profile_data.py`, `audit_analog_retention.py`, `diagnose_validation_slices.py` | Descriptive data/reference/slice analysis |
| Optional warehouse reconciliation | `render_bigquery_sql.py`, `reconcile_bigquery_profile.py` | Read-only template rendering and offline reconciliation; no default live scan |

Training, tuning, merge, probe and shell-chain scripts retain historical experiment recipes. Many assume an existing local artifact tree, specific filenames or the original Windows environment. They are not the default reproduction workflow or proof of a completed result. Some legacy submission helpers can publish/submit when deliberately executed; documentation maintenance and current research checks do not invoke them. No shell chain is run as part of setup.

Do not change frozen historical bundles to match current code. Declare a new experiment before modifying populations, model settings or test selection. [Original experiment log](../EXPERIMENTS.md).
