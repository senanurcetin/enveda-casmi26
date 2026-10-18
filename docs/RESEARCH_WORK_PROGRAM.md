# Enveda CASMI26 — Research-only work program

**Owner:** NEXUS open-source research portfolio  
**Working baseline:** `main` as of 2026-10-18; previous research frozen 2026-10-08.  
**Scope:** open scientific analysis, reliable experiments, and auditable Kaggle/BigQuery data engineering.  
**Out of scope:** commercialization, SaaS, deployment, paid services, monetization, customer accounts, authentication, API products, sales pages, and leaderboard claims without new evidence.

This is an execution plan **committed to the repository**, not a promise that a process is running in the background. Later stages require reproducible datasets, appropriate permissions and normal execution in a development environment. Never place private datasets, training derivatives, credentials or private leaderboard answers into the public repository.

## Prioritization and schedule

| Period | Priority | Work package | Acceptance criteria | Status on 2026-10-18 |
| --- | --- | --- | --- | --- |
| Sprint 0 · 18 Oct | P0 | Run comparability and provenance | Class C1/C2/C3 rank arrays paired by unique query ID; identical query spectra/exclusions; reordered run lists align; mismatched or legacy unverified runs cannot produce ACCEPT/REJECT | **Implemented** in validation code and offline tests |
| Sprint 0 · 18 Oct | P0 | Experiment file reproducibility | Explicit query seed, library and per-molecule spectra cap; report records query IDs/signatures; collision-resistant output filename, no silent overwrite | **Implemented** |
| Sprint 0 · 18 Oct | P0 | BigQuery-to-local count reconciliation | Compare each source library's spectrum and distinct-2D-structure counts, flag unknown/missing sources and total drift, record file hashes | **Implemented** as offline CLI; **live warehouse not queried** |
| Sprint 0 · 18 Oct | P0 | Automated QA | Offline unit tests and Python syntax checks on feature branch; merge only after green GitHub Actions | **Implemented** |
| Sprint 1 · 19–21 Oct | P1 | Authoritative data inventory | Record actual BigQuery project/dataset/field names, permissions, snapshot time, scan-cost dry run, Kaggle dataset/input version | **Pending evidence/access** |
| Sprint 1 · 19–21 Oct | P1 | Leakage and quality diagnostics | Verify public example overlap by exact spectra, not approximate fingerprints; validate structure-disjoint splitting and cross-library contamination | **Exact-comparison code implemented; full data rerun pending** |
| Sprint 2 · 22–25 Oct | P1 | Evaluation stress tests | Preserve molecule-level pairing; add blank-spectrum / shuffled-spectra / pool-membership-only negative controls; score all novelty classes and source/adduct slices | **Requires data and compute** |
| Sprint 2 · 22–25 Oct | P1 | Reproducible comparisons | Fixed seeds, scenario manifests, paired confidence intervals, explicit NOT PROVEN where uncertainty spans zero | **Partially enabled** |
| Sprint 3 · 26–29 Oct | P2 | One-at-a-time algorithm experiments | Compare one hypothesis against the frozen best baseline; regression tests, source manifest and independent validation before considering a submission | **Not started** |
| Sprint 3 · 26–29 Oct | P2 | Open research write-up | Document negative/positive outcomes, uncertainty, Kaggle rules, pipeline commands and reproduction costs | **Not started** |

Dates are **planned target windows**, not scheduled jobs or future autonomous execution.

## Research sprint update — statistical validity and scoring (2026-10-18)

- **Implemented:** A fixed-grid, approximately simultaneous 95% bootstrap
  interval for C1/C2/C3 mixtures, preserving per-molecule pairing across
  classes and guarding against treating 12 pointwise confidence intervals as
  a joint 95% statement.
- **Implemented:** Rank-slot-preserving SMILES parsing in the offline scorer.
  Invalid and blank guesses still occupy their original rank positions, so a
  malformed list cannot receive an artificially improved local MRR@25.
- **Implemented:** Comparison guard for incompatible query-signature algorithm
  versions; unverified historical comparisons retain an explicit unverified
  verdict, including class-specific numeric summaries.
- **Verification:** Synthetic regression tests and GitHub Actions checks
  exercise numerical behaviors; a complete competition-data rerun and live
  BigQuery reconciliation have **not** been performed.
- **Unchanged:** Models, competition submissions, portfolio score of 0.422
  dated 2026-10-08, and non-commercial scope. These changes improve research
  validity, **not** predictive accuracy.

## Scientific quality gates

### Gate G0 — population identity (**enforced in code**)
- Both reports contain exactly the same set of nonempty, unique `query_ids`; order may differ.
- `query_signatures` hash actual source spectra, precursor/adduct/energy and held-out reference row indices.
- C1/C2/C3 rank arrays each contain exactly one valid rank (0–25) per molecule.
- Any reported MRR must equal the MRR recomputed from its ranks.
- Legacy run files without query provenance are **rejected by default**; a dedicated override only produces `UNVERIFIED`, never `ACCEPT` or `REJECT`.

### Gate G1 — source integrity (**tools ready; live validation pending**)
- Source of truth: approved competition inputs for model data, BigQuery aggregates for warehouse quality, `reports/data_profile.json` for the **historical** snapshot only.
- Do not aggregate per-library distinct structure counts into a global distinct count.
- Validate BigQuery source schema, location, timestamp, costs and scan bytes before executing SQL.
- Reconcile BigQuery-derived counts and provenance hashes. Investigate all deviations before use.

### Gate G2 — independent model assessment (**requires data**)
- No held-out 2D structure may leak from any training library to the validation reference channels.
- Report sensitivity of novelty classes C1/C2/C3, adduct, ionization mode and spectral instrument.
- Run negative controls (blank/shuffled spectra and membership-only priors). A suspiciously high negative-control score is a failure.
- Use query-paired bootstrap differences; leaderboard improvements under the model/noise floor do not justify confident claims.

### Gate G3 — release-to-public-research (**no productization**)
- CI green; every new run has manifest, seed, experiment hypothesis and decision log.
- Existing `0.422` public MRR@25 is **8 October 2026 history** and must not be represented as a new result.
- No paid deployment, user tracking, API/product layer or new production app. The existing NEXUS panel remains a static research visualization only.

## Reproduce the newly implemented gates

Generate future signed validation reports on the **same source snapshot**:

```sh
python scripts/run_validation.py --analog --tag baseline \
  --library enveda-np-examples --max-queries 250 \
  --max-spectra 5 --query-seed 0

python scripts/compare_runs.py runs/BASELINE.json runs/EXPERIMENT.json
```

A historical run with no query IDs or spectra fingerprints can be explored with
`--allow-unverified-legacy-order`, but this **cannot** generate a valid scientific acceptance decision. Never splice ranks from different snapshots or split protocols.

Render and dry-run the existing read-only BigQuery query, export the aggregate CSV, then reconcile it:

```sh
python scripts/render_bigquery_sql.py --project your-project-id \
  --dataset casmi26 --table train_spectra --output runs/library_quality.sql
bq query --use_legacy_sql=false --dry_run < runs/library_quality.sql
# Only execute the read-only query after schema, permission and cost review.
python scripts/reconcile_bigquery_profile.py \
  --profile reports/data_profile.json \
  --aggregate runs/library_quality.csv \
  --json-out runs/bigquery_reconciliation.json
```

Checks:

```sh
python -m unittest discover -s tests -p '*_stdlib.py' -v
python -m py_compile scripts/run_validation.py scripts/compare_runs.py \
  scripts/reconcile_bigquery_profile.py src/casmi/validate/run_alignment.py
```

## Decision policy

Progress is evaluated from tested commits and evidence artifacts, **not** from a commitment to keep Kaggle/BigQuery workloads running. Do not trigger expensive training, new submissions, reprocessing, or public claims simply to complete a calendar milestone. The existing inference pipeline stays frozen unless a controlled independent experiment is demonstrably better.
