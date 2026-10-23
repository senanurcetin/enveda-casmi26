# NEXUS CASMI26 — Research Intelligence

**Publication review, 24 October 2026:** both owned Kaggle research notebooks are COMPLETE with matching source hashes. See [publication readback](../reports/kaggle_publications.json), [current documentation](README.md) and [model card](../MODEL_CARD.md). The dated host check below describes a different earlier environment; it is not the current publication status. Full corrected-pipeline evaluation remains pending.

**Status:** historical research dashboard + reproducible QA utilities.  
**Snapshot date:** 2026-10-08. This extension does **not** claim a new Kaggle run, a new model score, or a new live BigQuery connection.

## Deliverables

| Asset | Purpose | Source of truth |
| --- | --- | --- |
| `site/index.html` | Responsive NEXUS 2.0 research dashboard | `reports/data_profile.json`, README, experiment and handoff logs |
| `site/assets/nexus-official.svg` | Approved transparent NEXUS identity artwork | Google Drive: Logo / NEXUS_Identity_2_0 |
| `scripts/audit_submission.py` | Exact test-ID coverage and ranked-SMILES preflight | Submitted CSV and official test IDs |
| `sql/bigquery/library_quality.sql` | Read-only source-library and data-quality aggregate | BigQuery table populated from official train.parquet |
| `scripts/render_bigquery_sql.py` | Validate and insert BigQuery identifiers safely | Operator-supplied GCP project/dataset/table |
| `tests/test_*_stdlib.py` | Offline reproducibility and regression smoke tests | Checked-in code, historical JSON |
| `.github/workflows/nexus-research-checks.yml` | CI for offline tests | Repo source |

## Updated research interface and Kaggle readiness (20 October 2026)

The responsive NEXUS 2.0 dashboard includes a dedicated **Research readiness**
section separating completed repository integrity safeguards from pending
Kaggle authentication, licensed input access, independent evaluation and
Notebook submission. The panel is intentionally read-only; it never invents
live scores or marks unperformed runs as completed. The existing official
NEXUS SVG identity and historical experiment graphics remain intact.

A new local, non-destructive readiness gate checks official file presence,
the submission CSV header, required source scripts and (if installed)
Parquet footers/columns without loading the training dataset into memory:

```bash
python scripts/kaggle_submission_readiness.py \
  --data-dir /path/to/authorized/competition-data \
  --json-out runs/kaggle-input-readiness.json
```

Expected official files: `train.parquet`, `test.parquet`,
`sample_submission.csv`. A PASS is only local file/schema readiness,
not proof of competition-rule acceptance, GPU/runtime viability,
successful notebook execution or leaderboard evaluation.

**SentinelX host check on 20 October 2026:** Windows agent reachable;
on its current execution environment, Kaggle CLI and Python commands
were not on PATH, `KAGGLE_API_TOKEN` was not set, and the expected
CASMI dataset directory was absent. No private credential was printed
or persisted. This check describes the connected host environment
rather than the user's Kaggle account status.

## View the NEXUS dashboard

From the **repository root**:

```sh
python -m http.server 8000
```

Open http://localhost:8000/site/ .

The site is self-contained: no tracking scripts, third-party CSS/fonts, API keys, BigQuery credentials or Kaggle private data. It dynamically reads the **historical** `reports/data_profile.json` from the same repository; the stored leaderboard numbers come from the documented README/HANDOFF experiment snapshot. Opening `site/index.html` directly as a `file://` URL can prevent the browser from loading the profile. Serve over HTTP instead.

Design references: NEXUS 2.0 official transparent SVG and midnight palette `#0B0F19 / #061524 / #00E0FF / #4EA8FF / #8F9BFF / #6E55C6 / #B9C5CE`. Branding belongs to NEXUS; the competition itself belongs to its organizers.

## Audit a Kaggle submission

The competition evaluates **one ranked list per molecule**, capped at **25 SMILES**. Run the structural audit before executing a submission:

```sh
python scripts/audit_submission.py runs/submission.csv \
  --test-parquet data/test.parquet \
  --json-out runs/submission_audit.json
```

If you have a newline-delimited list of expected molecule IDs, use `--expected-ids ids.txt` instead. That path needs only Python's standard library.

For molecule-level 2D equivalence checks, use the pinned RDKit environment:

```sh
python scripts/audit_submission.py runs/submission.csv \
  --test-parquet data/test.parquet --chemistry \
  --fail-on-warnings --json-out runs/submission_audit.json
```

Checks include:
- Exactly `molecule_id,smiles` columns in that order, no duplicate IDs, no missing/extra test molecules (when IDs supplied).
- Nonempty ranked slots and at most 25 candidates.
- Duplicate literal candidates by default.
- Optional invalid/unparseable SMILES and duplicate canonical InChIKey14 structures (including tautomer equivalence).

The `--chemistry` switch uses the repository's `canonical_ik14` implementation, which depends on the official RDKit version. An audit result with `PASS` does **not** mean the SMILES are correct or the MRR@25 is high. Duplicate predictions use ranked positions in the official metric, so duplicates are **warnings** rather than silently dropped. Structural checks are stricter than some scorer behaviors by design, and never mutate the submission CSV.

## BigQuery: read-only data-quality profile

The existing project is described as Kaggle/BigQuery-connected, but this change cannot validate actual connected table names or issue live queries. The query is therefore a documented template, **not** a deployed BigQuery view.

It assumes that a BigQuery table such as `my-casmi-project.casmi26.train_spectra` contains the original training columns `ingest_lib, inchikey14, precursor_mz, precursor_error_ppm, adduct, ionization_mode, ms2_mzs, ms2_normalized_intensities`. Arrays must be loaded as repeated numeric fields.

Render the table identifier safely:

```sh
python scripts/render_bigquery_sql.py \
  --project my-casmi-project --dataset casmi26 --table train_spectra \
  --output runs/library_quality.sql
```

Review estimated scan cost before running:

```sh
bq query --use_legacy_sql=false --dry_run < runs/library_quality.sql
bq query --use_legacy_sql=false --format=prettyjson < runs/library_quality.sql
```

The query returns per-library counts and issues: invalid molecular keys, invalid precursor masses, missing/mismatched fragment arrays, empty spectra, >10 ppm precursor mismatch, positive ion share, and test-adduct coverage. Distinct structure counts **must not** be summed across libraries, because the same structure can appear in more than one source. No raw training spectra or private labels are exported to the dashboard.

For production BigQuery integration, first verify dataset location, access, source schema, partitioning, storage/billing permissions and Kaggle redistribution rules. Then build a governed aggregate export and give it a source timestamp, query hash and snapshot identifier. Never embed service-account JSON, Kaggle tokens, or raw leaderboard submissions in the site.

## Scientific guardrails / important results

1. **The public example test file is unusable as a true holdout.** The checked-in historical profile recorded 1,213 / 1,213 matches via rounded precursor/peak-count/first-peak buckets. Such coarse matches are potential, not proven, duplicates. The updated profiler now checks full m/z and intensity arrays; rerun it before asserting exact overlap. In any event, do not use known training-overlap examples as a trustworthy holdout.
2. **MRR@25 is evaluated by canonical 2D structure equivalence**, not SMILES-string equality or exact stereochemistry.
3. **The hidden test has novelty classes** (library-known, known-without-spectra, novel). Validation must preserve structure-level holdouts and consider candidate-pool membership shortcuts.
4. **Instrument shift matters.** Training data span heterogeneous libraries/instruments while hidden test spectra come from timsTOF. Analyze source/adduct/energy slices rather than only total mean rank.
5. **Leaderboard increments can be noise.** The recorded 0.422 vs 0.420 fusion difference is historical; don't claim it is statistically significant without an appropriate independent test.
6. **Never retroactively rewrite experiment history.** Compare independent query-level paired bootstrap confidence intervals and log rejected changes too.

## Research-only continuation (2026-10-18)

This repository is an open scientific research/portfolio project, **not** a commercial product. See [the committed research work program](RESEARCH_WORK_PROGRAM.md) for dated priorities, scientific gates, completed code changes and evidence-dependent stages. Future validation runs now carry query identity, actual spectrum-input signatures and reproducible query sampling; the paired comparison tool refuses unverified legacy pairing by default. The new metadata-only BigQuery reconciliation CLI can test source-library counts after a genuine read-only warehouse export, without inventing connection state.

## Next science sprint (evidence gate before retraining)

- Confirm BigQuery/Kaggle credential **access** without revealing keys; inventory authoritative source tables and training-to-evaluation lineage.
- Validate `reports/data_profile.json` against a read-only BigQuery aggregate, including exact row count and source distribution; write the date and input hash into a provenance record.
- Compare C1/C2/C3 holdout validity; run negative controls with blank spectra, shuffled candidate features and candidate-pool-only priors.
- Quantify uncertainty with molecule-level paired bootstrap; prioritize **transfer performance**, not a small public leaderboard increment.
- Only after the previous gates pass, test one causal change at a time (pool membership controls, timsTOF-aware calibration, class-specific rank diagnostics) and record both gains and failures.

## Local quality checks

```sh
python -m unittest discover -s tests -p '*_stdlib.py' -v
python -m py_compile scripts/audit_submission.py scripts/render_bigquery_sql.py
```

These focused checks require **no Kaggle download, RDKit, torch, BigQuery client or cloud credentials**. Full scientific pipeline tests remain in `tests/test_foundations.py` and `tests/test_spectra.py` with the repository dependencies.

## Reproducibility notice

The original project was closed on 2026-10-08. This branch extends its **observability, auditability and brand presentation**, not its recorded experimental performance. Existing Kaggle submission settings, models and scoring functions are unchanged.
