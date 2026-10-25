# Enveda CASMI26 — Research-only work program

**Owner:** NEXUS open-source research portfolio  
**Working baseline:** `main` as of 2026-10-18; previous research frozen 2026-10-08.  
**Scope:** open scientific analysis, reliable experiments, and auditable Kaggle/BigQuery data engineering.  
**Out of scope:** commercialization, SaaS, deployment, paid services, monetization, customer accounts, authentication, API products, sales pages, and leaderboard claims without new evidence.

This is an execution plan **committed to the repository**, not a promise that a process is running in the background. Later stages require reproducible datasets, appropriate permissions and normal execution in a development environment. Never place private datasets, training derivatives, credentials or private leaderboard answers into the public repository.

## Completed data-evidence gate — 24 October 2026

The [full-content public-example audit](PUBLIC_EXAMPLE_OVERLAP.md) completed on Kaggle: all 2,539,608 training rows scanned; all 1,213 public spectra / 400 molecules matched exact arrays and acquisition mode/adduct. Two coarse candidate pairs were rejected. Input/code hashes and authenticated publication readback are committed. This closes the public-example exact-overlap task for that input snapshot.

## Completed model/control gate — readback 25 October 2026

The [corrected-protocol retrospective study](../reports/revalidation_review.md) completed on 24 October; authenticated source/status/artifact verification completed on 25 October. All 250 eligible V-in molecules were retained across D1, D1+D2 and frozen v16 under C1/C2/C3, blank fragments and intensity shuffling. Full-source rebuilds exposed and explicitly quarantined 1,710 mass-conflicting identities and 13 unparseable candidates; 1,071,863 fingerprints and a cache-free D2 index were verified. One absent standard truth key remained in the denominator.

v16 C1/C2/C3 MRR is 0.8324/0.6241/0.2374. Blank fragments markedly reduce v16 performance but improve simple-model C2 ranking; both favorable and adverse diagnostics are published with paired uncertainty. Input/model/code hashes and source-matched public readback are committed. No retraining, competition submission, model-promotion decision or untouched-holdout claim follows.

The current implementation and retrospective execution gates are complete. Independently attested new training, untouched external/scaffold validation, broader source/error slices and optional live BigQuery reconciliation are separate research extensions. Historical progress notes below retain their original dates and are superseded by this completed-run evidence where applicable.

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

## Additional P0 research correction — 2026-10-18

The scientific audit found that C1 queries capped with `max_spectra` previously excluded only the sampled rows, leaving unselected same-structure library references available. `build_queries` now preserves all same-key source-library reference row IDs independently of selected query spectra, C1 also hides metric-equivalent tautomer variants, and the `vin` feature-table builder applies the same exclusion rule. The query signature format was versioned to v2. The synthetic CI regression suite verifies this behavior; historical scores **must not** be recast as corrected without a full data-dependent rerun. See [scientific assessment](CASMI26_RESEARCH_ASSESSMENT.md).

## P0 correction — D2 analog shortlist (2026-10-18)

A follow-up code audit identified that D2 discarded knocked-out reference rows *after* approximate top-k, potentially hiding eligible spectral analogs. A tested eligible-first prefilter now applies knockout exclusions, mass and adduct constraints before selecting candidates. This is an **evaluation selection-bias correction**, not a claim of a new model score. Remaining risk: reference subsampling during offline index construction can still favor the withheld source library; investigate by measuring cross-library index retention on the permitted local training artifacts. See [D2 scientific audit](D2_PREFILTER_KNOCKOUT_AUDIT.md).

## P0 research extension — D2 cross-source index coverage (2026-10-18)

- **Implemented:** opt-in `source_balanced` D2 index selection, selecting a high-quality reference plus qualified cross-library spectra before filling remaining per-structure/adduct slots; the historical `legacy` default remains unchanged.
- **Implemented:** independent, provenance-checked cache for the alternative index; reference selection shared with the metadata-retention audit; same-population experimental CLI flags and synthetic regression tests.
- **Pending:** run the metadata audit against actual licensed CASMI store data, rerun model validation on one fixed query set, then determine whether any MRR improvement is statistically supported. No new competition-data experiments have been run for this change.
- See [predeclared D2 protocol](D2_SOURCE_BALANCED_EXPERIMENT.md). Do not mix cached index variants, claim new Kaggle results, or commercialize the scientific work.

## P1 spectral negative controls — 2026-10-18

- **Implemented:** `blank` and `intensity_shuffle` diagnostics for selected query spectra, preserving source metadata, molecule sampling and novelty-class reference knockout.
- **Implemented:** deterministic per-molecule seed, original-versus-perturbed spectral signatures in run reports, and a fail-closed gate that prohibits scientific model-improvement ACCEPT/REJECT judgments for control runs.
- **CI:** synthetic controls and unverified-verdict tests, not real MS/MS data execution.
- **Pending:** run controls on approved real datasets, audit residual C1/C2/C3 scores and investigate library/pool membership bias. See [negative-control protocol](SPECTRAL_NEGATIVE_CONTROLS.md).

## P1 extension — verified spectral-control diagnostics (2026-10-19)

- **Implemented and covered by synthetic tests:** `compare_control_runs.py` verifies that negative-control and unperturbed reports share exactly the same original query spectra, molecule IDs, exclusion scope, model configuration and sampling parameters. The transformed query fingerprint is allowed to differ by design, and at least one actual perturbation is required.
- **Implemented:** descriptive molecule-paired per-class MRR contrasts and a fixed-grid approximately simultaneous bootstrap band for the diagnostic effect. Output is always `DIAGNOSTIC_ONLY`, with **no model ACCEPT/REJECT** decision. Historical `compare_runs.py` continues to reject altered-spectrum reports.
- **Not performed:** full licensed data rerun, Kaggle/BigQuery live verification or claims of a new predictive score. These remain required evidence gates.

## P0 data-source provenance gate (2026-10-19)

- **Implemented:** explicit streaming SHA-256 inventories for user-authorized local spectrum metadata, peak memmaps, candidate pools and fingerprint artifacts. Manifests contain only relative paths, file lengths and hashes.
- **Implemented:** `run_validation.py --source-manifest` aborts before scoring if any inventoried input differs; the report records the composite inventory hash, and diagnostic control pairings require equal hashes when supplied.
- **Validated with synthetic files:** same-byte integrity, same-size modification detection, path traversal/symlink protection and malformed manifests.
- **Not completed:** verify actual Kaggle dataset lineage, BigQuery schema, source snapshot, or rerun models against real competition artifacts. See [local snapshot instructions](INPUT_SNAPSHOT_MANIFEST.md).

## P0 — verified source snapshots for scientific decisions (2026-10-19)

- **Implemented:** optional local input snapshots now require all static files consumed by validation: four spectrum-store artifacts and the active candidate-pool table, plus both D2 fingerprint artifacts when analog scoring is enabled.
- **Implemented:** `compare_runs.py` distinguishes query-ID/spectrum pairing from full local input provenance. A report without a matching, coverage-checked SHA-256 snapshot is **exploratory only**, even when the bootstrap effect excludes zero.
- **Regression evidence:** the synthetic end-to-end comparison no longer produces `ACCEPT` from identical query fingerprints alone; altered or missing inventory SHA-256 values are rejected.
- **Remaining limits:** approved Kaggle/BigQuery files have not been independently retrieved; local inventory cannot authenticate upstream sources, changes in code dependencies, learned model weights or mutable caches. See [snapshot protocol](INPUT_SNAPSHOT_MANIFEST.md).

## P1 descriptive error slices — 2026-10-19

- **Implemented:** `run_validation.py` persists per-query, label-free metadata bins for adduct, ionization mode, query spectra count and median precursor m/z from original query inputs.
- **Implemented:** `diagnose_validation_slices.py` validates query identities/spectra and metadata correspondence before computing per-class descriptive MRR differences; cells with fewer than 15 molecules are suppressed by default.
- **Prohibited:** interpreting slices as independent significant effects, issuing model ACCEPT/REJECT, or presenting synthetic fixtures as real CASMI results. See [slice protocol](ERROR_SLICE_DIAGNOSTICS.md).
- **Pending:** actual data-dependent slice profiles, instrument-family quality review and independently replicated hypotheses.

## P0 — structural A1 peak-store safeguards (2026-10-19)

- **Implemented:** store loading and post-build verification reject invalid offset counts/order, mismatched binary m/z/intensity arrays, truncated/unaligned files and inconsistent per-row `n_peaks`.
- **Implemented:** an optional read-only `audit_store_integrity.py` CLI and strict row index bounds checks.
- **Test evidence:** numerical synthetic corrupt-file fixtures plus mocked metadata constructor tests in CI. The full real CASMI spectrum store was **not** rerun or checked.
- **Next evidence gate:** run [A1 structural audit](SPECTRUM_STORE_INTEGRITY.md), then the full SHA-256 source manifest and dataset-backed C1/C2/C3 benchmarks on the authorized scientific artifacts.

## P0 — finite MS/MS similarity validation (2026-10-19)

- **Implemented:** entropy and cosine scoring reject mismatched arrays and return zero (not NaN) for invalid, non-finite and zero-intensity MS/MS inputs; valid spectra retain the established scoring path.
- **CI:** conventional positive-intensity invariants, ordinary peak cleaning and new numerical edge cases are exercised together in the pinned Numba scientific environment.
- **Pending:** count/quantify bad-spectrum prevalence on actual authorized datasets; diagnose source/library-specific failures and rerun corrected validation. See [numeric safety audit](SPECTRAL_SIMILARITY_NUMERIC_SAFETY.md).

## P0 — optional full numeric peak-content integrity (2026-10-19)

- **Implemented:** optional bounded-memory peak scan for invalid/non-finite m/z or intensity, reporting aggregate counts with a nonzero exit status when a requested scan finds bad values.
- **Implemented:** NumPy and Numba peak cleaning now reject `m/z <= 0` to preserve physically plausible fragment bins and protect D2 sparse vectors.
- **No automated CI requested:** changed files and synthetic tests were committed with GitHub `[skip ci]` directives; existing research workflow remains available for future deliberate verification.
- **Data gap:** the actual large CASMI store and Kaggle/BigQuery pipelines have not been rerun. See [peak-content research audit](PEAK_CONTENT_QUALITY_AUDIT.md).

## P0 — D1 reference ionization-mode alignment (2026-10-19)

- **Implemented:** D1 groups mixed-mode query acquisitions and evaluates each only against the corresponding ionization-mode reference rows.
- **Tests:** synthetic mixed-mode, row exclusion, reference cap and missing-mode regressions; isolated local assertions passed without invoking CI.
- **Pending:** quantify mixed-mode query prevalence and any rank effects on independently verified, molecule-disjoint CASMI research datasets. See [D1 ion-mode audit](D1_IONIZATION_MODE_AUDIT.md).

## P0 — prevent destructive A1 store rebuilds (2026-10-19)

- **Implemented:** new A1 builds publish only after post-write structural verification and refuse to overwrite any existing store directory.
- **Implemented:** isolated `--output-dir` and `--train` command options for versioned research datasets. Temporary staging directories are removed on exceptions.
- **Verification:** standard-library synthetic invariants and four isolated local checks; no CI or data-dependent full-source build was run.
- **Open evidence:** [store immutability protocol](A1_STORE_BUILD_IMMUTABILITY.md). A separate full-source build must retain its own source manifest and evaluation signature.

## P0 — COCONUT neutral-parent structure identity (2026-10-19)

- **Implemented:** de-salted COCONUT candidates now carry a consistent parent SMILES, standard InChIKey14, exact mass and molecular formula. Malformed or charged parents are excluded instead of replaced by placeholder molecules.
- **Implemented:** directly supplied neutral parent records take precedence over salt-derived duplicates; synthetic regression tests cover structure identity and source annotation selection.
- **Testing:** focused local RDKit/pandas checks passed; **CI was not invoked**, in accordance with the user's instruction.
- **Pending:** rebuild an isolated authorized COCONUT candidate-pool variant, regenerate aligned fingerprints and rank features, record a source manifest, and test on molecule-disjoint C1/C2/C3 populations. [Scientific identity audit](COCONUT_PARENT_IDENTITY_AUDIT.md).

## P0 — Kaggle-source access and D2 fingerprint identity (2026-10-19)

- **Implemented:** local `KAGGLE_API_TOKEN` availability and optional authenticated file-listing preflight; API secrets never enter arguments, logs, reports or the public repository.
- **Implemented:** versioned fingerprint provenance checks tie the exact ordered pool key, SMILES and mass to packed Morgan bytes, eliminating silent reuse of stale fingerprints after COCONUT parent-identity corrections.
- **Implemented:** source pool and fingerprint builders refuse to overwrite frozen variants; `run_validation.py --require-fingerprint-provenance` rejects stale or missing D2 artifacts when explicitly required.
- **Verification status:** synthetic offline tests were added; Kaggle account authorization, rules acceptance, actual dataset download and full model evaluation remain **unverified**. No CI was requested.
- See [Kaggle connection and research snapshot protocol](KAGGLE_ACCESS_AND_FP_PROVENANCE.md).

## P0 — D2 fingerprint decision eligibility and candidate mass domain (2026-10-19)

- **Implemented:** `align_runs` requires verified, source-bound Morgan fingerprint provenance for every D2 report before an `ACCEPT`/`REJECT` verdict; historical reports without a sidecar remain exploratory and contradictory fingerprint evidence raises an error.
- **Implemented:** matched D1-only comparisons remain eligible without unused fingerprints, and D1–D2 ablations can share a complete superset source snapshot when the D2 sidecar is verified.
- **Implemented:** `CandidatePool` rejects present nonfinite/nonpositive masses and validates positive finite query mass plus nonnegative finite ppm tolerance. Intentional missing-mass filtering is preserved.
- **Local checks:** isolated synthetic NumPy/Pandas and D2 evidence smoke cases passed. Per user instruction, no CI or full scientific suite was triggered.
- **Unverified:** Kaggle dataset access, C1/C2/C3 rerun and metric improvement. See [scientific data-quality audit](D2_FINGERPRINT_EVIDENCE_GATE.md).

## P0 — reject malformed Morgan inputs and audit frozen vectors (2026-10-19)

- **Implemented:** fingerprint construction now uses the same `CandidatePool` invariants as inference and fails closed for invalid RDKit SMILES, zero-heavy-atom molecules and all-zero Morgan fingerprints.
- **Implemented:** a read-only, batched `audit_pool_fingerprints.py` checks candidate-key order, packed-vector shape, zero-vector prevalence and optional provenance manifest on immutable existing research artifacts.
- **Tests:** synthetic generator failure/no-output guards and aggregate-quality regressions added. Local isolated RDKit and NumPy smoke checks passed; full scientific suite and GitHub CI were **not run** per instruction.
- **Pending:** run the audit against authorized full candidate pools, rebuild only separately tagged variants, and independently evaluate C1/C2/C3. See [fingerprint quality protocol](FINGERPRINT_ZERO_VECTOR_AUDIT.md).

## P0 — multi-provider source identity and immutable pool creation (2026-10-19)

- **Implemented:** before grouping by standard InChIKey14, drop incomplete *whole source identity rows* and reject malformed keys, invalid present masses and serious within-key mass conflicts (threshold max[0.01 Da, 20 ppm]).
- **Implemented:** publish the new pool Parquet only after a staged write and row-count footer check, with exclusive destination creation; old candidate pools cannot be overwritten by a racing writer.
- **Local validation:** five focused synthetic NumPy/pandas/stdlib tests passed; a separate real Parquet readback test was blocked by missing `pyarrow` in the isolated local environment. The full suite and GitHub CI have not been run.
- **Remaining source gate:** authorized train/COCONUT/optional-provider data must be reconciled, conflicts reviewed, and separately tagged pool/fingerprints validated before any scientific score claim. See [source audit protocol](CANDIDATE_SOURCE_RECONCILIATION.md).

## P0 — generated-candidate fingerprint and identity safety (2026-10-19)

- **Implemented:** D3 overlay registration requires balanced key/SMILES inputs, RDKit-valid molecular structures, the correct generated 'G'+InChIKey14 identity and nonzero packed Morgan fingerprints.
- **Implemented:** new generated rows are staged before modifying the pool overlay; malformed later rows leave existing state unchanged. `packed_fps` validates entire source fingerprint row count and 2048-bit uint8 shape.
- **Verification:** 5 isolated local RDKit/NumPy cases passed; additional synthetic tests committed. No GitHub CI or full licensed data rerun.
- **Pending:** separately signed D3-enabled novelty evaluation and negative controls on authorized data. See [generated candidate integrity protocol](GENERATED_CANDIDATE_INTEGRITY.md).

## P0 — D2 analog-reference cache integrity (2026-10-19)

- **Implemented:** opt-in `run_validation.py --fresh-analog-index` bypasses all D2 sparse-matrix cache reads/writes and computes an attested digest over selected reference IDs, CSR matrix and model configuration.
- **Implemented:** D2 statistical model comparisons now require independent cache-free index evidence, matching source input inventory and fingerprint sidecar; historical cached D2 reports remain exploratory.
- **Coverage:** synthetic cache-free/cache-loaded and contradictory provenance tests; existing D1 and D1-vs-D2 comparisons remain usable when the D2 run is verified.
- **Remaining:** authorized full dataset index rebuild, file manifest comparison, model-independent replication and full scientific regression suite; **CI deliberately not invoked**. See [D2 index evidence protocol](D2_FRESH_INDEX_EVIDENCE.md).

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
