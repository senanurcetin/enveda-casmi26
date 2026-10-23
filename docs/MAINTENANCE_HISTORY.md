# Archived research maintenance notes

Preserved from the pre-23 October 2026 README. These are dated implementation notes, not a new evaluation record.

### Research reliability update (18 October 2026)

- The scientific run comparator now reports **approximately simultaneous paired-bootstrap confidence intervals** over the fixed C1/C2/C3 novelty-mixture grid, rather than interpreting several pointwise 95% confidence intervals as one joint 95% statement. The correction does **not** address repeated experiment selection or distribution shift; see [the original experiment log and caveat](EXPERIMENTS.md).
- Invalid or blank SMILES rank slots are preserved in local MRR scoring, with independent regression tests; no malformed candidate can be silently removed to improve its rank.
- Validation comparisons reject mismatched spectrum-signature algorithm versions, and legacy unverified runs cannot establish positive/negative scientific conclusions.
- These are **non-commercial, research-only** improvements. No new training, Kaggle submission, competition score or live warehouse data is claimed.

## Verified fresh D2 index comparisons (2026-10-19)

A [cache-free D2 analog-reference index protocol](docs/D2_FRESH_INDEX_EVIDENCE.md) prevents models using unverified historical sparse-index caches from receiving conclusive ACCEPT/REJECT judgments. New research runs can opt into `--fresh-analog-index` with `--analog`, verified fingerprints and a complete source manifest. This rebuild is resource-intensive and does not alter frozen leaderboard results.

## Generated molecular candidate safeguards (2026-10-19)

The [D3 generated-candidate integrity audit](docs/GENERATED_CANDIDATE_INTEGRITY.md) prevents invalid SMILES or mismatched generated InChIKey14 values from being stored as empty Morgan fingerprints. Generated overlay updates are now batch-atomic, and candidate fingerprint matrices must be correctly row-aligned uint8 arrays. Historical model scores and stored source candidates are untouched.

## Cross-source candidate identity consistency (2026-10-19)

The [candidate-source reconciliation audit](docs/CANDIDATE_SOURCE_RECONCILIATION.md) guards against combining SMILES and neutral masses from different provider records and flags large same-InChIKey14 mass discrepancies. The A3 candidate pool is now staged and structurally checked before exclusive publication, preserving frozen research artifacts. No independent leaderboard improvement or productionization is claimed.

## Fail-closed Morgan fingerprint generation and audit (2026-10-19)

[Read-only fingerprint quality auditing](docs/FINGERPRINT_ZERO_VECTOR_AUDIT.md) is now available to detect previously signed but chemically invalid all-zero candidate vectors. New fingerprint builds reject malformed SMILES, zero-heavy-atom structures and empty Morgan vectors before output files are written, and use the same candidate mass ordering as inference. Historical artifacts are untouched; no Kaggle performance improvement is claimed.

## D2 scientific decision gating and candidate mass checks (2026-10-19)

[Verified D2 fingerprint evidence and candidate mass validation](docs/D2_FINGERPRINT_EVIDENCE_GATE.md) now prevents ACCEPT/REJECT judgments from D2 runs lacking source-bound fingerprint records, even when query signatures and a general source inventory match. D1-only and verified D1–D2 ablations remain supported. Candidate mass windows reject impossible values instead of propagating invalid ppm calculations. No CI or new competition-score claim.

## Kaggle access and candidate-fingerprint provenance (2026-10-19)

The [secure Kaggle access and fingerprint provenance protocol](docs/KAGGLE_ACCESS_AND_FP_PROVENANCE.md) adds a local, read-only Kaggle CLI preflight without embedding credentials. Rebuilt COCONUT candidate pools and Morgan arrays now have a content-bound SHA-256 fingerprint manifest; historical pool/fingerprint files cannot be silently overwritten. D2 validation can require this provenance for new research runs. No Kaggle download, submission or leaderboard experiment was performed here.

## COCONUT source identity correction (2026-10-19)

The [COCONUT neutral-parent identity audit](docs/COCONUT_PARENT_IDENTITY_AUDIT.md) corrects a research data-quality issue: removing salt fragments previously changed candidate SMILES and mass without updating the original compound InChIKey. Parent SMILES, standard InChIKey14, molecular formula and exact mass are now derived from one neutral molecular graph. Frozen candidate pools and leaderboard results are not modified.

## Immutable A1 spectrum-store builds (2026-10-19)

A [staged research-store build safeguard](docs/A1_STORE_BUILD_IMMUTABILITY.md) now refuses to overwrite existing MS/MS peak stores. New builds are written into isolated temporary directories, structurally verified and only then published, preserving historical input artifacts and experiment reproducibility. No CI, productization or new leaderboard data is implied.

## D1 ionization-mode separation (2026-10-19)

A [D1 scientific audit](docs/D1_IONIZATION_MODE_AUDIT.md) corrects reference matching for mixed positive/negative ionization acquisitions. Each spectrum now uses its own ionization-mode reference group; missing-mode references contribute zero rather than creating cross-mode matches. The research-only change has synthetic tests but no new CASMI leaderboard result.

## Numeric peak-content auditing (2026-10-19)

A [read-only chunked peak-content audit](docs/PEAK_CONTENT_QUALITY_AUDIT.md) now detects invalid m/z and intensity values in the stored MS/MS arrays. The standard store-layout check remains fast; full peak-content scanning is optional. Both peak cleaners also reject nonpositive fragment m/z. Historical scores are unchanged and no data-dependent accuracy improvement has been claimed.

## Spectral-similarity numerical safeguards (2026-10-19)

The [MS/MS entropy/cosine numerical audit](docs/SPECTRAL_SIMILARITY_NUMERIC_SAFETY.md) ensures invalid, zero-intensity and non-finite spectra cannot silently yield NaN channel similarities. Previously untested numerical boundary cases are now covered in the pinned scientific CI job. This is data-quality protection, not a proven model-score improvement.

## Peak-store structural integrity (2026-10-19)

The [A1 spectrum-store structural audit](docs/SPECTRUM_STORE_INTEGRITY.md) now checks binary offset ordering, m/z/intensity file alignment, terminal peak counts and per-spectrum metadata before reference arrays are used. Invalid row indexing fails closed. This is a read-only data-integrity safeguard, not a new prediction result or complete scan of peak contents.

## Scientific subgroup diagnostics (2026-10-19)

A [molecule-level slice analyzer](docs/ERROR_SLICE_DIAGNOSTICS.md) now supports label-free adduct, ionization-mode, precursor-mass and spectra-count summaries, preserving query pairing and suppressing sparse subgroups. All effects are descriptive only; no new validated model scores or production features.

## Input-complete model comparison safeguard (2026-10-19)

Scientific paired model comparisons now require a verified *complete minimum* local input inventory as well as matching molecule/spectral signatures before returning an `ACCEPT`/`REJECT` decision. Earlier runs without these manifests remain exploratory. See [source-inventory evidence rules](docs/INPUT_SNAPSHOT_MANIFEST.md). No Kaggle/BigQuery model score has been recomputed.

## Research input snapshot integrity (2026-10-19)

The [local SHA-256 snapshot protocol](docs/INPUT_SNAPSHOT_MANIFEST.md) allows researchers to record and verify the approved spectrum store, candidate pools and fingerprints without committing dataset bytes or absolute paths. Validation can optionally fail closed on changed source artifacts. This is reproducibility infrastructure, not proof of Kaggle/BigQuery authenticity or improved prediction scores.

## Source-verified negative-control comparisons (2026-10-19)

A [diagnostic-only control comparator](docs/SPECTRAL_NEGATIVE_CONTROLS.md) now checks original molecule/spectrum identities before calculating paired bootstrap differences between unperturbed and blank/shuffled spectra. Results cannot produce model-improvement ACCEPT/REJECT decisions or new Kaggle claims. This remains non-commercial research work.

## Spectral evidence diagnostics (2026-10-18)

A [non-commercial negative-control protocol](docs/SPECTRAL_NEGATIVE_CONTROLS.md) adds blank-fragment and deterministic intensity-shuffle modes to the knockout validation runner. Control reports preserve original source fingerprints and may **never** be promoted to ordinary model-improvement ACCEPT/REJECT verdicts. These are scientific diagnostic tools only; no new Kaggle results are claimed.

## Controlled D2 cross-library index experiment (2026-10-18)

The [source-balanced D2 experiment protocol](docs/D2_SOURCE_BALANCED_EXPERIMENT.md) introduces an **opt-in** reference sampler that preserves the highest-quality spectrum while reserving independent source-library coverage when the index cap permits. The baseline policy is unchanged. A distinct, provenance-checked cache keeps the new experiment separate from historical inference artifacts; the metadata audit can compare both index-selection policies. This is a validation and reproducibility improvement, **not** a new MRR result. Productization remains excluded.

## Scientific integrity update — D2 analog shortlisting (2026-10-18)

The [D2 knockout and prefilter audit](docs/D2_PREFILTER_KNOCKOUT_AUDIT.md) documents a second validation flaw: excluded self/reference spectra were removed **after** approximate top-k selection, so invalid hits could crowd out eligible independent analogs. The corrected algorithm applies knockout, adduct and mass eligibility **before** top-k and is covered by numerical plus direct `AnalogIndex.search` regression tests. Historical competition metrics remain unchanged. No new Kaggle run, submission or product work is implied.

## Scientific evaluation update — C1 reference leakage (2026-10-18)

The [non-commercial research assessment](docs/CASMI26_RESEARCH_ASSESSMENT.md) identifies and fixes a C1 holdout flaw: when only a subset of a molecule's source spectra was used as a simulated query, its *unselected* same-structure reference spectra could remain in the search index. New exclusion bookkeeping now hides every same-library reference, including tautomer-equivalent standard keys, while retaining legitimately independent cross-library references. Query signatures have been versioned to `casmi-query-spectra-v2`, and regression tests cover sampled queries.

**Existing public and historical validation scores have not been recomputed.** This is a methodology fix, not a newly demonstrated model-score increase. Project scope remains scientific research only; no productization.

## NEXUS Research Intelligence (2026-10-18 extension)

A [NEXUS 2.0 research dashboard](site/index.html) now presents the **historical 8 October** model results, experimental provenance, and the source-library profile in `reports/data_profile.json`. It uses the approved NEXUS SVG and works locally with `python -m http.server 8000` (open `/site/`). **It does not represent live Kaggle or BigQuery results and does not change the frozen model or leaderboard claims.**

**Research-only scope:** No productization, service deployment, monetization, or live dashboard pipeline is planned. The staged implementation and scientific acceptance gates are in [RESEARCH_WORK_PROGRAM.md](docs/RESEARCH_WORK_PROGRAM.md). The current sprint adds query-identity and source-spectrum fingerprints to validation reports, disallows unsupported paired bootstrap claims for legacy runs, and provides offline BigQuery profile reconciliation. These updates do not change any model score or trigger a Kaggle submission.

New quality utilities include [`scripts/audit_submission.py`](scripts/audit_submission.py) for submission shape, exact test molecule coverage, candidate budget and optional canonical structure checks; a read-only [BigQuery data quality SQL template](sql/bigquery/library_quality.sql) with a safe renderer; and dependency-free regression tests in GitHub Actions. See [implementation details and limitations](docs/NEXUS-CASMI26.md).
