# Enveda CASMI26 — Scientific codebase assessment

Assessment date: **2026-10-18**. This review concerns the existing **non-commercial scientific research codebase**. It is not a market, product or deployment review.

**Current evidence update — 25 October 2026:** the [full-source frozen-model/control study](../reports/revalidation_review.md) has completed, with authenticated publication readback, all 250 V-in molecules retained, verified fingerprints and a freshly rebuilt D2 index. It closes the historical review's retrospective rerun/negative-control gap for this population. v16 C1/C2/C3 MRR is 0.8324/0.6241/0.2374; blank fragments improve simple-model C2, an adverse diagnostic reported explicitly. Historical leaderboard scores remain unchanged. Full weight-training lineage, untouched external validation and live BigQuery reconciliation remain distinct limitations. The dated assessment below records the original review state.

## Scope and evidence

Inspected: `README.md`, `EXPERIMENTS.md`, the then-current historical handoff (retained in Git history), C1/C2/C3 knockout logic, query sampling, rank-data construction, analog/library indices, the paired-bootstrap comparator and repository CI. The latest archived leaderboard figures are from **2026-10-08**, not new independent model evaluations. No new Kaggle training/leaderboard result or live BigQuery verification was performed in this review.

| Research criterion | Assessment | Evidence |
| --- | --- | --- |
| Scientific architecture | Mature experimental pipeline | D1 library, D2 analog, D3 fingerprint and D4 fragment channels; ranker and candidate generation |
| Experimental traceability | Strong relative to a notebook-only entry | Detailed positive and negative experiments in `EXPERIMENTS.md`; fixed reported baseline snapshots |
| Validation transfer | Significant unresolved limitation | V-in / V-real weight selection previously disagreed with Kaggle public results; novelty classes and instrument shift remain |
| Pairing and uncertainty | Improved, not fully solved | Verified query IDs, per-query spectra fingerprints, simultaneous fixed-grid paired confidence bands |
| Data quality | Partially implemented | Source-library quality SQL and offline reconciliation; actual BigQuery dataset and input snapshot not independently confirmed |
| Leakage resistance | **Critical sampled-query C1 problem corrected in code** | See source-library exclusion regression below; full dataset evaluation still required |
| Reproducibility | Improved CI, incomplete end-to-end | Dependency-free checks, synthetic statistics, pinned official RDKit metric tests; complete GPU/Kaggle pipeline not rerun |

**Research readiness:** The foundation supports further controlled experiments, but evidence does not yet support claims that newly added code improves prediction quality.

## Critical finding — C1 leakage under spectral subsampling

Before this sprint, `build_queries(..., max_spectra=N)` kept only sampled spectrum rows in `Query.spectrum_rows`. The C1 knockout used that same set as the reference exclusion list. For a structure with more acquired spectra than `N`, its unselected, still-present spectra remained eligible as library/analog references. This can artificially increase direct/self-spectrum matching on the validation set and make comparisons unreliable. A second risk is leaving distinct standard-key representations with the same tautomer-canonical metric key in the query's own library.

### Remediation implemented

- `Query.spectrum_rows` keeps **actual selected input spectra**, while the new `Query.source_library_rows` records **all source-library acquisitions for the same standard structure key**, including spectra not eligible for query sampling.
- The `knockout_for` C1 exclusion combines the sampled rows with all source-library reference rows. C2/C3 also combine them with the indexed cross-library matches.
- `include_equivalent_source_rows` adds tautomer-equivalent standard-key variants to C1 exclusions in `scripts/run_validation.py`, without deleting their query input spectra or unintentionally dropping valid references from other libraries.
- `scripts/build_rank_data.py --set vin` also hides same-library tautomer-equivalent D1/analog references in C1 ranker-validation tables; existing explicit cross-library training logic remains.
- The query provenance recipe is advanced to `casmi-query-spectra-v2` and signs **both** selected query spectra and complete source-library exclusion rows. Comparisons with v1 query signatures are rejected by the existing version check.

### Regression evidence

Synthetic data include two selected-source-compatible spectra, an additional same-key spectrum excluded from query sampling, a tautomer-equivalent source-library reference, and an allowed other-library reference. With `max_spectra=1`, tests require C1 to hide all four own-library references, leave the other-library reference available, and require C2/C3 to hide the complete set. Additional tests verify malformed row references, altered signatures and deterministic, unchanged query input spectra.

**Caveat:** The regression suite demonstrates correctness for its cases. It does **not** quantify the impact on historical MRR scores; doing that requires rerunning the competition-provided training inputs, candidate pools and channel models under the corrected protocol. Historical metrics must remain unchanged and explicitly dated until those reruns occur.

## Prioritized next research experiments

1. **P0 — validation rerun and contamination audit:** Rebuild C1/C2 reference features with the corrected knockout, preserve the same query population and log the new v2 manifest. Measure how many queries previously retained own-library references and the per-query change in MRR@25.
2. **P0 — data provenance:** Confirm actual BigQuery table schema and Kaggle input versions, snapshot hashes, instrument labels and candidate pool lineage. Compare historical profile counts with a read-only BigQuery aggregate when accessible.
3. **P1 — negative controls:** Run blank-spectrum, shuffled-peak and membership-only baselines for all novelty classes. Examine direct-match versus analog propagation under intentionally withheld structures.
4. **P1 — class/source diagnostics:** Report C1/C2/C3, adduct, source library, collision energy and timsTOF versus other instruments separately. Use structure-disjoint holdouts; select on independent validation, not the public leaderboard.
5. **P2 — single-hypothesis model work:** Only after evidence gates pass, test one change at a time against the frozen baseline with matched query identities, simultaneous bootstrap bands and a separate holdout.

## Non-commercial boundary

There is no product road map, application commercialization, billing, customer workflow, deployment target or paid service in this project. NEXUS branding remains a static research identity for provenance and communicating historical experiments.

## Local verification

```sh
python -m pytest tests/test_knockout_leakage.py -q
python -m pytest tests/test_foundations.py -q
python -m unittest discover -s tests -p '*_stdlib.py' -v
```

The first two tests need the pinned scientific environment. GitHub Actions executes them in its official-metric regression job, plus existing offline smoke tests. Data-dependent accuracy improvement is explicitly **not** claimed.
