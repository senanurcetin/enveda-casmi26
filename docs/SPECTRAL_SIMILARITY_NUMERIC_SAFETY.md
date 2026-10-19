# Spectral-similarity numerical reliability audit

**2026-10-19 | Enveda CASMI26 public scientific research**

## Finding

Entropy and cosine similarity assume correctly shaped, positive-intensity
MS/MS spectra. The historical implementation normalized every nonempty
intensity vector by its sum, then performed logarithms or square roots
without guarding against a zero total, non-finite values, negative
intensities or unequal m/z/intensity lengths.

While the standard cleaning path normally filters invalid peaks, the
similarity functions themselves are used for validation and library/analog
comparisons. If malformed arrays get through, a division by zero or
`0 * log(0)` could create `NaN` spectral similarities and contaminate
candidate ranking or experiment summary statistics.

## Implemented safeguards

- Shared Numba-compiled numerical input validation checks equal-sized 1D
  arrays, finite nonnegative peak m/z, and finite nonnegative intensities.
  Malformed array lengths are rejected rather than silently truncated.
- A spectrum with no fragment peaks, no positive intensity, non-finite
  values or negative values yields similarity **0**, not `NaN`.
- The entropy-weighting function guards against zero or non-finite
  normalization. Matched zero-intensity peaks contribute zero without
  evaluating invalid `log(0)` products.
- Cosine similarity returns zero when either spectrum has a zero
  normalization denominator.
- Existing valid positive-intensity spectra continue through the same
  entropy-weighting and peak-matching algorithms; no ranking/weighting
  strategy is intentionally changed.

## Verification

Synthetic numerical regressions cover both scoring functions for empty
arrays, zero-intensity spectra, NaN, infinity, negative values, malformed
lengths, shared zero-intensity peaks and ordinary valid-spectrum invariants.
The pinned scientific CI job now runs the **existing** peak-cleaning and
entropy/cosine tests in addition to the new cases.

**Scientific limits:** Returning zero for an invalid spectrum is a
conservative fail-safe, not a substitute for counting data-quality issues
at ingestion. This does **not** prove the absence of corrupted spectra in
the actual original CASMI dataset; the dataset itself was not rerun.
Validation on the full store and source/instrument-specific quality reports
remain required.

No new Kaggle submission, BigQuery query, model training, leaderboard score,
productization, SaaS, service or commercial feature was introduced. The
historical public MRR@25 **0.422 from 2026-10-08** is unchanged.
