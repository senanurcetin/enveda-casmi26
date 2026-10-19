# Peak content quality audit — Enveda CASMI26

**Research-only update — 2026-10-19**

This non-commercial scientific change addresses a gap between the structural
integrity of the A1 spectrum store and the numeric correctness of its contents.

## Finding

The existing `inspect_store_layout` validates offsets, file lengths and
metadata peak counts, without examining peak values. A store can pass these
checks while containing physically impossible fragment m/z values
(`m/z <= 0`), non-finite numbers or nonpositive intensities.

The peak-cleaning implementation previously removed non-finite values and
nonpositive intensities, but did not explicitly exclude nonpositive m/z.
Negative m/z values are not physical MS/MS fragments and can also create
invalid sparse feature-vector bin indices in D2 analog retrieval.

## Implementation

- Both `clean_peaks` (NumPy reference) and `clean_peaks_fast` (Numba) now
  remove fragments where `m/z <= 0`, in addition to their existing filters.
- `scan_store_peak_content` uses **bounded-memory streaming** of the two
  float32 peak arrays, after checking file-length/offset consistency. It
  reports only counts for non-finite m/z, nonpositive m/z, non-finite
  intensities, nonpositive intensities and unique invalid peak pairs.
- Existing `audit_store_integrity.py` supports optional
  `--scan-peaks` and `--chunk-peaks`. Its **default remains the fast
  structural audit**. A requested full scan uses more disk I/O but no
  full-dataset in-memory peak array, and exits nonzero if invalid numeric
  peak pairs are detected.
- Synthetic tests cover both cleaner implementations, chunk boundaries,
  empty files, NaN/infinity, nonpositive values and invalid chunk settings.

## Read-only procedure

```sh
python scripts/audit_store_integrity.py \
  --store /path/to/artifacts/store \
  --scan-peaks --chunk-peaks 262144 \
  --output runs/peak_content_quality.json
```

Use only a locally authorized A1 store. The report is aggregate and contains
no raw spectrum or held-out truth label. The output uses exclusive file
creation, so a previous research record is not silently overwritten. The scan
caps each batch at **1,048,576 peak pairs** to bound memory even when a user
passes an unusually large `--chunk-peaks` value.

**Interpretation:** `VALID_CLEANED_PEAKS` means the inspected stored peak
values were finite and strictly positive for both m/z and intensity. It
does not establish instrument calibration, chemical accuracy, absence of
duplicate spectra, cross-library independence, or Kaggle/BigQuery provenance.

## Reproducibility constraints

Changing the cleaner's admissible input domain affects only spectra
containing invalid nonpositive fragment masses; it does not change normal
positive-m/z fragment spectra by design. Nevertheless, historical A1 stores
are immutable experiment artifacts: **do not overwrite** the original store
or reinterpret earlier leaderboard / validation MRR scores as corrected.
Rebuild only in a separately versioned research run, verify an input
SHA-256 manifest and evaluate C1/C2/C3 on identical held-out molecules.

The source dataset was not made available for a real data scan in this
development turn. This is code and synthetic-test coverage, **not** a
measured data quality result or new Kaggle score. Public MRR@25 0.422
(2026-10-08) is historical and unchanged.

## CI usage during this sprint

At the user's direction, **automatic CI was not invoked for this work**.
Commit messages use the GitHub Actions `[skip ci]` directive. The existing
research workflow is retained for explicitly authorized future checks.
