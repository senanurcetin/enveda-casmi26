# MS/MS peak-store integrity — fail-closed structural validation

**Research audit, 2026-10-19.** This work belongs to the public, non-commercial
Enveda CASMI26 scientific codebase. It does not modify a trained model, deploy
an inference endpoint, use private Kaggle answers or claim a new leaderboard
result.

## Finding: incomplete validation of the A1 store

Before this change, `SpectrumStore` loaded the metadata table, flat peak
files and offset vector, then asserted only that the number of offsets equaled
the metadata row count plus one. This did not detect:

- A nonzero initial offset, negative or non-monotonic offsets
- A terminal offset beyond the end of the peak arrays, or a truncated array
- Inconsistent m/z and intensity array lengths
- Binary files truncated in the middle of an int64 offset or float32 value
- Per-spectrum metadata `n_peaks` counts disagreeing with actual peak slices
- Negative/out-of-range spectrum row IDs accidentally resolving to the last
  spectrum through Python negative indexing

A structural fault can silently misassociate spectrum peaks and molecules,
invalidate C1/C2/C3 evidence and taint model comparisons even when individual
query identities appear correctly signed.

## Controls implemented

`src/casmi/spectra/store_integrity.py` validates the file structure
*before* memory mapping peak data. It verifies offset/peak byte alignment,
consistent lengths, a starting offset of zero, nonnegative and monotone
offsets, a final offset equal to flat-peak count, metadata row count and
the `n_peaks` of every spectrum.

The validator reads the int64 offset vector, but **does not scan the
full peak arrays**. This keeps the structural check proportionate to
the number of spectra rather than the total number of fragmentation
peaks. The store constructor now fails closed with `ValueError` on
broken layouts instead of a single optional Python `assert`.
Spectrum access also rejects negative and out-of-range indices, and
a valid zero-length peak array no longer causes a memmap exception.

The builder `scripts/build_store.py` rechecks newly persisted files
after writing metadata, helping catch incomplete A1 production earlier.

## Reproducible read-only audit

With the already-approved local artifact tree available and the Parquet
runtime installed:

```sh
python scripts/audit_store_integrity.py \
  --store /path/to/artifacts/store \
  --output runs/store_layout_audit.json
```

The output contains the number of spectra, total peaks, empty spectra,
maximum peaks per spectrum, and a structural verification status.
No raw peak data, molecular labels, private answers, absolute source paths
or credentials are included. The output file is created exclusively and
will not overwrite a prior experiment record.

## Scientific scope and limits

- **Structural integrity is not scientific accuracy.** This does not
  verify peak m/z ranges, intensity values, chemical formula correctness,
  instrument calibrations, overlap/leakage or source metadata truth.
- It also does **not** replace the full streaming SHA-256 input manifest:
  equal lengths and offsets do not prove byte-for-byte data identity.
  Use the [input snapshot protocol](INPUT_SNAPSHOT_MANIFEST.md) for that.
- Local successful CI tests use small **synthetic** corrupted files and
  mocked metadata; the original multi-million-spectrum CASMI store was not
  opened or revalidated during this change.
- Historical public MRR@25 **0.422 from 2026-10-08** remains unchanged.
  No live Kaggle or BigQuery remeasurement was performed.

## Acceptance tests

The research-assets CI job runs NumPy-only tests for offset corruption,
empty files, unaligned binary lengths, contradictory metadata counts and
truncated arrays. The scientific stack CI job exercises the actual
`SpectrumStore` constructor and peak access on synthetic raw files.
