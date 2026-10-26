# Executed public-example overlap audit — 24 October 2026

**Finding:** every spectrum in the mounted public example `test.parquet` has an exact stored numerical-content match in the competition training input. This rules out treating those examples as an independent unseen-spectrum validation set. It does not establish contamination of the hidden evaluation set.

## Executed evidence

| Measurement | Result |
|---|---:|
| Training rows scanned | 2,539,608 |
| Public example spectra | 1,213 |
| Exact precursor + complete ordered m/z/intensity matches | 1,213 / 1,213 |
| Also matching adduct and ionization mode | 1,213 / 1,213 |
| Molecules with all public spectra matched | 400 / 400 |
| Exact-precursor candidate training rows | 14,620 |
| Coarse-bucket candidate pairs after precursor filtering | 1,215 |
| Rejected coarse-bucket pairs | 2 |
| Source library of the exact matches | enveda-180 |

The [Kaggle notebook](https://www.kaggle.com/code/senanuretin/casmi26-exact-public-example-overlap-audit) completed with CPU execution and internet disabled. Its downloaded source matches the tracked bundle. [Aggregate result and complete input/source hashes](../../reports/public_overlap_audit.json) · [Authenticated publication readback](../../reports/public_overlap_publication.json).

Run timestamp: 23 October 2026, 22:03 UTC / **24 October 2026, 01:03 Europe/Istanbul**. The run used Python 3.13.15 and PyArrow 23.0.1 on Kaggle. Local Parquet regression tests use the installed PyArrow 24.0.0; the executed environment is recorded separately, not inferred from the local environment.

## Protocol

The public examples are indexed by rounded precursor, peak count and first fragment only to shortlist comparisons. Arrow scans the complete training input in batches of 8,192 rows. Exact unrounded precursor equality is a necessary prefilter; all selected candidate pairs then compare every ordered fragment m/z and normalized intensity without tolerance or rounding. Missing, nonfinite, empty or misaligned public arrays fail the audit. Source-library attribution is required for any exact match.

Counts distinguish target spectra from matching reference pairs. Multiple copies of one spectrum cannot inflate the number of matched public targets. The molecule denominator includes all targets, and a molecule is fully matched only when all its public spectra match. Adduct and ionization mode are checked separately from numeric content.

Only aggregates, software versions and SHA-256 hashes are published. Raw identifiers, structures and peak arrays are absent from the report. The 3,033,286,496-byte training file and 4,848,729-byte public example file remain on Kaggle; file hashes identify this specific mounted snapshot.

## Interpretation limits

- Exact stored numerical content is stronger evidence than the historical coarse keys. The two rejected coarse candidates demonstrate why that distinction matters.
- Numeric identity and metadata agreement do not independently authenticate common acquisition provenance.
- This does not measure generalization of any model and does not validate the earlier public leaderboard scores.
- Hidden evaluation inputs were not read. No model was trained, rescored or submitted.
- Changed precision, peak order or preprocessing can produce non-matches even for related acquisitions; this is an exact-content audit, not approximate spectral search.

## Reproduce

```bash
python scripts/build_overlap_notebook.py
# Run the prepared bundle on Kaggle with the existing competition attachment.
python scripts/check_overlap_publication.py
```

The readback command refuses incomplete, private, source-mismatched or differently attached notebooks, downloads only the aggregate JSON, validates count invariants and code hashes, and records the artifact hash. It does not publish, execute a new run or submit predictions.

For already-authorized local files:

```bash
python scripts/audit_public_overlap.py --train data/train.parquet --public data/test.parquet --output runs/public-overlap-new.json
```

An existing output is never overwritten. The default workflow uses no paid API, cloud warehouse or GPU.
