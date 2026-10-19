# Generated candidate overlay integrity (CASMI26)

**2026-10-19 | Non-commercial scientific research**

## Finding

`CandidatePool.register_generated`, used by the optional D3
biotransformation branch of `BaselineEngine`, previously zipped
generated keys and SMILES without checking equal lengths. An invalid
SMILES could be saved as an apparently valid generated candidate
with an **all-zero Morgan fingerprint**.

Moreover, generated identity and fingerprints were added row by row.
If a later row in a batch was invalid, the overlay could be left in a
partially updated state. Although the regular source-pool fingerprint
builder now rejects invalid SMILES, that safeguard did **not** previously
cover the generated-candidate overlay.

## Fix

`register_generated` now:
- Checks key and SMILES list lengths before processing.
- Requires a parseable RDKit molecule with at least one heavy atom.
- Requires the generated key to equal **'G' + standard InChIKey14**
  derived from that specific molecular structure, as used in
  `BaselineEngine._generated`.
- Requires nonempty Morgan radius-2, 2048-bit fingerprints.
- Stages every newly generated key/SMILES/packed-vector tuple in memory,
  then publishes the entire batch only if all records validate.
- Reports only the invalid **row index** on failure. Molecular strings,
  formulas and held-out labels are not emitted.

`packed_fps` now rejects mismatched row counts, non-uint8 vectors
and shapes other than `[n_pool, 256]` before mixing frozen candidate
vectors with generated overlays. Requests for unknown keys raise an
explicit error rather than returning an unrelated or zero vector.

## Scientific boundaries

The baseline uses frozen source-pool records for C1/C2/C3 and an
optional generated overlay for class-3 hypothesis testing. This change
does not regenerate source fingerprints or alter any historical
candidate artifacts. Invalid generated candidates now fail closed
instead of contributing misleading similarity features.

Because this is a model-behavior change on formerly invalid input,
new D3-enabled experiments must be compared on identical independently
signed queries, verified source inventory, and controlled model settings
before claiming an improvement.

Local isolated RDKit/NumPy smoke tests: **5 passed** (valid
candidates, invalid late row without partial update, identity mismatch,
unbalanced input lengths and malformed fingerprint array shapes).
Additional synthetic repository tests were committed. The full scientific
test suite and GitHub CI were **not run**, respecting the user's
instruction to avoid CI.

Actual Kaggle/BigQuery source data and leaderboard evaluation were
not accessed. Historical public MRR@25 **0.422 (2026-10-08)** remains
unchanged. No productization, API, service or commercialization.
