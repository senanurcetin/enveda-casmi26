# Candidate-source identity reconciliation and immutable pool publication

**Enveda CASMI26 — research-only, 2026-10-19.** This is a scientific
integrity safeguard, not a new model, production feature or Kaggle result.

## Finding: cross-provider records could be spliced or misidentified

The A3 pool aggregates training, COCONUT, and optional PubChemLite,
ChEBI, LIPID MAPS and PubChem-derived records by standard
`InChIKey14` connectivity block. It previously called
`groupby.first()` separately on `smiles` and `mass`.

Pandas `first` skips nulls on each column separately: when one source row
has SMILES but no mass, a later provider's mass could be attached to the
first provider's SMILES. More importantly, materially different masses
reported under the same connectivity block were silently collapsed
under source priority. This can corrupt neutral mass windows, molecular
fingerprint interpretation and C1/C2/C3 candidate selection.

## New source-identity gate

Before aggregation, `reconcile_source_rows`:

- Requires complete *same-row* `pool_key`, nonblank `smiles` and `mass`
  identity fields. Incomplete source identity records are counted and
  excluded. The original source priority among remaining rows is retained.
- Rejects malformed standard 14-uppercase-letter key blocks, nonnumeric,
  infinite, zero and negative **present** molecular masses.
- Computes each shared key's minimum and maximum reported neutral masses.
  A spread exceeding **the greater of 0.01 Da or 20 ppm relative to the
  minimum mass** fails the build. These are conservative engineering
  tolerances for inconsistently annotated source identities, not an
  instrument-specific measurement uncertainty model.
- Reports only aggregate completeness counts, duplicate-key counts and
  maximum within-key mass spread. It does not print licensed SMILES or
  row-level source identifiers.

InChIKey14 collapses stereochemical and isotopic distinctions. Therefore a
mass conflict may be a legitimate isotopologue or an incorrect annotation;
**it must be scientifically reconciled**, not silently discarded or forced
onto an arbitrary first record. This gate does not recompute SMILES from
InChI or guarantee that a remaining key is chemically correct.

## Exclusive staged Parquet publication

The candidate builder now writes to a uniquely named sibling staging file.
The Parquet footer is inspected for the expected number of rows, and only
then is the file published via a same-filesystem **exclusive hard link**.
A previously existing candidate pool cannot be overwritten, even if it
appears during the build. Staging files are cleaned on failures. Filesystems
without hard-link support fail closed rather than falling back to unsafe
replacement. This is not a distributed transaction or a guarantee of
durable writes after sudden power loss.

The output `artifacts/pool{CASMI_POOL_TAG}.parquet` remains immutable.
Create a newly tagged research variant rather than modifying the frozen
baseline; regenerate fingerprints and dependent candidate features under
that tag and record a source SHA-256 manifest.

## Verification

Synthetic regression tests cover same-key roundoff, major mass conflicts,
missing source rows, invalid keys, empty or malformed SMILES, bad masses,
existing-output protection, staged-file cleanup, concurrent output creation
and new-target publication.

**Local isolated NumPy/pandas/standard-library smoke checks:** 5 passed.
The environment lacked `pyarrow`, so a sixth full Parquet footer/readback
smoke check could not run here. The complete repository scientific suite
and the actual multi-provider CASMI dataset build were **not run**.
Per user instruction no CI was triggered; branch and merge commits use
`[skip ci]`.

## Scientific limits

The historical public MRR@25 **0.422**, recorded **2026-10-08**, remains
unchanged. No Kaggle/BigQuery data was retrieved, no model was retrained and
no leaderboard prediction was submitted. The magnitude or prevalence of
actual source mass conflicts remains unknown until authorized provider
data is available. This work is public research only, not productization.
