# Morgan fingerprint zero-vector scientific audit

**Enveda CASMI26 | 2026-10-19 | research-only**

## Finding

The prior fingerprint generation loop treated an invalid RDKit SMILES as an
empty Morgan vector by leaving a preallocated row of zeros and continuing.
Those vectors could then be written alongside otherwise valid candidates
and included in a versioned SHA-256 fingerprint provenance sidecar.

A hash binds a vector to a source snapshot, but it **does not establish
that the fingerprint is chemically meaningful**. Consequently, a
cryptographically consistent artifact might still contain chemically
invalid zero-vector candidate rows.

## Corrections

- Fingerprint generation now fails before artifact publication when
  RDKit cannot parse a candidate, the parsed structure has zero heavy
  atoms, or a generated Morgan fingerprint has no set bits.
- Errors report the **candidate row number only**. Source SMILES and
  other potentially licensed molecular data are not copied into logs.
- The generator uses `CandidatePool` itself for mass validation,
  stable sorting and key deduplication. Previously it implemented
  those steps a second time, risking row-order drift after revisions
  to the inference pool.
- The historical `pool_fp2048*.npy`, `pool_fp_keys*.parquet` and
  sidecar files are never overwritten. An invalid candidate aborts
  before any of these output artifacts are created.

## Auditing existing frozen fingerprint artifacts

A separate **read-only**, chunked audit can inspect historical vectors
without regenerating or rewriting them.

```sh
# Point at an authorized snapshot (do not overwrite the original).
export CASMI_ARTIFACTS=/approved/casmi/artifacts
export CASMI_POOL_TAG=_parent_v2

python scripts/audit_pool_fingerprints.py \
  --batch-rows 65536 \
  --output runs/parent_v2_fingerprint_quality.json
```

The output reports the candidate count, number and fraction of all-zero
packed Morgan vectors, row alignment and whether a matching content
manifest was verified. The audit exits with a nonzero status if any
zero-vector candidate is found. It also rejects malformed vector arrays,
mismatched key order, and empty candidate pools.

The audit reads fingerprint arrays using NumPy memory mapping, limiting
working batch size. Manifest checking additionally reads the fingerprint
bytes to compute SHA-256. It reports **aggregate numbers only**, not
molecule rows, SMILES, sources or private competition labels.

**Interpretation:** `PASS_NONZERO_FINGERPRINTS` says that every inspected
packed vector has at least one set bit and the key order matches.
It does not prove that those bits were generated from the correct
structure; this requires the independently generated source-binding
fingerprint manifest. Even both checks do not independently authenticate
Kaggle's original files, source licenses or the true test distribution.

## Research evaluation and operational constraints

A frozen historical pool should never be rebuilt in place. Instead:

1. Select an approved, **new** `CASMI_POOL_TAG`.
2. Rebuild the corresponding pool and fingerprints together.
3. Audit the resulting fingerprint files and verify their provenance.
4. Record a full source snapshot and evaluate C1/C2/C3 on the
   identical independently selected query population.
5. Train a new ranker if fingerprint-derived features have changed;
   never compare mismatched historical model artifacts.

Tests cover valid Morgan vectors, invalid SMILES, wildcard/empty molecule
inputs, failure before publication, source-pool mass validation, vector
shape, row alignment and aggregate zero-row counting.

**Verification limitation:** Focused local RDKit smoke checks passed.
No GitHub CI was invoked per the user's instruction. The full scientific
test suite and real licensed dataset fingerprint rebuild remain unrun.

The historical public MRR@25 **0.422 (2026-10-08)** is unchanged.
No productization, deployment, Kaggle submission, BigQuery query or new
leaderboard claim is part of this correction.
