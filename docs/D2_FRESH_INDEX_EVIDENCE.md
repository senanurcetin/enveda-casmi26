# D2 analog-reference index: cache-free scientific evidence

**Enveda CASMI26 | Research-only correction | 2026-10-19**

## Finding

The historical D2 analog-reference index (`analog_index.npz`) includes a
sparse fragment/neutral-loss matrix and the selected library row IDs.
Its **legacy** cache format does not carry content-bound identity for
the cleaned MS/MS source arrays, reference-selection configuration or
matrix semantics. It can therefore be reused by a new experiment even
when its underlying source store or index-building parameters have
changed. The opt-in source-balanced cache has stronger selection/file
stamp checks, but its stamp is not a cryptographic assertion of the
actual sparse matrix derivation.

The earlier decision gate checked matching source SHA-256 manifests and
Morgan fingerprints, but a D2 run using the historical unverified cache
could still receive an **ACCEPT/REJECT** model-comparison verdict.

## New validation policy

The frozen legacy workflow is retained for exploratory/replication use.
No old cache is overwritten or deleted.

For a new **decision-grade D2** experiment, execute the scientific
validation with all of these options:

```bash
python scripts/run_validation.py \
  --analog \
  --require-fingerprint-provenance \
  --fresh-analog-index \
  --source-manifest runs/verified-input-snapshot.json \
  --library enveda-np-examples \
  --max-queries 250 --max-spectra 5 --query-seed 0 \
  --tag fresh-d2-reference
```

- `--fresh-analog-index` constructs the D2 analog reference matrix
  **from currently loaded SpectrumStore files**, bypassing both legacy
  and source-balanced `*.npz` cache reads/writes. This can be **slow**
  and resource-intensive on millions of spectra. It is an explicit
  opt-in rather than a silent change in baseline behavior.
- The run records `analog_index_provenance`: index schema/version,
  `REBUILT_WITHOUT_CACHE` status, selected-reference count,
  sparse nonzero count, selected policy, and a SHA-256 digest over
  selected row IDs, matrix structure/content and configuration.
- The comparator requires this evidence for **every D2 run** in addition
  to signed matching query spectra, complete matching local input
  manifests and independently verified pool-bound Morgan vectors.
- Historical D2 reports lacking this new index evidence remain
  **exploratory** even if their older source/fingerprint evidence is
  otherwise valid. Malformed/misrepresented cache-free status or matrix
  hashes are rejected.
- A D1-only comparison does not require analog-index evidence.
  A D1-vs-D2 experiment can still be conclusive when the D2 run is
  verified and both share a complete identical source inventory.
  Two freshly built D2 indexes using different intentionally selected
  reference policies can be compared without requiring the two
  *model-specific index* hashes to match.

## Evidence limits

This hash identifies the **matrix actually created by the code**, not
the correctness of source chemistry, machine learning scores, or
upstream Kaggle/BigQuery data. `run_validation.py` must verify its
local source artifact SHA-256 inventory before building the index.

Neither a SHA-256 checksum nor a self-reported experiment JSON is a
cryptographic signature or independent third-party audit. Learned
model weights, mutable metric-key caches, instrument metadata validity,
software environment and held-out sample independence still need
separate controls.

The change does not alter the historical leaderboard score or
replace the public baseline. No licensed dataset, Kaggle Notebook,
BigQuery warehouse job or new scientific model run was executed here.
Recorded historical public MRR@25 **0.422 (2026-10-08)** is unchanged.

## Synthetic regressions

- A cache-free rebuild generates a deterministic index content digest.
- Rebuilt indexes with different reference policies may yield different
  digests.
- Cached, previously saved or loaded indexes cannot claim the
  cache-free provenance status.
- Comparison rejects false cache-free declarations, broken digest/count
  fields and policy inconsistencies.
- Historical D2 reports without index attestations are exploratory.

**User-requested CI policy:** no GitHub Actions runs were triggered;
commits use `[skip ci]`. The full repository scientific suite and
dataset-backed validation remain pending.
