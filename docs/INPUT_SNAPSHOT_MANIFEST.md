# Input snapshot integrity for CASMI26 research (2026-10-19)

The CASMI validation reports already identify their query molecules, selected
spectra and C1/C2/C3 reference exclusions. These fingerprints **do not identify
all underlying input artifacts**, such as candidate-pool files or the memmapped
reference spectrum store. The local input snapshot manifest introduced here
makes permitted dataset/feature artifact provenance independently auditable.

## Explicit, read-only SHA-256 inventory

The manifest writer streams selected files and records their **relative**
paths, exact byte sizes, SHA-256 digests and a deterministic composite file-list
digest. It does not copy source files, store their contents or expose absolute
machine paths, credentials, Kaggle test answers or BigQuery records.

Use **only** the locally approved, available input files, e.g. in the standard
artifact directory:

```bash
python scripts/manifest_artifacts.py \
  --root /YOUR/LOCAL/ARTIFACT_DIR \
  --file store/meta.parquet \
  --file store/offsets.i64 \
  --file store/mz.f32 \
  --file store/it.f32 \
  --file pool.parquet \
  --file pool_fp2048.npy \
  --file pool_fp_keys.parquet \
  --output runs/casmi-source-snapshot.json
```

The SHA-256 inventory is **not zero-cost**: hashing every multi-GB memmap
requires reading every byte. Run only when there is sufficient disk I/O and
user-authorized compute. Include the correct tagged pool and fingerprint
artifacts when `CASMI_POOL_TAG` is set. Optional trained-model weight files
should be included only if relevant, and their presence cannot be inferred
from this manifest unless explicitly listed.

Verify independently before expensive model scoring:

```bash
python scripts/manifest_artifacts.py \
  --root /YOUR/LOCAL/ARTIFACT_DIR \
  --verify runs/casmi-source-snapshot.json

python scripts/run_validation.py --analog \
  --source-manifest runs/casmi-source-snapshot.json \
  --library enveda-np-examples --max-queries 250 \
  --max-spectra 5 --query-seed 0 --tag baseline-verified
```

The validation runner now verifies all listed local artifacts **before**
loading the spectrum store; any missing, changed, corrupted, path-escaping or
invalid manifest aborts the validation. The report records
`source_manifest_sha256` and `source_manifest_file_count`, or `null`
when the run was intentionally executed without a manifest. Blank/shuffled
diagnostic pairings reject unequal snapshot digests when at least one side
declares a snapshot.

## Evidence limits

- A matching file hash proves **local byte identity**, not that the dataset
  came from Kaggle's intended immutable version or a correctly configured
  BigQuery table. External provenance still requires provider-side validation.
- Inclusion is explicit; a manifest with one file **cannot** attest to the
  integrity of omitted files. Do not describe an incomplete inventory as a
  complete dataset snapshot.
- Manifests are **not tamper-proof digital signatures**, and full SHA-256 does
  not replace query identity, structure-disjoint validation or leakage checks.
- The original competition public MRR@25 `0.422` is dated 2026-10-08. No
  new CASMI, Kaggle or warehouse workload was run to implement this utility.
- Keep generated inventory files in `runs/` or another private location.
  Do not commit restricted metadata or raw datasets into the public repo.
- No productization is authorized or included.
