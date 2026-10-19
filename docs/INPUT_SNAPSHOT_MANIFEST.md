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

## Complete-source comparison gate (2026-10-19)

The first version of the experiment comparator verified *query* identities and
selected MS/MS spectra but allowed source-pool artifacts to differ silently.
Now the two levels of evidence are tracked separately:

1. **Query-paired:** Molecule IDs, actual query spectra and C1/C2/C3
   reference exclusions match. Old runs with these signatures remain
   scientifically useful as **exploratory comparisons**.
2. **Decision-eligible:** In addition, **both** reports declare an identical
   verified SHA-256 input manifest with documented coverage of the exact files
   consumed by their pipeline. Only this level can emit a statistical
   model-comparison `ACCEPT` or `REJECT`.

When `--source-manifest` is supplied, the validation runner now checks
**required file coverage** in addition to byte hashes. D1 requires
`store/meta.parquet`, `store/offsets.i64`, `store/mz.f32`,
`store/it.f32` and the active `pool{CASMI_POOL_TAG}.parquet`. D2 also
requires `pool_fp2048{CASMI_POOL_TAG}.npy` and
`pool_fp_keys{CASMI_POOL_TAG}.parquet`. Extra files are allowed. A manifest
that omits any required file **fails before** the expensive scoring run.

The report records `source_manifest_schema`,
`source_manifest_pool_tag`, `source_manifest_coverage` and the
`source_manifest_required_inputs` list, together with the existing SHA-256
and checked file count. `scripts/compare_runs.py` checks that both reports
have matching hashes and valid coverage for their own active channels;
incompatible hashes/missing manifests fail closed, while legacy reports
with no manifest or incomplete coverage produce only `UNVERIFIED`.
This is a deliberate research-policy tightening, not an improved model score.

**Important:** This attests to a minimum *input set*, not to the algorithm's
Git commit, derived D2 cache, mutable metric-key memoization or learned model
weights. Explicitly inventory additional relevant artifacts and preserve
environment/dependency provenance when making a claim about those components.
Even complete local SHA-256 manifests do not authenticate the upstream
Kaggle/BigQuery source.

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
