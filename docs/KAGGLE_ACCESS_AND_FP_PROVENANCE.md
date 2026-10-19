# Kaggle CASMI26 access and fingerprint provenance

**Research note (2026-10-19).** This remains public, non-commercial
scientific work. No source datasets, credentials, Kaggle submissions or
private labels belong in this GitHub repository.

## Official competition context

- Competition: `enveda-CASMI26-molecule-id-mass-spectra`
- Website: https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra
- Source files: `train.parquet`, `test.parquet`,
  `sample_submission.csv`. The official dataset is about **3.04 GB**.
- The published `test.parquet` may contain examples drawn from training;
  the official hidden evaluation file replaces it during the Kaggle
  Notebook rerun. Do **not** treat public example overlap as independent
  test-set evidence.
- Competition submissions must use Kaggle Notebooks with internet disabled
  during execution, respect the official runtime rules, and generate
  `submission.csv`.
- Official source data is designated CC BY-NC 4.0 and competition file
  access requires acceptance of the competition rules in the Kaggle account.

Always check the current Rules and Data pages before a new experiment.

## Credentials: no token in source code

Kaggle's CLI supports the current `KAGGLE_API_TOKEN` environment variable.
Set it **locally**, using a secret manager or an untracked process
environment. Do not paste tokens into the repository, command-line
arguments, commit messages, notebook outputs, public issues or shared logs.

```bash
python -m pip install --upgrade kaggle

# Configure KAGGLE_API_TOKEN in your local secure environment.
# Do not paste its value into public scripts or terminal recordings.
python scripts/check_kaggle_access.py

# Optional read-only authenticated file listing after accepting rules:
python scripts/check_kaggle_access.py --check-remote
```

The preflight prints **only sanitized statuses**. It never downloads
data, accepts Kaggle rules, generates submissions or uploads predictions.
If access is unconfirmed, check that the token is active, the Kaggle CLI
is installed, the competition rules have been accepted, and the account
has permission. Remote access is not claimed until a real authenticated
listing succeeds in the user's environment.

The credential supplied in a chat should be **rotated** once work has
been moved to a secure local token store. The assistant did not store
the token in files, commits, or environment variables.

## Critical source correction: COCONUT parent -> fingerprint mismatch

[The COCONUT neutral-parent identity fix](COCONUT_PARENT_IDENTITY_AUDIT.md)
changes some candidate SMILES, standard InChIKey14 values, masses and
formulas. Existing Morgan fingerprints **cannot be assumed compatible**
with a pool rebuilt from corrected parent identities.

The new `fingerprint_provenance.py` binds a fingerprint array to the
*ordered* candidate (pool key, SMILES, mass) values and exact packed
fingerprint bytes using SHA-256. It does not expose those candidate rows.

A safe new experiment creates a **separately tagged pool** and its
fingerprints; the build scripts refuse to overwrite frozen variants:

```bash
# In an authorized local scientific environment with source data available:
export CASMI_POOL_TAG=_parent_v2
python scripts/build_pool.py
python scripts/build_fingerprints.py
```

The generator now emits:
- `pool_parent_v2.parquet`
- `pool_fp2048_parent_v2.npy`
- `pool_fp_keys_parent_v2.parquet`
- `pool_fp_provenance_parent_v2.json`

The sidecar specifies schema version, row count, ordered pool-row digest
and packed fingerprint-byte digest. It also checks the exact row order.

When running D2 knockout validation, optionally **require** the
sidecar and fail closed on stale pool SMILES, mass, ordering or fp bytes:

```bash
export CASMI_POOL_TAG=_parent_v2

python scripts/run_validation.py --analog \
  --require-fingerprint-provenance \
  --library enveda-np-examples --query-seed 0 \
  --max-queries 250 --max-spectra 5 \
  --tag parent-v2-d2
```

For complete experimental evidence, also pass an input SHA-256
`--source-manifest` covering the store, candidate pool and exact fp
artifact files. Keep comparisons to the same query-signature version,
candidate-pool policy and licensed source snapshot.

**Legacy compatibility:** Historical D2 runs lacking the new sidecar
are still executable without `--require-fingerprint-provenance`;
they produce a visible warning and do not acquire retroactive source
verification. If a sidecar exists, an incompatible array or pool **always
fails closed**, even without the optional flag.

## Verification and limitations

- New offline tests use **synthetic candidate rows** and false test
  credentials only. They check identity, mass and fingerprint drift,
  byte corruption, row ordering and read-only Kaggle preflight behavior.
- The Kaggle account was not authenticated here. No official 3 GB
  dataset was downloaded, no submission was made, and no new leaderboard
  result is claimed.
- A digest establishes *local artifact consistency*, not correctness of
  external sources, RDKit chemistry or leaderboard predictions. A
  separately tagged pool still requires full feature/regressor
  rebuilding and held-out C1/C2/C3 validation.
- The historical public MRR@25 `0.422` recorded on 2026-10-08
  remains unchanged. No CI workflows were invoked in this work cycle.
