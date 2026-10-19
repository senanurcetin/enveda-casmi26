# D2 fingerprint evidence and candidate mass integrity

**2026-10-19 | Scientific research only, no CI requested**

## Finding 1: source-manifest equivalence was not enough for D2 model decisions

`run_validation.py` can now validate a Morgan fingerprint sidecar that
binds each packed bit-vector to the ordered candidate pool's exact SMILES and
mass. Yet `align_runs` previously treated a matching general source-file
manifest as sufficient to issue a statistical **ACCEPT/REJECT** verdict, even
when a D2 run did not have a verified fingerprint sidecar.

A byte-identical set of files does not by itself establish that fingerprint
rows were generated from the correct SMILES/mass identities. This is
especially relevant following the COCONUT parent/salt identity corrections.

**Fix:** `fingerprint_evidence.py` adds a standard-library validation gate
for D2 report evidence. A decisive paired model comparison now requires:
- The existing molecule, query-spectra and reference-exclusion signatures;
- Two matching complete local input-source manifests;
- For **each D2 run**, a verified fingerprint provenance record with the
  supported schema, positive fingerprint-row count and valid SHA-256 digests
  for the ordered pool rows and fingerprint bytes;
- If both runs use D2, matching fingerprint identities across those runs.

Historical D2 reports missing a verified sidecar remain **exploratory**,
without ACCEPT/REJECT. Conflicting D2 fingerprints or malformed positive
claims are rejected outright. D1-only comparisons do not require D2
fingerprints. A D1–D2 ablation is eligible if both runs use the same
**superset input snapshot**, and the D2 run has its independent verified
fingerprint record. This preserves scientifically useful model ablations.

These checks confirm documented *verification evidence* in run reports.
They are not tamper-proof digital signatures and cannot authenticate the
underlying Kaggle publisher or correctness of the external dataset.

## Additional D2 index requirement (2026-10-19)

A content-bound candidate fingerprint sidecar does **not** prove that an
`analog_index.npz` reference spectrum matrix came from the currently
verified source store. The comparator now also requires
`analog_index_provenance` from a cache-free rebuild whenever either
side uses D2. Use `run_validation.py --analog --fresh-analog-index`
together with the fingerprint sidecar and local input manifest.
Older cached-index reports remain **exploratory**. Full details:
[D2 cache-free reference provenance](D2_FRESH_INDEX_EVIDENCE.md).

## Finding 2: invalid candidate masses could enter ppm ranking

`CandidatePool` used to drop null values and sort by mass, but could accept
`mass <= 0` or infinity. Its `window()` divides by candidate mass to form
candidate ppm, which could create infinite or undefined values. Invalid
query masses or tolerances could also yield undefined or empty result
windows without an explanatory error.

**Fix:** Candidate pools now reject finite-but-nonpositive and non-finite
mass entries, and reject nonnumeric present masses with a clear error.
Missing mass records are still dropped per the existing data preparation
behavior. The mass-window search requires positive finite query mass and
finite nonnegative ppm tolerance, while correctly handling array/set
exclusion collections, including empty exclusions.

## Verification and status

- Added standard-library D2 report evidence regressions: missing sidecars,
  changed fingerprint bytes or pool digest, malformed verified claims,
  baseline D1–D2 ablations and original source-snapshot mismatches.
- Added synthetic NumPy/Pandas candidate mass and ppm query-window tests.
- Updated the existing D2 signed-report test fixture to provide the
  fingerprint evidence now required for a scientific decision.
- The user requested **no CI**; commits use `[skip ci]`. No GitHub Actions
  workflow was triggered. New checks still require execution in the full
  scientific dependency environment before claiming suite-wide success.
- No data download, Kaggle API execution, BigQuery job, model retraining or
  new public leaderboard submission took place.

The historical **0.422 MRR@25 on 2026-10-08** is not a new result. Research
models and immutable historical artifacts remain unchanged. No SaaS,
commercialization, inference product or customer-facing deployment is part
of this work.

## Reference commands for an authorized local environment

```bash
python -m unittest discover -s tests -p '*_stdlib.py' -v
python -m pytest tests/test_candidate_pool_mass_validation.py -q

# Scientific inference: use a complete matching source snapshot and
# explicitly authenticated fingerprint provenance when D2 is active.
python scripts/run_validation.py --analog \
  --require-fingerprint-provenance \
  --source-manifest runs/ACTUAL_SOURCE_SNAPSHOT.json \
  --max-queries 250 --max-spectra 5 --query-seed 0 \
  --tag verified-d2-reference
```
