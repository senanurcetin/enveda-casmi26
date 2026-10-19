# Enveda CASMI26 — molecule-level diagnostic error slices

**Status (2026-10-19):** Research-only error-characterization code and
synthetic regression tests. No licensed CASMI data or Kaggle leaderboard
rerun has been performed. No commercialization or model improvement claimed.

## What the existing aggregate MRR hides

C1/C2/C3 MRR@25 summaries average over molecule inputs with different
fragment evidence. A higher mean score can conceal failure for particular
adducts, modes, precursor-mass windows or query-spectrum densities.
This addition makes those heterogeneities observable **without requiring
molecule labels in a published report**.

For each unperturbed input query, `run_validation.py` persists a
`query_context` object with four **predeclared** dimensions:

| Dimension | Categories | Data source |
| --- | --- | --- |
| Adduct | Exact metadata adduct string, or `mixed` | Selected original query spectra |
| Ionization mode | Exact reported mode, or `mixed` | Selected original query spectra |
| Query spectra | `1`, `2-3`, `4+` | Selected input spectrum count |
| Precursor m/z | `<300`, `300-600`, `>=600`, `unknown` | Median finite precursor m/z in Da |

The context is recorded **before** any blank-spectrum or
intensity-shuffle treatment. No held-out truth SMILES, original
spectrum arrays or row-level query IDs are emitted in the diagnostic
output. Source reports already contain query IDs for pairing.

## Offline descriptive comparison

Generate baseline and candidate validation reports with the same query
seed, query sampling, library, frozen source snapshot and per-query
signatures. Prefer `--source-manifest` on both runs to permit *complete*
source provenance verification as described in
[INPUT_SNAPSHOT_MANIFEST.md](INPUT_SNAPSHOT_MANIFEST.md).

```bash
python scripts/diagnose_validation_slices.py \
  runs/ACTUAL_BASELINE_REPORT.json runs/ACTUAL_EXPERIMENT_REPORT.json \
  --min-n 15 --output runs/error_slices.json
```

The tool:
1. Aligns the candidate's ranks and metadata by the actual molecule IDs,
   rejecting different query spectra/exclusion sets, duplicate IDs and
   manipulated-spectrum negative controls.
2. Requires exactly matching context metadata for each molecule.
3. Suppresses all MRR numbers for cells with fewer than `--min-n` molecules
   (default: **15**); output contains only the cell count and
   `SUPPRESSED_SMALL_CELL` status in those cases.
4. For sufficiently populated cells, reports C1/C2/C3 baseline MRR,
   candidate MRR and paired difference as **descriptive** values.
5. Distinguishes source-inventory completeness from query pairing in
   `source_evidence`, avoiding claims of fully verified source data
   where input manifests were missing.

A descriptive difference is **not statistically significant evidence**:
groups and novelty classes share molecules, multiple cells may be inspected
post-hoc, and no family-wise correction or independent holdout is performed.
The tool intentionally emits **no ACCEPT/REJECT field or confidence
intervals**. Use predeclared hypotheses and held-out replication before
calling any subgroup change a demonstrated improvement.

## Restrictions and next experiments

- The query library is selected per validation run. A single-library
  report cannot support source-library *comparisons*; perform separately
  matched runs per valid library if needed.
- Instrument family and collision-energy groupings are not yet emitted:
  check exact metadata schema, missing values, privacy and sampling counts
  before adding these slices.
- A slice's MRR only measures the given validation population and C1/C2/C3
  simulation protocol; it does not establish real Kaggle test distribution.
- Recorded historical public MRR@25 **0.422 (2026-10-08)** remains unchanged.
- The current CI exercises synthetic fixtures, not actual licensed spectra.

## Tests

```bash
python -m unittest discover -s tests -p '*_stdlib.py' -v
```

Coverage includes reordering candidate molecules, independent context
mutation, privacy-preserving aggregate output, small cells, invalid bins,
and rejecting changed query spectra or negative controls.
