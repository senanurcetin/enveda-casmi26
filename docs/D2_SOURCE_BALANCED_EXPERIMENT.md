# D2 source-balanced index — controlled research experiment

**Scope:** CASMI26 scientific validation only. No SaaS, production service,
paid deployment, Kaggle submission or unverified leaderboard claim.

## Problem and hypothesis

The old D2 reference index takes up to three top spectra per
`(standard_2D_key, adduct)`, preferring timsTOF spectra regardless of source
library. During C1 cross-library knockout, all three retained spectra can
belong to the withheld query source even when qualified cross-library
reference spectra existed in the underlying store. PR #6 fixed *query-time*
knockout ordering; it could not restore references never indexed.

**Testable hypothesis:** At the same per-key/adduct reference budget,
source-aware selection will improve the **cross-library reference retention
rate** without reducing the quality of the top-ranked reference. Whether
this helps MRR@25 on any novelty class remains unknown.

## Two explicitly selectable policies

- `legacy` **(default, historical baseline):** all three slots selected by
  global timsTOF/peak-count priority with stable ties, matching the previous
  `AnalogIndex._select` behavior.
- `source_balanced` **(research opt-in):** retain the top-scoring reference,
  reserve the highest-priority independent non-timsTOF source when the first
  is timsTOF and the budget allows, then choose other distinct source
  libraries before filling leftover slots by legacy quality priority.

Neither policy increases the per-key/adduct cap. When a pair has only one
qualified source library, no cross-library recovery is possible.

## Cache and reproducibility

A `source_balanced` index uses a separate cache file, for example
`artifacts/analog_index_source_balanced_v1.npz`, never the historical
`artifacts/analog_index.npz`. The new cache validates its selection-policy
format, selected row IDs, index configuration and selected-row metadata
fingerprint; it also incorporates the raw peak file sizes and modification
times when accessible. Corrupt/mismatched nonlegacy caches are rebuilt.

The file stamp is **not** a cryptographic full-spectrum content hash.
Independent data provenance should additionally record input snapshots and
checksums. The unchanged legacy cache is *not* upgraded retroactively; do
not compare old provenance to new outputs without explicitly recording the
policy and source snapshot.

## Controlled comparison without leaderboard probing

First audit how many (standard-key, adduct) pairs can retain an independent
cross-source reference using an approved, locally available metadata file:

```bash
python scripts/audit_analog_retention.py \
  --meta artifacts/store/meta.parquet \
  --query-library enveda-np-examples \
  --per-structure-adduct 3 \
  --selection-policy legacy --compare-policies \
  --json-out runs/d2_source_retention.json
```

Then run the same frozen query population, spectra and seed twice:

```bash
python scripts/run_validation.py \
  --analog --analog-selection-policy legacy \
  --library enveda-np-examples --max-queries 250 \
  --max-spectra 5 --query-seed 0 --tag d2-legacy

python scripts/run_validation.py \
  --analog --analog-selection-policy source_balanced \
  --library enveda-np-examples --max-queries 250 \
  --max-spectra 5 --query-seed 0 --tag d2-balanced

python scripts/compare_runs.py runs/EXACT_LEGACY_REPORT.json \
  runs/EXACT_BALANCED_REPORT.json
```

Replace report filename placeholders with the actual generated files. Do
**not** compare different `query_signature_version` values or source
snapshots, and do not assume the change must help every novelty mixture.

For research-only training feature tables, isolate variants under a suffix:

```bash
python scripts/build_rank_data.py --set vin \
  --analog-selection-policy source_balanced --suffix _d2balanced
```

Training a separate ranker would need a corresponding new feature set,
identical molecule-disjoint splits and an untouched test. **Do not use a
historically trained ranker to claim a controlled D2 benefit.**

## Success and failure criteria

1. The synthetic selector and D2-index regressions pass and prove the old
   selection is unchanged for the default policy.
2. The source-balanced index retains at least one cross-source row when a
   suitable source exists under the configured cap and reservoir conditions.
   No claim is made for all possible source combinations with more source
   types than the cap.
3. The metadata-only audit explicitly lists legacy and balanced lost-pair
   counts; no private spectra or labels appear in committed artifacts.
4. **Algorithm acceptance** is independent of leaderboard performance.
   **Scientific effectiveness** requires new, signed, query-aligned
   validation reports and simultaneous confidence bands; otherwise
   verdict is `NOT PROVEN`.

## Status (2026-10-18)

- Opt-in policy, index selection, cache isolation, synthetic tests and
  metadata-audit comparison: **implemented**.
- Live Kaggle/BigQuery evidence, full MS/MS rerun, learned-ranker retraining,
  independent model performance evaluation: **not performed**.
- Recorded public MRR@25 `0.422` from 2026-10-08: **historical, unchanged**.
