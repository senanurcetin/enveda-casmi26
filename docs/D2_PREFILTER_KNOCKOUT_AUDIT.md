# D2 analog top-k knockout — scientific audit (2026-10-18)

**Scope:** existing non-commercial CASMI26 research pipeline. No deployed application, model retraining, competition submission, or inference score claim.

## Root cause

The previous implementation of `AnalogIndex.search` first masked incompatible precursor mass and adduct, then selected `prefilter_top` candidate rows based on approximate spectral similarity. **Only after the top-k operation** did it remove reference rows listed in `excluded_rows`.

For C1/C2/C3 structure knockout experiments, this ordering is incorrect. If the best approximate hits belong to excluded structures, they take slots away from legitimate independent analogs, which are never advanced to exact entropy similarity. A query with only excluded references in its first shortlist may report zero D2 analogs despite having eligible reference spectra with slightly lower prefilter scores.

This is a **false-negative / evaluation selection bias**, not direct reference leakage in the final D2 similarity output. The C1 correction in PR #5 addressed its *reference-exclusion population*; this sprint addresses *when* that exclusion is applied in D2.

## Implemented correction

- `casmi.channels.prefilter.eligible_top_k` creates a shortlist exclusively from positive, finite, eligible scores. Ties are deterministic by original reference index and the input scores/masks are not modified.
- `AnalogIndex.search` evaluates mass/adduct restrictions and structure-exclusion rows **before** the approximate top-k step. References belonging to held-out query structures cannot consume shortlist slots.
- An empty index, all-excluded search and non-finite scores return an empty shortlist without invalid indexing.
- Tests reproduce the old failure with synthetic high-scoring self-reference rows: valid analogs are now selected even when the excluded rows have better approximate scores.
- GitHub Actions includes lightweight NumPy algorithm checks and a direct `AnalogIndex.search` regression suite under a pinned scientific environment.

## Metadata-only index-retention audit

The scientific pipeline can now quantify the remaining offline selection bias
without reading MS/MS spectra, private labels or Kaggle answers:

```sh
python scripts/audit_analog_retention.py \
  --meta artifacts/store/meta.parquet \
  --query-library enveda-np-examples \
  --per-structure-adduct 3 \
  --json-out runs/analog_retention_audit.json
```

The output estimates how many **standard-key/adduct** pairs have at least one
qualified reference in the query library and in another library, but lose
every cross-library alternative during index subsampling. This reports
eligible spectra, indexed spectra, pairs lost, the loss fraction, unique
structures affected and a SHA-256 of the metadata file. These are **data
quality counts**, not per-molecule recall or MRR. The audit is currently
code- and fixture-tested; no real-store output has been collected.

## Known remaining limits

1. The on-disk index is itself constructed by retaining the best `per_structure_adduct` reference rows **before a validation knockout exists**. If those indexed rows all belong to the withheld library, cross-library references omitted during index construction cannot be recovered by a correct runtime prefilter. Fixing this would require an explicitly versioned/benchmarked index-construction experiment and cache invalidation.
2. Candidate pool, learned rankers and D3/D4 predictions have not been revalidated in this change.
3. C1/C2/C3 and timsTOF domain shift require an actual dataset-backed rerun after these integrity fixes. The recorded historical 2026-10-08 public MRR@25 of 0.422 remains unchanged.

## Acceptance criteria

```sh
python -m unittest discover -s tests -p '*_integration.py' -v
python -m pytest tests/test_analog_index_knockout.py -q
```

The first command requires NumPy; the second needs NumPy, SciPy, pandas, Numba and pytest. Both are exercised in the NEXUS GitHub Actions research workflow.

**Decision:** CI-tested algorithmic correction may merge; experimental predictive improvement is **NOT PROVEN** until an independent molecule-level validation rerun is available.
