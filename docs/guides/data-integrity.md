# Data integrity and provenance

These controls protect the association between molecular identities, candidate structures, fingerprints and spectrum arrays. They establish computational consistency; they do not independently authenticate chemical annotations or acquisition provenance. Dataset permissions and redistribution terms are described in the [data card](../../DATA_CARD.md).

## Spectrum stores

`SpectrumStore` reads a metadata table, int64 offsets and flat float32 m/z/intensity arrays. Before memory mapping, structural validation checks byte alignment, metadata row counts, zero initial offset, monotone nonnegative offsets, terminal peak count and agreement with every `n_peaks` value. Negative or out-of-range spectrum indices fail explicitly. Empty peak stores have defined behavior.

`build_store.py` writes into a new destination, validates the persisted layout and publishes only after successful construction. Existing stores are not overwritten. A default integrity audit checks layout; a full numeric scan is an explicit bounded-memory, full-I/O operation:

```bash
python scripts/audit_store_integrity.py --store /path/to/artifacts/store
python scripts/audit_store_integrity.py --store /path/to/artifacts/store --scan-peaks
```

Peak cleaning excludes nonfinite or nonpositive m/z and intensity values. Cleaner and similarity regressions cover invalid inputs, empty spectra and numerical boundaries.

## Candidate-source reconciliation

The candidate builder reconciles training structures, COCONUT and configured optional sources without mixing one provider's structure with another provider's mass. Conflicting standard identities are rejected by default. Completed studies explicitly opted into whole-identity mass-conflict quarantine and selected-structure usability screening, with counts retained in their reports. Quarantine never removes evaluation molecules from scoring denominators.

COCONUT parent identity is checked after the configured parent-structure transformation. These rules do not establish stereochemical identification or agreement with an external chemical authority.

## Fingerprints and generated candidates

Pool fingerprints must match candidate identities, row ordering, uint8 packing and the expected 2,048-bit Morgan representation. Builds reject invalid/no-heavy-atom structures and zero/error fingerprints. Provenance records bind source identities, configuration and stored vectors; a signed but all-zero vector is still invalid.

Generated overlays are registered atomically. Their identity binds standard InChI14 and the SHA-256 of canonical isomeric SMILES, because metric-equivalent tautomers can have different representation-sensitive fingerprints, fragments and NP-scores. Ambiguous reuse of a legacy generated key fails closed. The [fresh-study amendment](../../protocols/advanced-internal-v1-amendment.json) records the outcome-blind correction. The frozen retrospective notebook retains its actual earlier implementation.

## Input manifests

`manifest_artifacts.py` streams explicit files and records relative POSIX paths, byte sizes, per-file SHA-256 values and a composite inventory digest. Traversal, duplicates, missing files and changed contents are rejected. Hashing a large store reads every byte; it is an I/O operation rather than a metadata-only check.

```bash
python scripts/manifest_artifacts.py --root /path/to/artifacts \
  --file store/meta.parquet --file store/offsets.i64 \
  --file store/mz.f32 --file store/it.f32 \
  --file pool.parquet --file pool_fp2048.npy --file pool_fp_keys.parquet \
  --output runs/source-snapshot.json
python scripts/manifest_artifacts.py --root /path/to/artifacts --verify runs/source-snapshot.json
```

Use the actual tagged pool names and include every input required by the chosen experiment. General manifests do not imply that unlisted weights or optional feature files were verified. D2 comparisons also require [fresh-index evidence](evaluation.md).

## Data access and retained evidence

Licensed spectra, structures, trained weights and private per-query outputs stay outside Git. Credentials remain in the operator's existing local authentication mechanism. The [report directory](../../reports/README.md) contains aggregates, hashes, software versions and bounded publication readbacks. The executed [overlap audit](../studies/public-example-overlap.md) confirms numeric public-example overlap for one hashed snapshot; hidden evaluation data and live BigQuery provenance remain outside that finding.
