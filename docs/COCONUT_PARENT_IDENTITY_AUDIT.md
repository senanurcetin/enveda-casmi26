# COCONUT neutral-parent identity integrity audit

**Date:** 2026-10-19  
**Scope:** Public, non-commercial Enveda CASMI26 scientific research only.

## Scientific issue

The original COCONUT ingestion procedure removed salts/counterions from
multi-component `canonical_smiles` by retaining the largest heavy-atom
fragment. It then recomputed exact mass from that fragment, **but retained
the original multi-component `standard_inchi_key` as `pool_key`**.

As a result, one candidate could contain:

- SMILES describing the neutral parent,
- a mass describing that parent,
- an InChIKey14 describing the *original salt or multi-component mixture*, and
- a source molecular formula describing the original mixture.

This inconsistency matters for structure deduplication, matching reference
spectra to candidate structures, and C1/C2/C3 novelty/knockout logic.

## Correction

For every multi-component COCONUT record eligible for neutral-source
ingestion, `neutral_parent_identity` now:

1. Parses the SMILES with RDKit and chooses its largest heavy-atom fragment.
2. Rejects malformed inputs, fragments with zero heavy atoms and **charged
   parent fragments**; no placeholder molecule is substituted.
3. Derives canonical parent SMILES, **standard non-tautomer InChIKey14**,
   exact monoisotopic mass, and molecular formula from the *same* RDKit
   molecular graph.
4. Discards invalid parents before mass-window filtering and key
   deduplication.
5. Prefers directly supplied single-component neutral parents over
   salt-derived duplicates of that same parent, preserving annotations
   from the unmodified parent record.

The use of *standard* InChIKey14 here is intentional. The candidate pool
uses standard structural identity, whereas later leaderboard-metric
equivalence is handled separately by tautomer-canonical metric keys. We
must **not** silently substitute the metric key for the pool's standard key.

## Synthetic regression and local focused checks

Tests include a synthetic aspirin + sodium/chloride mixture and a directly
represented neutral aspirin molecule. The following properties are checked:

- The compound/salt InChIKey14 differs from the neutral parent key.
- After stripping salts, parent SMILES, InChIKey14, formula `C9H8O4`
  and exact mass approximately 180.0422587 describe the same molecule.
- Both salt and parent records deduplicate under the same parent key.
- Direct neutral-source annotations win when an identical salt-derived
  parent is also present.
- Invalid SMILES, all-inorganic salt fragments and charged organic
  parents are excluded instead of mapped to an unrelated fallback.

A focused local RDKit/pandas smoke check passed for the parent identity,
mass, formula and key deduplication. **No CI was invoked**, as requested.
The full source-data ingestion and full repository test suite have not
been run for this change.

## Artifact, data and metric boundaries

This change **does not modify existing `pool.parquet`**, fingerprint
arrays or historical training artifacts in place. A researcher who has
authorized COCONUT source files may build a separately versioned research
pool, then regenerate the *matching* fingerprints and candidate features
under that experiment's own provenance manifest.

Do not combine a newly normalized candidate pool with a previously cached
`pool_fp2048.npy`/key mapping or a historical ranker without an
independent compatibility check. Doing so could silently misalign candidate
rows and structure identities.

Historical public **MRR@25 0.422**, recorded 2026-10-08, remains unchanged.
No new Kaggle leaderboard result, BigQuery data verification, independent
holdout evaluation, product, API or monetization claim is made.
