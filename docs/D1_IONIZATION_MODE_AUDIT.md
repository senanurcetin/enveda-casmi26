# D1 ionization-mode alignment — scientific reliability note

**Research-only correction (2026-10-19).** This change addresses direct-library
MS/MS comparison in mixed-ionization acquisitions. It introduces no production
features, retraining, external data access or leaderboard submissions.

## Issue and mechanism

Previously `LibraryIndex.score` selected the reference library spectra for
the ionization mode of **the first query spectrum**, then compared every
query spectrum against those reference spectra. If a query contained both
positive- and negative-mode acquisitions, one mode could receive a spuriously
high spectral similarity against references collected in a different mode,
while actual same-mode references were never considered.

This is a **scientific validity and ranking concern**, not demonstrated
improvement in a measured dataset. The original query grouping generally
contains same-mode spectra; the correction explicitly addresses the
heterogeneous case without changing the intended same-mode reference logic.

## Corrected selection

- Partition query spectra by recorded `mode`.
- Apply D1 source exclusion, quality thresholds and `max_refs_per_candidate`
  **within each mode**, prioritizing adducts from queries in that mode.
- Compare each query spectrum only with matching-mode references.
- A query with no references in its mode contributes **zero** to the
  mean-of-best-per-spectrum similarity, rather than being compared with an
  incompatible reference.
- Count distinct reference rows used across all modes. Empty inputs or
  candidates with no references return `(0.0, 0.0, 0)`.

The scoring function preserves the conventional maximum and mean of
per-spectrum similarities for ordinary single-mode queries. For mixed-mode
queries, the final result now respects instrumental mode separation.

## Evidence and remaining limits

Synthetic tests cover mixed positive/negative spectra, reversed query order,
reference caps, knockout exclusions, fallback reference counts, unmatched
mode, empty queries and missing candidate keys.

A **local isolated functional check passed 8 assertions**. The user's
instruction to avoid CI was respected: GitHub Actions was not invoked, and
commits use `[skip ci]`. This is not a claim that the full project's
chemistry suite or original CASMI dataset was run.

The historical 2026-10-08 public MRR@25 of **0.422 remains unchanged**.
Before declaring any predictive benefit, use an approved local source
snapshot, fixed molecule-disjoint evaluation and signed C1/C2/C3 reruns.
The research repository is public; restricted source data must remain private.
