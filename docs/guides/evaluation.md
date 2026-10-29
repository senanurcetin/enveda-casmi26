# Evaluation methodology

The evaluation unit is a molecule. Multiple spectra and metric-equivalent molecular representations must be handled together. Full design details and outcome boundaries live with each [completed study](../README.md#completed-studies).

## Ranking and identity

MRR@25 uses the first tautomer-canonical InChIKey14 match within 25 ranked slots under the pinned RDKit metric implementation. Missing predictions contribute zero; invalid, duplicate or blank candidate slots must not be compressed into better ranks. Standard or tautomer identity alone does not establish scaffold independence.

Current numerical summaries and paired bootstrap accept nonempty one-dimensional integer rank arrays from 0 through 25; 0 denotes a miss. Nonfinite, floating, negative and out-of-range ranks fail closed. Bootstrap draws must be integers of at least 100 with a nonnegative integer seed. Tests verify valid fixed-seed results remain unchanged.

## Reference knockout

| Condition | Reference exclusion | Candidate handling |
|---|---|---|
| C1 | All own-library acquisitions for the query's metric-equivalent identity, including unsampled spectra | Retain ordinary mass-compatible candidates |
| C2 | All metric-equivalent truth references across libraries | Retain ordinary mass-compatible candidates |
| C3 | C2 reference exclusion | Exclude metric-equivalent truth before model feature generation and generation |

Selected input spectra and complete own-library exclusion rows are recorded separately. Query signatures bind both. The fresh and retrospective model-study runners recompute the appropriate inference features for each condition. Removing truth from an already-scored C2 table does not establish C3 behavior for generation-dependent models.

## D1 and D2 evidence

D1 reference eligibility includes ionization compatibility. D2 adduct/mass/reference exclusions are applied before top-k selection. The historical `legacy` reference-selection policy remains the default; `source_balanced` is a research option for retaining cross-source references within the same per-identity/adduct budget. It has no general MRR benefit established by the completed studies.

Decision-grade D2 comparisons require verified pool fingerprints, matching complete input manifests and a freshly constructed analog index. Legacy cache files and source-balanced file stamps do not establish cryptographic matrix derivation. A reproducible local run can explicitly request:

```bash
python scripts/run_validation.py --analog --require-fingerprint-provenance \
  --fresh-analog-index --source-manifest runs/source-snapshot.json \
  --library enveda-np-examples --max-queries 250 --max-spectra 5 \
  --query-seed 0 --tag fresh-d2-reference
```

Use the declared study runner and frozen bundle to reproduce a model study; the generic validation CLI has its own baseline scope and does not reproduce fresh neural/ranker fitting.

## Spectral controls

`blank` removes fragment arrays while retaining precursor and acquisition metadata. `intensity_shuffle` permutes intensities among fixed fragment positions using molecule/spectrum/seed-dependent randomness. Both operate on query copies, leaving the source store intact. Candidate generation and downstream features can respond to the changed evidence.

Controls measure dependence and expose shortcuts. They do not certify causality or generalization. Blank fragments reduce learned-model performance in the completed studies but improve some simple C2 baselines; those adverse findings remain part of the reports. `compare_control_runs.py` is the descriptive control comparator; ordinary model comparison must not promote a manipulated-spectrum run as an improved model.

## Pairing, uncertainty and information boundaries

Comparisons require aligned unique query identities, matching original spectra/exclusions and compatible source manifests. Legacy unsigned or unmatched runs cannot support a conclusive comparison. Molecule-level bootstrap preserves pairing; fixed-grid simultaneous intervals cover only their declared grid. The fresh study's 27 secondary contrasts share a separate approximate max-deviation bootstrap family. Neither procedure corrects unknown historical selection or domain shift.

The fresh study separates selected supervised query families and metric identities before fitting. Evaluation chemistry can still occur as training negatives or in task-defined references; the fixed NP-score is an external chemical prior. The 600/150/250 partition is internal and family-balanced, with substantial source concentration. It is not an external laboratory replication.

## Coverage and error slices

Candidate truth coverage bounds attainable MRR; ranking cannot recover a missing candidate. The [loss report](../../reports/candidate_loss_review.md) decomposes ideal MRR into observed score, within-coverage ranking loss and absent-candidate loss. Conditional metrics retain the covered-subset denominator and do not replace the full-population score. Aggregate-only evidence cannot recover conditional bootstrap intervals or causal effects of changing generation.

Source, adduct, ionization mode, mass, spectrum-count and instrument slices are descriptive. The fresh study suppresses cells below 20 molecules. Its joint secondary intervals do not provide multiplicity protection for exploratory subgroup selection.

## Separate purged-reference estimand

The [purged-reference experiment](../../reports/purged_external_review.md) uses disjoint internal and public-source external test cohorts. Both fixed generator budgets are measured even when development recommends one. Its primary contrast is paired internal C3 truth coverage, not MRR; all 11 secondary comparisons share joint max-deviation intervals. Molecules are paired within each cohort and resampled independently across cohorts. A positive primary coverage result does not establish a ranking gain. Original-snapshot novelty feasibility and laboratory authentication are separate questions.
