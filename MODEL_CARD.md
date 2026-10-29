# Model card — CASMI26 retrieval and ranking

**Owner:** Senanur Cetin / NEXUS. **Scope:** noncommercial molecular identification research. **Reviewed:** 29 October 2026.

## Intended use and components

Rank up to 25 candidate structures per molecule from MS/MS spectra. The pipeline combines mass-compatible retrieval, D1 reference matching, D2 spectral analogs, D3 learned fingerprints, D4 fragment explanations, candidate generation and LightGBM LambdaRank. Rank scores are not calibrated identification probabilities. [Implementation](ARCHITECTURE.md) · [Evaluation semantics](docs/guides/evaluation.md).

MRR@25 uses the first tautomer-canonical InChIKey14 match within the original 25 slots. Missing predictions score zero; invalid, duplicate and blank slots are not compressed. This metric does not establish stereochemical identification or laboratory identity.

## Historical competition attribution

| System | Author-reported public MRR@25 | Attribution |
|---|---:|---|
| Own `ranker-v16-gen2` pipeline | 0.309 | Original retrieval/ranking work |
| External v4n base | 0.420 | Third-party solution |
| `casmi26-fusion-ours-d` | 0.422 | External base plus own rank-list fusion |

These are [8 October 2026 author logs](EXPERIMENTS.md), not a current leaderboard readback or final private score. The +0.002 fusion difference has no established statistical significance or transfer guarantee. The dashboard serves no trained model.

## Completed model evaluations

| Experiment | Fitting and evaluation boundary | Principal result | Full evidence |
|---|---|---|---|
| Frozen-model retrospective study | Rebuilt sources; original v16 weights; all 250 eligible V-in queries | v16 C1/C2/C3 MRR 0.8324/0.6241/0.2374 | [Report and controls](reports/revalidation_review.md) |
| Fresh internal family benchmark | 600/150/250 train/development/test molecules; two fresh D3 networks and eight three-seed rankers; 26 learned files sealed before test features | Full C1/C2/C3 MRR 0.8665/0.8388/0.0373 | [Report, all 27 secondary contrasts and slices](reports/advanced_review.md) |
| Separate purged-reference experiment | Selected external families/identities removed from the full training/reference corpus; 300/100/150 internal partitions and 150 external test molecules; eight fresh learned files sealed before test | Expanded external C2/C3 MRR 0.4454/0.0466; internal C3 coverage 13/150 to 17/150 | [Report and all 11 secondary contrasts](reports/purged_external_review.md) |

The fresh internal study's primary C2 full-minus-D1+D2 difference is +0.6322, paired 95% CI [+0.5816, +0.6829]. Mass-only C2 exceeds D1+D2, which limits interpretation of that baseline comparison. C2 D1/D2/D4/generator-feature contributions are supported within the fixed secondary family; incremental D3 and NP-score effects remain unresolved.

The separate purged-reference study's primary internal C3 coverage change is +0.0267, paired 95% CI [+0.0067, +0.0533], with four gains and zero losses. Expanded-minus-legacy ranking effects remain unresolved in the joint 11-contrast intervals. External C3 truth coverage is only 23/150. Coverage is not achieved ranking accuracy. Different populations, source corpora and fitting prevent a causal before/after comparison between studies.

## Evidence boundaries and failure modes

- Original-snapshot novelty feasibility failed: MassBank 0/491, two-library GNPS 7/85 and expanded GNPS 7/936 novel molecules. The two GNPS stages overlap and are not independent replications. The separate family-purged result measures generalization relative to a rebuilt corpus; it does not reverse these findings. [Source methods](docs/studies/external-source-validation.md).
- Selected supervised query families/metric identities are separated, but task-defined candidate knowledge remains shared. Evaluation structures may occur as negatives in other training queries. The fixed NP-score prior is not freshly fitted and its upstream overlap is unverified.
- All 1,213 public-example spectra match training content. They cannot establish independent generalization. [Content audit](docs/studies/public-example-overlap.md).
- Public-source identity and acquisition annotations lack independent laboratory authentication. Restricted chemistry/instruments and unavailable numeric eV encoded as zero limit external transfer; NCE/percent values are not converted to eV.
- The original internal study lacks the C3 truth for 224/250 molecules. Its fixed-list oracle MRR ceiling is 0.1040, versus achieved MRR 0.0373. Conditional performance on covered queries cannot replace the full population denominator. [Loss accounting](reports/candidate_loss_review.md).
- Frozen-model controls diagnose the original weights; their full training lineage is not independently attested. Current source corrections do not silently rewrite historical executed snapshots.
- Isomers, incorrect adducts, contamination, missing references and generated-structure errors can produce confident-looking incorrect rankings. Negative controls are diagnostics; some simple baselines improve under blank spectra.

These results do not establish state-of-the-art performance, independent laboratory replication, production readiness or suitability for consequential chemical decisions. Source code, datasets, third-party models and the NEXUS identity have separate rights and attribution. [Data card](DATA_CARD.md) · [Dated verification](docs/verification.md).
