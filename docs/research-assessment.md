# Research assessment

Reviewed 29 October 2026. CASMI26 is an advanced, noncommercial research software and data science portfolio for molecular identification from MS/MS spectra. The assessment concerns reproducibility, experimental discipline and evidence boundaries; it does not claim state-of-the-art accuracy or deployment readiness.

| Area | Current evidence | Remaining limit |
|---|---|---|
| Modeling | Mass-constrained retrieval, D1/D2/D3/D4 channels, generated candidates and LambdaRank | Rank scores are not calibrated identification probabilities |
| Fresh fitting | Two new fingerprint networks and eight three-seed rankers; 26 learned-file hashes sealed before test features | Small fixed architecture/sample; fixed external NP-score prior |
| Internal partitioning | 600 training / 150 development / 250 test queries; selected family/metric identities separated | Shared task-defined candidate/reference knowledge; source concentration |
| Statistical evaluation | Predeclared primary contrast, 27 joint secondary contrasts, paired controls and descriptive slices | No correction for external domain shift or historical hypothesis selection |
| Data engineering | Validated immutable stores, candidate reconciliation, verified fingerprints, explicit manifests and fresh D2 indexes | Computational consistency does not authenticate provider annotations |
| Reproducibility | Hash-locked Python 3.12 profile, exact Kaggle inventories, frozen executed source and CI | CUDA/platform differences may prevent bitwise reproducibility |
| Communication | Source-linked reports, model/data cards and NEXUS research dashboard | Historical competition scores retain author-log attribution |

## Principal findings

The [fresh internal study](../reports/advanced_review.md) scores 0.8665/0.8388/0.0373 MRR@25 for C1/C2/C3. Its predeclared C2 full-minus-D1+D2 difference is +0.6322 with paired 95% interval [+0.5816, +0.6829]. Mass-only C2 exceeds D1+D2, so the declared baseline's weakness matters when interpreting the large difference. D1/D2/D4/generator-feature C2 contributions are supported within the fixed secondary family; D3 and NP-score incremental effects remain unresolved.

The [C3 diagnosis](../reports/candidate_loss_review.md) finds 224/250 truths absent from candidates. Only 26 are covered; 24 reach top 25 and 6 rank first. Perfect ordering of the same candidate lists could reach only 0.1040 MRR. This is descriptive loss accounting, not a new model result.

The [public-example audit](studies/public-example-overlap.md) confirms exact training overlap for all 1,213 public spectra. These examples cannot establish independent generalization. The [frozen-model study](studies/frozen-model-revalidation.md) supplies separate retrospective diagnostics; comparing its scores directly with the fresh study does not isolate improvement because population, fitting and architecture differ.

## Follow-on research

The [external-source feasibility study](studies/external-source-validation.md) scans every original training row. MassBank leaves zero novel molecules among 491 eligible identities; the expanded GNPS source leaves the same seven among 936. Neither reaches the declared minimum of 30. These are measured negative feasibility findings. A separately preregistered experiment now removes selected external families/identities from the full training/reference corpus and fits fresh weights. Its completed outcomes are reported below; original-snapshot novelty is not reclassified as successful.

External acquired spectra or independent laboratory replication would address evidence not supplied by internal splitting. C3 generation recall and broader cross-instrument samples are justified research priorities. Any follow-on experiment must declare its population and selection policy before results, use a separate evaluation set and preserve all favorable/adverse outcomes. The completed test must not become a tuning set.

There is no planned billing, hosted inference, commercial service or production deployment. [Verification scope](verification.md) · [Evidence registry](../reports/research_evidence.json).

## Completed family-purged study

The authenticated family-purged study removed 261,826 of 2,539,608 training spectra and fitted eight fresh learned artifacts. On 150 internal test molecules, expanded C3 truth coverage rose from 13 to 17 (paired change +0.0267, 95% CI [+0.0067, +0.0533]). On 150 external molecules, expanded C2/C3 MRR@25 was 0.4454/0.0466; expanded-minus-legacy ranking effects remain unresolved in the joint 11-contrast intervals. External C3 truth coverage remains only 23/150. Original-snapshot novelty gates remain unsuccessful; this is generalization relative to a rebuilt purged corpus, not independent laboratory replication. [Methods, all outcomes and limitations](../reports/purged_external_review.md).

The explicit CASMI_DATA snapshot takes precedence over automatic Kaggle discovery. All reference artifacts are rebuilt from filtered rows before fitting. Undefined generic-family labels quarantine the entire affected declared identity under the outcome-blind amendment; the failed first attempt remains preserved.
