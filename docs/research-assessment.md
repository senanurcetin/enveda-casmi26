# Research assessment

Reviewed 27 October 2026. CASMI26 is an advanced, noncommercial research software and data science portfolio for molecular identification from MS/MS spectra. The assessment concerns reproducibility, experimental discipline and evidence boundaries; it does not claim state-of-the-art accuracy or deployment readiness.

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

The [external-source feasibility study](studies/external-source-validation.md) has now tested two independent downloads against every training row: MassBank leaves zero novel molecules among 491 eligible identities, and two GNPS2 standard libraries leave seven among 85. Neither meets the predeclared minimum of 30. This is measured negative feasibility evidence; external model accuracy remains unmeasured. C3 expansion preserves historical candidates in implementation tests, but no real coverage-improvement claim is established yet.

External acquired spectra or independent laboratory replication would address evidence not supplied by internal splitting. C3 generation recall and broader cross-instrument samples are justified research priorities. Any follow-on experiment must declare its population and selection policy before results, use a separate evaluation set and preserve all favorable/adverse outcomes. The completed test must not become a tuning set.

There is no planned billing, hosted inference, commercial service or production deployment. [Verification scope](verification.md) · [Evidence registry](../reports/research_evidence.json).
