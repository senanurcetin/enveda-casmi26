# Fresh training and internal chemical-family holdout

Status: protocol registered before execution in commit `4f8ed74`. Runtime completion and results must be verified separately; implementation is not experimental evidence.

The [fixed protocol](../protocols/advanced-internal-v1.json) declares 600 training, 150 development and 250 final-test molecules. Generic Murcko ring topology, generic acyclic connectivity and metric identity prevent supervised query identities/families crossing partitions. Historical V-in families are excluded. One representative per family defines a family-balanced internal research population, not prevalence-weighted real-world accuracy.

Two compact fingerprint transformers are initialized randomly and trained on training-partition identities only. Fixed six-epoch training does not use development/test outcomes for checkpoint selection. Ranker-training queries receive the opposite fold's D3 model; development/test use the two-model ensemble. Eight prespecified three-seed rankers use development data for stopping. Every learned artifact, feature list and partition digest is sealed before any test features are computed.

The primary contrast is unperturbed test C2 MRR@25 of the fresh full model minus D1+D2. Secondary feature ablations and spectral controls receive joint, molecule-paired max-deviation bootstrap intervals over the prespecified comparison family. Source/adduct/mode/mass/instrument slices are descriptive and suppress cells below 20 molecules. Candidate-coverage ceilings separate retrieval misses from ranking misses.

Feature ablations remove a channel and its derived rank/gap/z-score context. They share the same generated candidate universe: `no_gen_features` removes generator features, not generated candidates. `membership_only` is a simulation-shortcut diagnostic and is never selected as a deployable model.

Partition disjointness concerns supervised fitting. The retrieval corpus remains available according to the task definition: C1 may retain same-structure cross-library references, while C2/C3 exclude metric-equivalent truth references. This is not a reference-free discovery benchmark or a newly acquired external/laboratory dataset. Generic topology groups are conservative computational partitions and do not guarantee distant chemical similarity.

The distinct [owned notebook](https://www.kaggle.com/code/senanuretin/casmi26-fresh-training-family-holdout) rebuilds source artifacts in scratch storage, attaches no pretrained-model dataset, and uses existing free Kaggle GPU execution with internet disabled. Only aggregates, hashes, training lineage and exact package versions may be exported. No raw molecular records, per-query ranks, weights, commercial infrastructure or competition submissions are part of this study.

Historical [retrospective v16 evidence](../reports/revalidation_review.md) is preserved against the actual source embedded in its executed notebook. Current implementation changes do not rewrite that experiment or retroactively attribute new code to its results.
