# External-source validation and C3 expansion

This follow-on study separates source eligibility, training-data independence and model performance. A source-preparation manifest establishes the downloaded inputs and filters. A completed overlap audit establishes whether a usable external population exists. Neither establishes model accuracy or a candidate-generation improvement.

## Registered design

The [original protocol](../../protocols/external-c3-v1.json), committed before new outcomes at `6164c5e`, fixes two generator budgets: the historical 15 analog seeds / 300 proposals and an expanded 30 seeds / 900 proposals with both transform-pair orders. Expansion preserves historical candidates and their intrinsic features. Candidate-list context is recomputed separately, and `audit_*` membership markers cannot enter fitted features. Ordinary inference keeps its historical default.

The planned experiment uses new 300/100/150 internal training/development/test families, excluding every selected family from the completed 600/150/250 study and historical V-in families. Two freshly initialized fingerprint networks and two three-seed rankers are specified. Development C3 coverage determines the retained policy before all weights, feature lists and populations are sealed. Both policies and all adverse test outcomes must be reported. [Execution details](../../protocols/external-c3-v1-execution.json).

External spectra must pass full-training standard identity, official tautomer-metric identity and exact stored-content comparisons. At least 30 novel molecules must remain after the relevant exclusions. No model prediction or candidate-coverage result may choose the source population. At most five accession-sorted spectra are retained per selected molecule. Unknown numeric collision energy is encoded as zero; this is a disclosed distribution-shift limitation, not an eV/NCE conversion.

## Executed source feasibility

| Source | Eligible molecules / spectra | Excluded molecules | Novel molecules | Conclusion |
|---|---:|---:|---:|---|
| MassBank 2026.03, explicit CC0 records | 491 / 4,567 | 491 | 0 | Independent population not established |
| Two GNPS2 direct-contribution standard libraries | 85 / 95 | 78 | 7 | Below the declared minimum of 30 |
| Expanded GNPS community plus the two libraries | 936 / 3,051 | 929 | 7 | Below the declared minimum of 30 |

All three source audits scan all **2,539,608 training rows** and report zero unparseable training labels. None found exact stored spectral-content matches. The expanded GNPS source includes the earlier two libraries; these are overlapping source stages, not three independent external replications. Lack of byte-equivalent spectra does not override molecular identity overlap. The MassBank source has 491 standard-identity overlaps; the two-library GNPS source has 78 standard and metric overlaps. These are negative feasibility findings, not failed model predictions.

Evidence: [MassBank aggregate](../../reports/external_audit.json), [authenticated publication](../../reports/external_audit_publication.json), [GNPS aggregate](../../reports/gnps_external_audit.json), [authenticated publication](../../reports/gnps_external_audit_publication.json). Source manifests retain pinned archive/record hashes and license attribution. Only directly contributed GNPS standard spectra use the GNPS CC0 policy; imported, propagated, MassQL-derived and suspect libraries are excluded.

## Source amendments and integrity

The [GNPS source addition](../../protocols/external-source-v2.json) was declared before model-performance inspection. Its first execution used the older connectivity prefilter and is superseded for independence decisions. The authenticated second execution uses the [declared isotope correction](../../protocols/external-source-v2-audit-correction.json): isotopic tautomer equivalents must survive prefiltering because metric InChIKey14 omits isotope layers. A synthetic regression verifies the issue. The current connectivity helper clears isotope labels before removing explicit hydrogen atoms. Frozen historical experiment sources remain preserved.

The [broader GNPS community source](../../protocols/external-source-v3.json) is an outcome-blind extension under the same novelty and minimum-population requirements. Its authenticated audit also leaves only the same seven novel molecules. The original-snapshot novelty design has not established a sufficient external population. [Expanded aggregate](../../reports/gnps_expanded_external_audit.json) and [readback](../../reports/gnps_expanded_external_audit_publication.json).

## Separate family-purged experiment

The [28 October protocol](../../protocols/purged-external-v1.json), registered at `1b7a6cc` before model outcomes, defines a new experiment. Choose at most 150 external molecules before predictions, one per generic-family component, and remove every matching family, standard identity and official metric identity from the full training/reference corpus. Block all rows sharing a removed declared training identity. Rebuild structures, stores, candidate fingerprints and the D2 index; independently verify zero remaining source identity/content overlap and zero retained selected families before fitting.

Historical internal query partitions are reconstructed on the original data first. Their exclusion labels remain protected during selection of fresh 300/100/150 internal partitions on the purged corpus. Two fresh fingerprint networks and two three-seed rankers retain the fixed budgets, development-only recommendation and 11 secondary contrasts. Both budgets and adverse outcomes remain reportable. [Execution details](../../protocols/purged-external-v1-execution.json).

This is external-family generalization relative to a newly purged training snapshot. It does not convert the failed original-snapshot novelty result into a success or independently authenticate source acquisition/annotations. Shared independent candidate knowledge and fixed NP-score priors remain explicit limitations. Attempt 1 failed during source preparation before fitting because generic-family grouping was undefined for some training labels. The [outcome-blind source amendment](../../protocols/purged-external-v1-source-amendment.json) quarantines every affected declared training identity and reports counts; external selection, budgets and all estimands remain unchanged. The [attempt history](../../reports/purged_external_execution_history.json) preserves source hashes and the sanitized failure. Repaired execution and model outcomes are pending.

## Interpretation boundary

C3 expansion is implemented and has chemistry, candidate-preservation and feature-context checks. Its real coverage gain and external model performance remain unmeasured until the declared experiment completes and its aggregate/source readback is authenticated. No result is promoted from synthetic tests or source preparation.

Even a passing external population would establish computational non-overlap against the pinned training snapshot, not independent laboratory authentication. Source annotations, restricted chemistry/instruments, unknown collision energies and shared candidate priors remain limitations. This is a noncommercial research project with no deployment, payment integration, paid compute or competition submission.
