# Corrected models and negative controls

This study evaluates frozen original model weights retrospectively on every eligible molecule in the designated `enveda-np-examples` population. It does not create a new untouched holdout, train a new model, evaluate public test examples or submit predictions.

**Completed execution:** 24 October 2026; authenticated source/status/artifact readback: 25 October 2026. All **250** molecules were retained. [Results, paired intervals and adverse control findings](../reports/revalidation_review.md) · [Aggregate evidence](../reports/revalidation_study.json) · [Publication verification](../reports/revalidation_publication.json).

## Source rebuild and quality policy

The owned [Kaggle notebook](https://www.kaggle.com/code/senanuretin/casmi26-corrected-models-and-negative-controls) rebuilds the complete spectrum store, candidate pool, Morgan fingerprints and D2 reference index from attached inputs. The official RDKit 2026.3.3 wheel is attached; computation has internet disabled. Models and licensed derived files remain in scratch storage. Exported evidence contains aggregates and hashes only.

The first execution stopped at the candidate-source quality gate: **1,710 shared InChIKey14 groups had conflicting neutral masses**, beyond the existing maximum of 0.01 Da or 20 ppm. The revised study explicitly opts into **whole-identity quarantine**. Every provider row for a conflicting identity is excluded. Training-source precedence cannot conceal the disagreement. Ordinary pool builds continue to reject these conflicts unless the explicit quarantine option is enabled.

Quarantine does not delete evaluation molecules. Ground-truth structures come from the independently retained training structure table. Metric-equivalent references are excluded using both that table and the retained candidate pool. A missing candidate truth remains in the denominator and can only score if the inference path actually recovers it. No labels are inserted into candidate generation.

The following execution reached fingerprint generation and stopped with a `ValueError`. The revised build also screens selected candidate structures before publication: unparseable SMILES, no-heavy-atom structures and zero/error Morgan vectors are explicitly quarantined and counted. The fingerprint builder remains strict and performs its own verification. This screening establishes computational usability, not provider-identity accuracy or chemical plausibility.

The completed build quarantined **13 unparseable selected candidates**, leaving **1,071,863** verified structures. Mass quarantine removed **3,891 provider rows** for the 1,710 conflicting identities. One query's standard truth key was absent from the retained pool; that molecule remained in every scoring denominator.

The fresh D2 index contained **1,062,995 reference spectra**. Complete source-library exclusions added **287 rows beyond the selected query inputs for 74 molecules**. Those counts describe exclusion scope, not an isolated estimate of how much the earlier bug changed accuracy.

## Paired experimental design

| Dimension | Prespecified settings |
|---|---|
| Population | All eligible designated V-in molecules; no maximum-query subsample |
| Query sampling | At most five spectra per molecule, seed 0 |
| Models | D1; D1+D2; frozen five-seed v16 ranker with both fingerprint models and required fragment/generation channels |
| Novelty simulation | C1 complete same-library exclusion; C2 all metric-equivalent reference exclusion; C3 C2 plus candidate exclusion before feature generation |
| Spectral conditions | Unperturbed; blank fragments; deterministic intensity shuffle |
| Absolute uncertainty | 2,000 molecule bootstrap draws |
| Paired diagnostics | 4,000 molecule bootstrap draws; fixed novelty-mixture grid |

The study calls the shared `Predictor.rank` implementation. Reference-count priors are cleared for withheld structures in C2/C3. C3 features and generated candidates are recomputed after exclusion; filtering C2 predictions afterward is insufficient.

Blank controls retain precursor, adduct, collision energy and mass-constrained candidate information. Intensity shuffling retains fragment positions. These diagnostics describe dependence on fragment evidence; they do not prove causality, training independence or external generalization.

## Verification

```bash
python scripts/build_revalidation_notebook.py
kaggle kernels push -p kaggle/audits/corrected-revalidation
python scripts/check_revalidation_publication.py
```

The readback command requires a completed public run with matching notebook/source hashes, the declared input attachments, disabled internet and no submitted predictions. It rejects incomplete summaries, inconsistent denominators and row-level query/rank exports. A completed readback writes `reports/revalidation_study.json` and `reports/revalidation_publication.json`.

## Interpretation boundary

The designated V-in population informed historical experiments. Frozen weights are not retrained, and their complete original training lineage is not independently attested. Rebuilt inputs and corrected knockouts cannot reverse historical model selection. Confidence intervals do not correct that bias or instrument/domain shift. The study makes no ACCEPT/REJECT promotion decision and replaces no historical public or private competition score.
