# Corrected models and negative controls

This study evaluates frozen original model weights retrospectively on every eligible molecule in the designated `enveda-np-examples` population. It does not create a new untouched holdout, train a new model, evaluate public test examples or submit predictions.

## Source rebuild and quality policy

The owned [Kaggle notebook](https://www.kaggle.com/code/senanuretin/casmi26-corrected-models-and-negative-controls) rebuilds the complete spectrum store, candidate pool, Morgan fingerprints and D2 reference index from attached inputs. The official RDKit 2026.3.3 wheel is attached; computation has internet disabled. Models and licensed derived files remain in scratch storage. Exported evidence contains aggregates and hashes only.

The first execution stopped at the candidate-source quality gate: **1,710 shared InChIKey14 groups had conflicting neutral masses**, beyond the existing maximum of 0.01 Da or 20 ppm. The revised study explicitly opts into **whole-identity quarantine**. Every provider row for a conflicting identity is excluded. Training-source precedence cannot conceal the disagreement. Ordinary pool builds continue to reject these conflicts unless the explicit quarantine option is enabled.

Quarantine does not delete evaluation molecules. Ground-truth structures come from the independently retained training structure table. Metric-equivalent references are excluded using both that table and the retained candidate pool. A missing candidate truth remains in the denominator and can only score if the inference path actually recovers it. No labels are inserted into candidate generation.

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
