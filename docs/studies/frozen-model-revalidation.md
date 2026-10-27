# Frozen-model retrospective revalidation

Completed on Kaggle on 24 October 2026, with source/status/artifact readback authenticated on 25 October. The study retains all 250 eligible designated V-in molecules and compares D1, D1+D2 and the frozen five-seed v16 ranker. It fits no new weights and makes no competition submission.

[Full results and paired controls](../../reports/revalidation_review.md) · [Aggregate artifact](../../reports/revalidation_study.json) · [Publication verification](../../reports/revalidation_publication.json)

## Design

| Dimension | Executed setting |
|---|---|
| Population | All 250 eligible `enveda-np-examples` molecules |
| Query sampling | At most five spectra per molecule, seed 0 |
| Models | D1; D1+D2; frozen v16 with attached original owned fingerprint/ranker weights |
| Conditions | C1/C2/C3, each unperturbed, blank and intensity-shuffled |
| Absolute intervals | 2,000 molecule-bootstrap draws |
| Paired comparisons | 4,000 draws and the declared novelty-mixture grid |

The shared inference path recomputes features and generated candidates after the relevant exclusions. Complete own-library exclusion includes unsampled acquisitions and metric-equivalent variants. Blank/shuffle controls preserve their stated precursor/acquisition information.

## Source reconstruction

The scratch-only build quarantines 1,710 mass-conflicting identities (3,891 provider rows) and 13 unparseable selected candidates, verifies 1,071,863 candidate fingerprints and rebuilds a cache-free D2 index containing 1,062,995 references. Labels are retained independently of candidate quarantine; one query's standard truth key is absent from the pool and remains in every denominator. Complete own-library exclusion adds 287 reference rows beyond sampled inputs for 74 molecules.

## Reproduction and provenance

```bash
python scripts/build_revalidation_notebook.py
python scripts/check_revalidation_publication.py
python scripts/build_revalidation_report.py
```

The completed bundle under `kaggle/audits/corrected-revalidation/` is frozen against its executed source. The checker authenticates public completion, matching source/input attachments and the bounded aggregate artifact. Reusing a URL or regenerating a report does not execute a new experiment. Fresh execution requires a separately declared source bundle and available licensed inputs.

## Limits

V-in informed historical model development. The original weights' complete training lineage is not independently attested. Corrected exclusions cannot undo model selection, establish an isolated effect for each fix or remove instrument/domain shift. The source snapshot retains its connectivity-only generated cache; the representation-sensitive correction applies to the later fresh-training study. Historical public scores remain unchanged.
