# MS/MS negative controls — preregistered scientific checks

**Purpose:** Determine whether C1/C2/C3 validation scores rely on genuine
fragmentary spectral evidence or on non-spectral priors (precursor mass,
candidate-pool membership and source-specific library effects).

**Scope:** Non-commercial, research-only CASMI26 validation. The controls
are **not** alternative prediction methods, leaderboard experiments or a
basis for model-improvement ACCEPT/REJECT claims.

## Treatments

| Mode | What changes | What stays fixed | Diagnostic interpretation |
| --- | --- | --- | --- |
| `none` | Nothing | All selected query spectra and knockout protocol | Reference only |
| `blank` | Every query's fragment m/z and intensity array becomes empty | Query molecule, precursor mass, adduct, ionization mode, collision energy, eligible candidates, C1/C2/C3 knockout scope | Residual ranking arises without fragment evidence; does **not** estimate a true competition score |
| `intensity_shuffle` | Intensities are permuted among the same fragment m/z values using a per-molecule deterministic seed | Original peak positions, precursor, adduct, collision energy, molecule population and knockout | Tests sensitivity to the intensity/m/z relationship; not all spectral evidence is removed |

Both transformations preserve source spectra in the original read-only store;
they operate only on **query copies**. Permutations are keyed to molecule ID,
spectrum index and user seed, so query-order changes do not silently change
a molecule's perturbation. Intensity arrays with identical values are
necessarily unchanged by a permutation.

## Running the diagnostics on approved local inputs

Use the **same** frozen source snapshot, source library, query seed and query
subset. A control run records original source query signatures and its own
transformed query signatures.

```bash
python scripts/run_validation.py --analog \
  --library enveda-np-examples --max-queries 250 \
  --max-spectra 5 --query-seed 0 --spectral-control none \
  --tag diagnostic-reference

python scripts/run_validation.py --analog \
  --library enveda-np-examples --max-queries 250 \
  --max-spectra 5 --query-seed 0 --spectral-control blank \
  --control-seed 0 --tag diagnostic-blank

python scripts/run_validation.py --analog \
  --library enveda-np-examples --max-queries 250 \
  --max-spectra 5 --query-seed 0 --spectral-control intensity_shuffle \
  --control-seed 0 --tag diagnostic-shuffle
```

Each run writes a separate, collision-resistant file under `runs/`.
Report fields include the negative-control mode, seed and original per-query
source signatures (for manipulated runs). Compare the **descriptive** C1/C2/C3
MRR and pool recall values to inspect dependence on real peak evidence.

**Important:** `scripts/compare_runs.py` rejects negative-control reports
by design, even if per-query IDs match. Its ACCEPT/REJECT decisions concern
models evaluated on *identical actual input spectra*. Intentionally
perturbing the query spectra violates that condition and must not be
misrepresented as ordinary model-improvement evidence. A dedicated control
analysis should verify original-source signature alignment explicitly before
quantifying uncertainty in treatment effects.

## Paired diagnostic report (source-verified, no model verdict)

The original control implementation recorded source spectra fingerprints but
did not provide an approved *diagnostic* comparator. The separate
`scripts/compare_control_runs.py` now validates:

- Exact same query IDs, irrespective of report row order, with no duplicates.
- Matching `casmi-query-spectra-v2` provenance. In the controlled report,
  `source_query_signatures` must correspond **exactly** to the original
  unperturbed query fingerprints, including excluded reference row scope.
- Identical model configuration and all non-control query-sampling parameters.
- Valid C1/C2/C3 ranks with coherent MRR summary, genuine changed treatment
  signatures, and an explicit valid perturbation mode and seed.

Only when all guards pass are molecule-aligned, descriptive
`negative_control_minus_unperturbed_reference` effects calculated. Each class
has a pointwise paired-bootstrap interval, with an approximately simultaneous
band across the twelve predefined novelty-class mixtures. The output is always
`DIAGNOSTIC_ONLY` and contains **no** ACCEPT or REJECT field. Deliberately
perturbed input spectra are **not** evidence for a model-improvement decision.

```bash
python scripts/compare_control_runs.py \
  runs/ACTUAL_REFERENCE_RUN.json runs/ACTUAL_BLANK_RUN.json \
  --bootstrap 4000 --seed 0 \
  --output runs/blank_diagnostic.json

python scripts/compare_control_runs.py \
  runs/ACTUAL_REFERENCE_RUN.json runs/ACTUAL_SHUFFLE_RUN.json \
  --bootstrap 4000 --seed 0 \
  --output runs/shuffle_diagnostic.json
```

These commands name placeholders: use actual generated validation paths.
Report creation uses exclusive file creation to avoid overwriting another
diagnostic. CI tests use synthetic reports only; they do not calculate
CASMI26 experimental control scores.

## Leakage and confounding cautions

- A strong `blank` score may indicate candidate-mass/pool priors, not
  necessarily forbidden label leakage. A suspicious effect motivates a
  separate source-membership and train/holdout contamination audit.
- Intensity shuffling retains fragment m/z values, adduct and collision
  energy; therefore a positive score after shuffling does **not** mean
  the model is insensitive to spectra.
- Control runs cannot remove contamination from a downstream ranker
  trained with leaked labels. C1/C2/C3 reference exclusion, out-of-fold
  training and data snapshot integrity must be audited independently.
- No GPU/Kaggle data execution or leaderboard measurements were performed
  when adding these controls. **No claims of model improvement are made.**

## Test gates

```bash
python -m unittest discover -s tests -p '*_integration.py' -v
python -m unittest discover -s tests -p '*_stdlib.py' -v
```

The added tests check deterministic perturbations, unchanged source inputs
and metadata, invalid controls, and inability of a manipulated run to
obtain a model ACCEPT/REJECT verdict from the standard comparator.
