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
