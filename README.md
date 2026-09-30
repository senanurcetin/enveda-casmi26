# Enveda CASMI 2026 — molecule identification from MS/MS spectra

Kaggle competition: [enveda-CASMI26-molecule-id-mass-spectra](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra).
Given the tandem mass spectra of an unknown natural product, rank up to 25 candidate structures;
scored by MRR@25 on the tautomer-canonical InChIKey (first block).

**Best public leaderboard: 0.309** (`ranker-v16-gen2`, notebook `casmi26-probe-gen3`).

## Approach in one paragraph

Retrieve, then rank. Candidates come from a structure pool (training structures + COCONUT +
PubChemLite, ~1.08 M) within ±10 ppm of the query's neutral mass, plus *generated* candidates:
biotransformations (methylation, glycosylation, homologs, acyl/sugar losses, …) applied to the
query's spectral analogs, for molecules that are in no database. Each candidate gets evidence from
four channels, and a LightGBM LambdaRank model trained on knockout-simulated queries ranks them:

| channel | evidence |
|---|---|
| D1 | direct library match (entropy similarity to the candidate's own spectra) |
| D2 | mass-shifted analog propagation: similarity to other structures' spectra × Tanimoto |
| D3 | spectrum → fingerprint network (set transformer), Bernoulli likelihood of the candidate fingerprint |
| D4 | fragment explanation: share of the spectrum explained by 1–2 bond cleavages of the candidate |

## What made the difference (leaderboard)

| step | public LB |
|---|---|
| direct match + analog rule | 0.190 |
| learned ranker, membership priors removed (they were simulation shortcuts) | 0.258 |
| + PubChemLite in the pool, ranker retrained on it | 0.273 |
| + D4 fragment explanation + Class-3 generated candidates | 0.298 |
| + extended transformations, 2× training queries | **0.309** |

Just as informative are the changes that lost on the leaderboard: bigger pools (ChEBI, LIPID MAPS,
PubChem literature), up-weighting Class 3, NP-likeness, tuned LightGBM, a larger D3, a learned bond
model and a learned spectral similarity. Each one is logged with its numbers in
[EXPERIMENTS.md](EXPERIMENTS.md).

## Validation

Two held-out query sets, both at structure (skeleton) level and never used for training:
**V-in** (the 250 instrument-matched `enveda-np-examples`) and **V-real** (600 harder, held-out
library structures). V-real tracks the leaderboard's direction for pool and channel changes;
effects below ~0.01 are beyond what either set resolves. Details and every decision in
[EXPERIMENTS.md](EXPERIMENTS.md); the design in [ARCHITECTURE.md](ARCHITECTURE.md).

## Reproduce

Competition data is not redistributed. Everything train-derived is rebuilt inside the Kaggle
notebook by [kaggle/submission_driver.py](kaggle/submission_driver.py) from `train.parquet` and the
public sources (COCONUT CC0, PubChemLite CC BY 4.0):

```
CASMI_RANKER=ranker-v16-gen2  CASMI_USE_PUBCHEMLITE=1  python kaggle/submission_driver.py
```

Locally, the same stages run as `scripts/build_*.py` → `scripts/make_submission.py`; ranker training
data comes from `scripts/build_rank_data.py --generate`, `scripts/merge_train_gen.py` and
`scripts/train_ranker.py`.

## Public notebooks from this project

- [CASMI26: Your Ranker Learns How You Built the Pool](https://www.kaggle.com/code/senanuretin/casmi26-your-ranker-learns-how-you-built-the-pool)
- [CASMI26: A Bigger Candidate Pool Lowered My Score](https://www.kaggle.com/code/senanuretin/casmi26-a-bigger-candidate-pool-lowered-my-score)
