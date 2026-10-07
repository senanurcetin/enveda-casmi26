# CASMI26 hand-off (first written 2026-10-06; update 2026-10-07 12:05 UTC at the end)

Goal: climb the public leaderboard of
[enveda-CASMI26-molecule-id-mass-spectra](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra)
(deadline 2026-12-14; top 3 is the aim). Our best own pipeline: public LB 0.309, rank 1475 of 2602.
Leaderboard top at the time of writing: 0.473, 0.466, 0.443.

## Credentials (not in this file)
Give the new session, as environment variables of its environment (never commit them):
- `KAGGLE_API_TOKEN` = the `KGAT_...` access token (rotate it: it was pasted into a chat), or
- `KAGGLE_USERNAME=senanuretin` and `KAGGLE_KEY=<legacy key>`.
Kaggle username: `senanuretin`. Competition slug for the CLI: `enveda-CASMI26-molecule-id-mass-spectra`.

## Kernels (all private, GPU T4, internet off)
| slug | what | state |
|---|---|---|
| `senanuretin/casmi26-base-v4n-pubchem` v1 | verbatim copy of the public notebook "CASMI26 v4n Fusion + PubChem join" (huseyinemreaksoy, author's LB 0.420) | COMPLETE; submission ref 56879149 PENDING since 11:15 UTC |
| `senanuretin/casmi26-fusion-ours` v2 | the same + our pipeline as a third RRF list (`OURS_BETA` 0.4, only molecules with lib_max < 0.7) | commit run COMPLETE (its assert proves our channel produced lists); submission ref 56880555 got "A system error, try resubmitting" after 20 min; resubmitted as ref 56880969 at 12:44 UTC, PENDING |

v1 of `casmi26-fusion-ours` failed immediately: the inserted cell had no `id` (nbformat 4.5). Fixed in v2.

## Pending checks
- BASE score expected ~0.42; ~0.39 would mean ICEBERG/GLACIER silently failed (Python 3.13 wheel issue).
- OURS vs BASE delta decides the next step: positive -> try beta 0.8; negative -> beta 0.2.
- A second OURS `ERROR` would point at our notebook (e.g. run time over the limit on the hidden test, disk, memory), not Kaggle.
- Daily limit: 3 submissions left on 2026-10-06 (resets 00:00 UTC). Never resubmit a duplicate while one is PENDING.
- Scheduled self check-in (routine `trig_01SZoytww6DXtVWBkuoBeJsb`) fires 13:40 UTC in the old session.

## How things are built
- Source of the base notebook: `kaggle kernels pull huseyinemreaksoy/casmi26-v4n-fusion-pubchem-on-public-0-421 -p DIR -m`.
  Third-party code is not committed.
- `notebooks/build_fusion_notebook.py SRC.ipynb OUT_DIR --meta KERNEL-METADATA.json --beta 0.4` writes the patched notebook and
  metadata; then `kaggle kernels push -p OUT_DIR`.
- `kaggle/ours_channel.py` runs our pipeline inside another notebook (fail-open).
- `scripts/run_kernel_submit.sh` waits for a kernel and submits it once (KNAME env var). It dies with its shell: launch it
  with nohup/setsid, or submit by hand with `kaggle competitions submit -c ... -k ... -v ... -f submission.csv -m ...`.
- Score and error text: `KaggleApi().competition_submissions(...)` has `error_description`; kernel status has no failure text.

## Network block
Kernel logs / outputs live on `www.kaggleusercontent.com:443`. The session proxy answers 403 to CONNECT
(`curl -sS "$HTTPS_PROXY/__agentproxy/status"` -> `recentRelayFailures`). That is the only host seen so far, but the
download stops at the first denied host, so after allowing it re-run
`kaggle kernels output senanuretin/casmi26-fusion-ours -p DIR` and look at the host in any new failure: Kaggle may
redirect the log to another storage domain (e.g. a `*.googleapis.com` bucket) that also needs allowing. Network changes
only apply to a NEW session. Do not route around the policy.

## Ideas, in order of expected value
1. Read OURS vs BASE; tune `OURS_BETA`. Our list is an independent second opinion, so a gain of +0.005..+0.02 is plausible, not more.
2. Knobs the public notebooks report as positive but unstacked: `GL_MH_ONLY`, `POP_MU` variants, `FILL_25`, larger `PC_JOIN_N`,
   DreaMS tail channel (dataset `hengck23/hengck23-dreams-enveda-casmi26`). One change per submission; the public LB is the only
   reliable judge (the notebook's own validation mode is optimistic because the public models were trained on all of train.parquet).
3. Top 3 needs ~0.45-0.47: beyond all public recipes. Do not promise it.

## Update 2026-10-07 12:05 UTC

### Public LB results (all on the plain public base, public datasets only)
| submission | public LB |
|---|---|
| BASE (copy of huseyinemreaksoy v4n fusion + PubChem join) | **0.420** |
| BASE + GL_MH_ONLY | 0.419 |
| BASE + POP_MU 0.35 (base uses 0.25) | 0.414 |
| BASE + FILL_25 | pending (ref 56908669) |
| OURS (our ranker as 3rd RRF list), 4 tries | all ERROR "A system error" |

### Why OURS failed, and the fix being tested
Kaggle log of a control commit: `ERRORED_MOUNTING_DATASET ... senanuretin/casmi-models ... dataset loading failed`;
one failure came mid-run, consistent with lazy reads of our datasets failing. Byte-identical fresh copies were
uploaded as `senanuretin/casmi-{code,models,coconut,pubchemlite}-b` (private, they hold train-derived artifacts) and
`senanuretin/casmi26-fusion-ours-d` uses only those (`notebooks/build_fusion_notebook.py --drop ... --add ...`).
`notebooks/build_knob_notebook.py` builds CFG-only variants of the base with public datasets.

### Research findings (web, unverified rules text: a human should read the Rules page)
- Public LB is ~33% of the hidden test (~130 molecules, SE ~0.035); private ~270 molecules. Public-LB differences
  under ~0.02 are noise, and the top-5 public order is about 1 SE wide. Choose the 2 final submissions by a
  validation that mirrors the test, plus one differently-risked hedge.
- External data: standard Kaggle wording (publicly available, free to all) per search snippets. Keep datasets
  derived from train private. Max 5 submissions/day, 2 final picks, team merge deadline 2026-12-07.
- Public design space saturates near 0.42 (Class 3 ~0). Best next ideas, gains are guesses: listwise
  same-formula-decoy head on an existing spectrum embedding (~1 day), candidate-InfoNCE head on frozen DreaMS +
  ChemBERTa (MSAlign recipe, ~3 days), NP prior from `dmitriigluzdov/casmi26-natural-product-knowledge-table`
  (LOTUS/NPAtlas, CC0/CC BY). Beware shortcut features (SMILES format, library membership).
- Useful public datasets: `hengck23/hengck23-dreams-enveda-casmi26` (DreaMS weights + Enveda-180 embeddings),
  `prvsiyan/casmi26-fp-models-v4`, `dmitriigluzdov/casmi26-natural-product-panel-benchmark`,
  `dmitriigluzdov/casmi26-fold-safe-fpnet` (needed to validate the public pipeline without leakage).
