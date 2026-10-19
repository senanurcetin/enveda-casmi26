# NEXUS CASMI26 research interface — design QA and provenance

**Reviewed 19 October 2026 · Research presentation only · No productization**

The existing `site/index.html` is the project's NEXUS-themed static research
intelligence surface. This design pass preserves the established midnight
background, cyan → blue → lavender/violet accents, restrained cards and
the **official transparent NEXUS SVG from the user's Google Drive identity
folder**. The repo's `site/assets/nexus-official.svg` content was checked
against the approved Drive SVG and matched exactly.

## Delivered layout and usability

- Fully responsive **research evidence readiness** cards distinguish
  verifiable repository work from Kaggle account access, authorized data
  reruns, deliberately unrun CI and the dated public leaderboard result.
- The six display states are loaded from **`site/evidence-status.json`**,
  a dated, human-reviewable static research ledger, not live Kaggle data.
  A missing or incompatible JSON source displays a caution rather than
  inventing a success status.
- Research work is separated from Kaggle's actual Notebook-based evaluation.
  Source links lead to public research protocols and the official competition.
- Keyboard **skip navigation**, strong visible focus outlines,
  active-section navigation, mobile layouts, high-contrast support,
  reduced motion (existing site) and an uncluttered print layout.
- JavaScript renders evidence labels and descriptions using
  `textContent`, not untrusted `innerHTML`, and makes no authenticated
  API calls.

## What the dashboard means

**0.422 public MRR@25** is a **historical snapshot from 2026-10-08**.
The dated card is not a new leaderboard measurement. Real Kaggle/BigQuery
access, new C1/C2/C3 scores, independent validation and full provider
provenance remain separate research tasks. A static site must not imply
those tasks have succeeded merely because safeguards exist in source.

The Kaggle official code competition currently requires a committed
Notebook with **internet disabled**, CPU/GPU runtime <=9 hours, and
`submission.csv`. The competition rules must be accepted in the
participant account; only the competition site can validate a submission.

## Reproduce locally

From repository root:

```sh
python -m http.server 8000
# Visit http://127.0.0.1:8000/site/
```

Run the static accessibility/provenance checks without CI:

```sh
python -m unittest discover -s tests -p 'test_nexus_site_design_stdlib.py' -v
```

**Verified:** Four standalone tests passed using local Windows
SentinelX Python 3.11. The embedded JavaScript also passed
`node --check`. Existing repository CSS/HTML and the approved NEXUS
SVG were checked for correct file references.

**Not independently verified:** a rendered browser screenshot from the
Windows service session could not be generated; production hosting is
not enabled for this repository. The files are available in `site/`
for local static preview. No claim of deployment, live telemetry or
real-data model success is made.

## Secure SentinelX / Kaggle runtime status

The connected Windows host can read the permitted research workspace,
and Python 3.11/pandas/pyarrow/Kaggle tooling is installed. The
SentinelX Windows **service account** cannot read the current user's
private `.kaggle/kaggle.json` credential file. Its Kaggle API
authentication is therefore **blocked**. The host is not given
broader access to the user's home folder. Repeated failed authorization
must never be reported as completed data acquisition.

A human-approved Kaggle login/credential delegation is required before
fetching the licensed files, rebuilding source artifacts and running
independent C1/C2/C3 evaluations. Keep the key in an approved secret
store; **do not put it in Git, screenshots, public logs or chat**.
