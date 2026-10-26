# Reproduce and review

## View without scientific dependencies

From the repository root, run `python -m http.server 8000`, then open `http://localhost:8000/site/`. Serve the repository root so the dashboard can load `reports/`. There are no external fonts, tracking scripts, credentials or paid APIs.

## Offline checks

```bash
python -m pytest tests/*_stdlib.py -q
python scripts/research_evidence.py
```

Pytest loads `src` through `pytest.ini`; on shells without wildcard expansion use `python -m pytest -q` for the full suite. For a dependency-free run set `PYTHONPATH=src` in your shell, then run `python -m unittest discover -s tests -p '*_stdlib.py' -v`. Windows skips the symlink-containment case only when the OS denies link creation; Linux CI executes it.

The evidence command verifies local source links, allowed evidence states, score ranges, unique IDs and source hashes; it writes `reports/research_review.md`. It cannot authenticate historical leaderboard claims.

## Scientific environment

Use **Python 3.12** and [the hash-locked portable profile](../locks/research-py312.txt). `requirements.txt` enables hash verification and includes exact versions for direct and transitive dependencies, including Matplotlib for standalone scientific figures; the source pins are in `requirements-research.in`. The universal lock has 69 conditional entries, of which 50 apply to the verified Windows environment. [Installed-environment and JUnit attestation](../reports/research_environment_validation.json). CI independently exercises Linux chemistry and fresh CPU neural/ranker training contracts.

This portable profile is distinct from Kaggle's execution image. Each executed notebook records its actual versions; the new fresh-training study additionally records installed package versions. A matching dependency profile is not a promise of bitwise-identical CUDA results across hardware. Synthetic-test success remains separate from dataset-backed experimental evidence.

```bash
python -m pip install --require-hashes -r requirements.txt
python scripts/check_research_lock.py
python -m pytest -q
```

Create and activate an isolated Python 3.12 virtual environment before these commands. To regenerate the lock deliberately: `uv pip compile requirements-research.in --python-version 3.12 --universal --generate-hashes --output-file locks/research-py312.txt`. Revalidate the resulting environment before replacing the attestation. Do not install the portable profile over an active Kaggle image merely to reproduce a historical notebook.

## Numerical rank and coverage contracts

Current numerical summaries and paired comparisons require nonempty integer rank arrays in 0–25 (`0` is a miss); out-of-range, floating, nonfinite and malformed inputs fail closed. This hardening does not rewrite completed experiment snapshots. `python scripts/build_coverage_report.py --tables-only` deterministically reconstructs aggregate-only candidate/ranking loss accounting after checking source/artifact hashes. Omit the flag in the locked environment to regenerate its PNG/SVG/PDF figures. Conditional coverage metrics preserve their explicit subset denominator and carry no reconstructed uncertainty claim.

## Historical source snapshots

Completed experiment bundles retain the source actually executed, embedded in their notebooks. `research_snapshot.py` parses that source inventory as literal data, without executing notebook code. Evidence validation binds results to the embedded snapshot and authenticated notebook hash rather than incorrectly requiring today's main branch to equal historical code. New experiments use distinct bundles and publication records; the completed v16 study is not rebuilt or relabeled when current code advances.

## Fresh training and family holdout

The [predeclared advanced protocol](ADVANCED_INTERNAL_STUDY.md) separates supervised query families before fitting any fresh weights. Build its distinct notebook with `python scripts/build_advanced_notebook.py` and authenticate a completed run with `python scripts/check_advanced_publication.py`. No archived models are attached. Execution status is established by the publication readback, not by the existence of these commands or the protocol.

## Kaggle publications

Run `python scripts/build_research_publications.py` to generate the two owned notebook publication bundles from the tracked notebook builders. [Publication operations](../kaggle/README.md) covers source hashes, remote status and readback. The source builders contain the editorial corrections and NEXUS chart defaults; generated notebook outputs are cleared before publication.

## Full model rebuild

Use licensed competition files and the matching public structure releases, then follow [the original architecture](../ARCHITECTURE.md), [fingerprint provenance](KAGGLE_ACCESS_AND_FP_PROVENANCE.md), [snapshot inventory](INPUT_SNAPSHOT_MANIFEST.md) and [fresh-index protocol](D2_FRESH_INDEX_EVIDENCE.md). Stage artifacts under new names. Check the input schema with `scripts/kaggle_submission_readiness.py` and audit ranked candidates with `scripts/audit_submission.py` before any separately chosen competition submission.

No competition submission is part of this update. Both **frozen-model retrospective evaluation and fresh supervised internal-family training have completed**. The explanatory publications and overlap audit use CPU; the model studies used free Kaggle GPU allocation with internet disabled. The static dashboard has zero infrastructure spend. Live BigQuery scans remain optional and unexecuted. Rebuild verified tables with `python scripts/build_advanced_report.py --tables-only`; omit the flag in the locked scientific environment to export standalone PNG/SVG/PDF figures. Plot files are checked visually as generated artifacts; no PC-screen interaction is required.

## Reproduce the completed model/control study

See [the protocol](CORRECTED_MODEL_REVALIDATION.md) and [verified results](../reports/revalidation_review.md). Generate the scratch-only bundle with `python scripts/build_revalidation_notebook.py`; authenticate its completed public execution with `python scripts/check_revalidation_publication.py`; render the report with `python scripts/build_revalidation_report.py`. The study rebuilds inputs, quarantines conflicting/unusable candidates explicitly, retains all eligible query molecules and records data/model/code hashes and runtime versions. Its approximately 41-minute ranking/control phase is a recorded duration for that run, not a runtime guarantee. Source rebuild time is additional.

## Reproduce the executed public-example content audit

Install PyArrow (tested with Kaggle 23.0.1 and local/CI 24.0.0), then run `scripts/audit_public_overlap.py` on already-authorized train/public inputs. It streams Arrow batches, exports only aggregates and hashes, and refuses to overwrite reports. Alternatively build the internet-disabled CPU bundle with `scripts/build_overlap_notebook.py`; its completed Kaggle readback is checked by `scripts/check_overlap_publication.py`. [Executed counts and methodology](PUBLIC_EXAMPLE_OVERLAP.md).
