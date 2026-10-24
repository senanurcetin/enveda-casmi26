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

Use Python 3.12 and the pinned versions in `.github/workflows/nexus-research-checks.yml` for official-metric regression. `requirements.txt` pins RDKit 2026.3.3; the remaining full-pipeline dependencies use lower bounds, so it is not a complete environment lock. Capture the environment with each real-data run. Do not interpret synthetic-test success as a reproduced competition result.

```bash
python -m pip install -r requirements.txt pytest
python -m pytest -q
```

## Kaggle publications

Run `python scripts/build_research_publications.py` to generate the two owned notebook publication bundles from the tracked notebook builders. [Publication operations](../kaggle/README.md) covers source hashes, remote status and readback. The source builders contain the editorial corrections and NEXUS chart defaults; generated notebook outputs are cleared before publication.

## Full model rebuild

Use licensed competition files and the matching public structure releases, then follow [the original architecture](../ARCHITECTURE.md), [fingerprint provenance](KAGGLE_ACCESS_AND_FP_PROVENANCE.md), [snapshot inventory](INPUT_SNAPSHOT_MANIFEST.md) and [fresh-index protocol](D2_FRESH_INDEX_EVIDENCE.md). Stage artifacts under new names. Check the input schema with `scripts/kaggle_submission_readiness.py` and audit ranked candidates with `scripts/audit_submission.py` before any separately chosen competition submission.

No new competition submission is part of this portfolio maintenance. Full GPU training and corrected-protocol evaluation remain pending. The static dashboard has zero infrastructure spend; hosted Kaggle execution is CPU-only for these publication bundles and subject to account quotas. Optional BigQuery scans and GPU retraining are not part of the default workflow.

## Reproduce the executed public-example content audit

Install PyArrow (tested with Kaggle 23.0.1 and local/CI 24.0.0), then run `scripts/audit_public_overlap.py` on already-authorized train/public inputs. It streams Arrow batches, exports only aggregates and hashes, and refuses to overwrite reports. Alternatively build the internet-disabled CPU bundle with `scripts/build_overlap_notebook.py`; its completed Kaggle readback is checked by `scripts/check_overlap_publication.py`. [Executed counts and methodology](PUBLIC_EXAMPLE_OVERLAP.md).
