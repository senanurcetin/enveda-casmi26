# Reproduce and review

Run commands from the repository root. The project has no hosted model service or paid infrastructure requirement. Data-dependent rebuilds require authorized competition files and the declared public source releases; a source bundle alone does not grant data access or free compute capacity.

## Inspect without scientific dependencies

```bash
python -m http.server 8000
# Open http://localhost:8000/site/
```

Serve the repository root so the static NEXUS dashboard can read committed reports. It loads no model weights and calls no Kaggle, BigQuery or private data API.

For the dependency-free checks, put both the repository root and `src` on `PYTHONPATH` (use `;` as the separator on Windows and `:` on POSIX), then run:

```bash
python -m unittest discover -s tests -p '*_stdlib.py' -v
python scripts/research_evidence.py
python scripts/check_documentation_links.py
```

The ledger command validates evidence references and regenerates source-file integrity hashes. The link checker validates repository-local destinations; it does not check remote availability or fragment anchors. Neither tool authenticates historical leaderboard scores.

## Locked scientific environment

Create and activate an isolated Python 3.12 environment, then run:

```bash
python -m pip install --require-hashes -r requirements.txt
python scripts/check_research_lock.py
python -m pytest tests -q
```

The [universal lock](../../locks/research-py312.txt) contains 69 conditional exact/hash-pinned requirements; 50 apply to the recorded Windows environment. [Dated installed-environment and JUnit attestation](../../reports/research_environment_validation.json). Kaggle uses a separate image: each experiment records its actual packages, hardware mode and attachments. Local test success is not dataset-backed model evidence or a promise of bitwise reproduction across platforms.

To change the lock deliberately, use `uv pip compile requirements-research.in --python-version 3.12 --universal --generate-hashes --output-file locks/research-py312.txt`, then validate the installed environment and complete tests before replacing its attestation. Do not overwrite a running Kaggle environment with the local profile.

## Regenerate committed reports

```bash
python scripts/research_evidence.py
python scripts/build_revalidation_report.py
python scripts/build_advanced_report.py --tables-only
python scripts/build_coverage_report.py --tables-only
python scripts/build_purged_external_report.py --tables-only
```

These commands validate stored aggregates and regenerate reviews; they do not fit models or supply new scientific evidence. Omit `--tables-only` in the locked plotting environment to export advanced, coverage and purged-study PNG/SVG/PDF figures. Inspect exported figures as files. Font/platform differences are not covered by deterministic table checks.

## Executed notebooks and authenticated readback

[The Kaggle index](../../kaggle/README.md) maps each publication to its source bundle, builder and readback command. Builders for completed studies verify frozen source instead of incorporating current code. `research_snapshot.py` reads embedded Python/JSON inventories as literal data without executing them. Do not modify or repush completed evidence to match today's code.

Authenticated `check_*publication*.py` commands inspect public completion, source/attachments and bounded aggregate artifacts. They update readback records and can refresh the stored aggregate from the authenticated publication; they do not start a notebook or submit predictions. Use existing local credentials and inspect changes before committing refreshed records.

The retrospective and fresh internal model studies used free Kaggle GPU execution. The separate purged-reference experiment used CPU; source/content audits and explanatory publications also use CPU. Runtime versions and inputs belong to their specific executions, not a single shared model environment.

## Full-data execution and historical recipes

Use the [architecture](../../ARCHITECTURE.md), [data contracts](data-integrity.md), [evaluation guide](evaluation.md) and the relevant [study design](../README.md#completed-studies). Record byte hashes, create new immutable artifact destinations, verify candidate fingerprints, and rebuild D2 references. The purged study additionally reconstructs historical exclusions before filtering the full reference corpus; its exact code and inputs are embedded in its completed bundle.

Completed tests must not become tuning sets. A new execution requires its own declared population/design, input permissions and available free compute. Original-snapshot novelty failure cannot be repaired by relabeling the completed purged-reference result.

`scripts/legacy/` contains archived machine-specific shell chains, including old publishing/submission automation. Root-level training and diagnostic helpers retain stable paths used by other code. These are not the default reproduction workflow. `runs/` is ignored local scratch space; archived CSV summaries are in `reports/historical/`. No shell chain, competition submission or live warehouse scan runs during setup or documentation checks.
