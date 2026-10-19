# A1 spectrum-store immutability and staged build protocol

**Open scientific research only (2026-10-19).** This change protects the
integrity of the Enveda CASMI26 A1 spectrum store. It does not create a
product, inference API, commercial service or new competition score.

## Failure mode discovered

The previous `scripts/build_store.py` used

```python
out = ARTIFACT_DIR / "store"
out.mkdir(parents=True, exist_ok=True)
f_mz = open(out / "mz.f32", "wb")
```

which **silently truncated existing peak files** before writing a new
store. A user attempting to rerun a build on an established research
workspace could destroy the frozen input evidence of historical model
results. An exception halfway through would leave partially written files
in the same location expected by the validation pipeline.

## Protection now implemented

- Build into a unique **staging directory on the same filesystem**,
  never inside the frozen existing `artifacts/store`.
- Verify offsets, file sizes and the `n_peaks` metadata using the existing
  post-build A1 integrity check **before publication**.
- Publish the completed directory with a directory rename only on
  successful exit; clean up staged files on raised exceptions.
- **Never overwrite** an existing destination directory, even if it is
  empty. Existing dangling symlinks are also rejected.
- An explicitly supplied `--output-dir` allows an independently named
  experimental store; `--train` selects the authorized source Parquet.
  The default target remains `ARTIFACT_DIR/store` but must not already
  exist.

### Commands

```bash
# First build, only when the default output does not already exist.
python scripts/build_store.py

# Explicitly version a new scientific build without overwriting the
# frozen original (use a permitted local source file).
python scripts/build_store.py \
  --train /approved/casmi/train.parquet \
  --output-dir /approved/casmi/artifacts/store_research_v2
```

If a target directory already exists, the command fails immediately and
asks for another new target. Do **not** delete the old store merely to make
the command pass: retain its source hash and prior experiment references.

This use of filesystem rename provides an atomic directory publication
on supported same-filesystem environments. It is **not a distributed
locking mechanism**; concurrent writers need an external coordination
lock. It is also not proof of chemical correctness or independent source
provenance.

## Testing and evidence

Synthetic standard-library tests cover successful one-time publication,
refusal to overwrite frozen files, cleanup after failures, concurrent-target
appearance, dangling symlinks and nested output paths.

A local focused check exercised **four staging and non-overwrite assertions**.
Per the user's instruction **no GitHub CI was launched**; the full Parquet
build with the original CASMI dataset was not executed.

Historical public MRR@25 `0.422` (2026-10-08) is unchanged. Any experiment
using a separately built store needs its own SHA-256 manifest, fresh
features/caches and identically sampled independent validation before scores
are compared.
