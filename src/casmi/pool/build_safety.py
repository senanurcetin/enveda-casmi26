"""Exclusive publication for a newly built immutable CASMI candidate pool.

A preflight `exists()` check does not guard against a second writer racing
to publish its file. Stage the Parquet on the same filesystem and atomically
create a hard link at the destination; unlike os.replace, os.link refuses
an existing destination. Clean the staging file on all outcomes.

Filesystems without hard-link support fail closed. Readers must use immutable
versioned pool tags rather than expecting an existing file to be replaced.
"""
from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
import tempfile


@contextmanager
def fresh_pool_file(destination):
    """Yield a temporary Parquet path; publish it only on successful exit."""
    destination = Path(destination)
    if not destination.name or destination.name in (".", ".."):
        raise ValueError("An explicit non-root pool output path is required")
    if os.path.lexists(destination):
        raise FileExistsError("Refusing to overwrite frozen candidate pool")
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(
        prefix=f".{destination.name}-building-", suffix=".parquet",
        dir=destination.parent,
    )
    os.close(fd)
    stage = Path(name)
    try:
        yield stage
        if not stage.is_file() or stage.stat().st_size == 0:
            raise ValueError("Candidate pool build produced an empty output file")
        # Hard-link creation, unlike rename or os.replace, fails if a second
        # process has created the final path during the build.
        os.link(stage, destination, follow_symlinks=False)
    finally:
        stage.unlink(missing_ok=True)
