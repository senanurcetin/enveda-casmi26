"""Exclusive, staged publication of immutable CASMI spectrum-store artifacts.

The A1 builder must never open an existing spectrum store for writing.
Writing four related output files directly into the public experiment's
historical store can silently truncate valid metadata/peak snapshots if
a build fails midway. This context manager stages a complete new directory
on the same filesystem and publishes it only after the caller's verification.

No external dependencies or CI environment required.
"""
from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
import shutil
import tempfile


@contextmanager
def fresh_store_target(destination):
    """Yield a new staging directory; publish once only on successful exit.

    Existing destinations, including empty directories and dangling symlinks,
    are refused. The staging directory is removed on errors. The caller must
    complete and verify every artifact before leaving the context.
    The same-file-system rename makes publication atomic on supported OSes;
    parallel writers should be managed externally with a project lock.
    """
    destination = Path(destination)
    if not destination.name or destination.name in (".", ".."):
        raise ValueError("An explicit, non-root output directory is required")
    if os.path.lexists(destination):
        raise FileExistsError(
            f"Refusing to overwrite an existing spectrum store: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}-building-", dir=destination.parent,
    ))
    try:
        yield stage
        if os.path.lexists(destination):
            raise FileExistsError(
                f"Output appeared during build; refusing to overwrite: {destination}"
            )
        stage.rename(destination)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
