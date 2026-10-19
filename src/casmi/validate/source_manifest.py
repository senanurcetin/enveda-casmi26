"""Privacy-preserving SHA-256 manifests for approved local CASMI artifacts.

No dataset payload, raw spectra, private labels, credentials or absolute paths
are written into a manifest. This is a provenance check, not a security
signature or an assertion that upstream Kaggle/BigQuery data are correct.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath

MANIFEST_SCHEMA = "casmi-local-artifact-manifest-v1"
CHUNK_BYTES = 1024 * 1024


def _safe_relative(root: Path, name: str) -> Path:
    if not isinstance(name, str) or not name or "\\" in name:
        raise ValueError("Manifest entries must be normalized relative POSIX paths")
    normalized = PurePosixPath(name)
    if (normalized.is_absolute() or name != normalized.as_posix()
            or any(part in (".", "..") for part in name.split("/"))):
        raise ValueError(f"Unsafe artifact path: {name!r}")
    path = (root / name).resolve(strict=True)
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Artifact resolves outside source root: {name!r}")
    if not path.is_file():
        raise ValueError(f"Expected an artifact file: {name!r}")
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_files_digest(files: list[dict]) -> str:
    payload = json.dumps(files, ensure_ascii=True, sort_keys=True,
                         separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def create_manifest(root: Path, relative_paths) -> dict:
    """Hash existing explicit input artifacts; no directory glob or traversal."""
    root = Path(root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("Source root must be a directory")
    if not isinstance(relative_paths, (list, tuple)) or not relative_paths:
        raise ValueError("Provide at least one explicit input-artifact path")
    if len(set(relative_paths)) != len(relative_paths):
        raise ValueError("Duplicate artifact manifest path")
    files = []
    for name in sorted(relative_paths):
        path = _safe_relative(root, name)
        files.append({
            "path": name, "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    return {
        "schema": MANIFEST_SCHEMA,
        "method": "streaming SHA-256 per local input file",
        "files": files,
        "files_digest_sha256": _canonical_files_digest(files),
        "source_assertion": (
            "local byte identity only; no independent Kaggle/BigQuery validation"
        ),
    }


def verify_manifest(root: Path, manifest: dict) -> dict:
    """Fail closed when a local artifact differs from the recorded snapshot."""
    if not isinstance(manifest, dict) or manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("Unsupported artifact manifest schema")
    entries = manifest.get("files")
    if not isinstance(entries, list) or not entries:
        raise ValueError("Manifest must contain artifact entries")
    paths = [x.get("path") for x in entries if isinstance(x, dict)]
    if len(paths) != len(entries) or len(paths) != len(set(paths)):
        raise ValueError("Invalid or duplicate manifest paths")
    for row in entries:
        if type(row.get("bytes")) is not int or row["bytes"] < 0:
            raise ValueError("Invalid manifest file size")
        sha = row.get("sha256")
        if not isinstance(sha, str) or len(sha) != 64 or any(c not in "0123456789abcdef" for c in sha):
            raise ValueError("Invalid manifest file hash")
    if entries != sorted(entries, key=lambda e: e["path"]):
        raise ValueError("Manifest entries must be sorted")
    digest = manifest.get("files_digest_sha256")
    if digest != _canonical_files_digest(entries):
        raise ValueError("Manifest file list/digest mismatch")
    fresh = create_manifest(root, paths)
    if fresh["files"] != entries:
        raise ValueError("Local input artifacts differ from the signed-inventory snapshot")
    return {
        "status": "VERIFIED_LOCAL_ARTIFACTS",
        "schema": MANIFEST_SCHEMA,
        "files_digest_sha256": fresh["files_digest_sha256"],
        "files_checked": len(entries),
    }
