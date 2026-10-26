"""Validate hashed exact pins and optionally attest an isolated installed environment."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import re
import sys
import xml.etree.ElementTree as ET
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def lock_entries(path):
    from packaging.requirements import Requirement
    entries, chunks = [], []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        chunks.append(value.removesuffix("\\").strip())
        if value.endswith("\\"):
            continue
        parts = " ".join(chunks).split(" --hash=sha256:")
        chunks = []
        req = Requirement(parts[0])
        pins = list(req.specifier)
        if req.url or len(pins) != 1 or pins[0].operator != "==" or "*" in pins[0].version:
            raise ValueError("Every locked requirement must have one exact version")
        if len(parts) < 2 or any(not re.fullmatch(r"[0-9a-f]{64}", digest) for digest in parts[1:]):
            raise ValueError("Every requirement needs valid SHA-256 distribution hashes")
        entries.append(req)
    if chunks or not entries:
        raise ValueError("Empty or truncated lockfile")
    return entries


def check_installed(entries, *, environment=None, version_of=None):
    from packaging.markers import default_environment
    from packaging.utils import canonicalize_name
    environment = default_environment() if environment is None else environment
    version_of = importlib.metadata.version if version_of is None else version_of
    selected = {}
    for req in entries:
        if req.marker and not req.marker.evaluate(environment):
            continue
        name = canonicalize_name(req.name)
        version = version_of(name)
        if not req.specifier.contains(version, prereleases=True):
            raise ValueError("Installed package differs from lock: " + name)
        if name in selected and selected[name] != version:
            raise ValueError("Overlapping conditional requirement pins")
        selected[name] = version
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure-only", action="store_true")
    parser.add_argument("--junit", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    path = ROOT / "locks/research-py312.txt"
    entries = lock_entries(path)
    from packaging.requirements import Requirement
    pinned = {req.name.lower(): str(req.specifier) for req in entries}
    for line in (ROOT / "requirements-research.in").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            req = Requirement(line)
            if pinned.get(req.name.lower()) != str(req.specifier):
                raise ValueError("Direct requirements and transitive lock differ")
    result = {"status": "HASHED_EXACT_PINS_VERIFIED", "profile": "portable Python 3.12 research environment",
              "checked_at": datetime.now(UTC).isoformat(), "conditional_locked_entries": len(entries),
              "lock_sha256": hashlib.sha256(path.read_text(encoding="utf-8").encode()).hexdigest()}
    if not args.structure_only:
        if sys.version_info[:2] != (3, 12):
            raise ValueError("Installed-profile attestation requires Python 3.12")
        result.update({"status": "INSTALLED_HASHED_PROFILE_VERIFIED", "installed_packages": check_installed(entries),
                       "python": platform.python_version(), "system": platform.system(), "machine": platform.machine()})
    if args.junit:
        suites = ET.parse(args.junit).getroot()
        suite = suites if suites.tag == "testsuite" else suites.find("testsuite")
        if suite is None or int(suite.attrib.get("failures", 0)) or int(suite.attrib.get("errors", 0)):
            raise ValueError("JUnit evidence is missing or contains failed tests")
        result["junit"] = {field: suite.attrib[field] for field in ("tests", "failures", "errors", "skipped", "time")}
        result["junit_sha256"] = hashlib.sha256(args.junit.read_bytes()).hexdigest()
        result["junit_scope"] = "Pytest JUnit cases may include subtests; not a line/branch coverage metric"
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"{result['status']}: {len(entries)} conditional pins; {len(result.get('installed_packages', {}))} installed packages attested")


if __name__ == "__main__":
    main()
