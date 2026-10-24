"""Validate the public research ledger and export a source-hashed review report."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATES = {"historical", "implemented", "measured", "pending"}


def validate(ledger: dict, root: Path) -> None:
    if ledger.get("schema_version") != 1:
        raise ValueError("Unsupported evidence schema")
    seen = set()
    for collection in ("experiments", "claims"):
        if not ledger.get(collection):
            raise ValueError("Empty evidence collection: " + collection)
        for row in ledger[collection]:
            if not row.get("id") or row["id"] in seen:
                raise ValueError("Missing or duplicate evidence ID")
            seen.add(row["id"])
            source = (root / row["source"]).resolve()
            if not source.is_relative_to(root.resolve()) or not source.is_file():
                raise ValueError("Evidence source must be an existing repository file")
            if collection == "claims" and row.get("state") not in STATES:
                raise ValueError("Unknown evidence state")
            if collection == "experiments":
                score = row.get("score")
                if isinstance(score, bool) or not isinstance(score, (int, float)) or not 0 <= score <= 1:
                    raise ValueError("MRR must be finite and between zero and one")
                if row.get("family") not in {"own", "external", "fusion"}:
                    raise ValueError("Unknown model attribution")


def report(ledger: dict, root: Path) -> str:
    validate(ledger, root)
    lines = ["# CASMI26 research review", "", f"Reviewed: {ledger['reviewed_on']}. Score snapshot: {ledger['score_snapshot']}.",
             "", ledger["scope"] + ".", "", "These are source-file integrity hashes, not external verification of the claims.",
             "", "| Experiment | Public MRR@25 | Attribution |", "|---|---:|---|"]
    for row in ledger["experiments"]:
        lines.append(f"| {row['label']} | {row['score']:.3f} | {row['family']} |")
    for row in ledger["claims"]:
        # Git can check out CRLF on Windows and LF on Linux. Hash canonical text.
        canonical = (root / row["source"]).read_text(encoding="utf-8").encode("utf-8")
        digest = hashlib.sha256(canonical).hexdigest()
        lines.extend(["", f"## {row['title']} — {row['state']}", "", row["detail"], "",
                      f"Source: [{row['source']}](../{row['source']}). SHA-256 (UTF-8 / LF): `{digest}`."])
    return "\n".join(lines) + "\n"


def main() -> None:
    ledger = json.loads((ROOT / "reports/research_evidence.json").read_text(encoding="utf-8"))
    output = ROOT / "reports/research_review.md"
    output.write_text(report(ledger, ROOT), encoding="utf-8")
    print("Evidence ledger valid; wrote reports/research_review.md")


if __name__ == "__main__":
    main()
