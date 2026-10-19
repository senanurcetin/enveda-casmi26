"""Read-only CASMI candidate Morgan fingerprint numeric quality audit.

Example:
  python scripts/audit_pool_fingerprints.py \
    --output runs/parent_v2_fingerprint_quality.json

Use CASMI_ARTIFACTS and CASMI_POOL_TAG to point at the intended immutable
research dataset. No Kaggle access, secret handling, remote calls, downloads
or modifications to candidate artifacts.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import (  # noqa: E402
    ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE, POOL_TAG,
)
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.pool.fingerprint_audit import fingerprint_quality_summary  # noqa: E402
from casmi.pool.fingerprint_provenance import verify_fingerprint_provenance  # noqa: E402


def audit(artifacts=ARTIFACT_DIR, *, pool_file=POOL_FILE,
          fingerprint_file=POOL_FP_FILE, key_file=POOL_FP_KEYS_FILE,
          pool_tag=POOL_TAG, batch_rows=65536) -> dict:
    """Check the inference pool order and the numerical fingerprint content."""
    artifacts = Path(artifacts)
    pool = CandidatePool(pd.read_parquet(
        artifacts / pool_file, columns=["pool_key", "smiles", "mass"]
    )).table
    keys = pd.read_parquet(artifacts / key_file, columns=["pool_key"]).pool_key.to_numpy()
    vectors = np.load(artifacts / fingerprint_file, mmap_mode="r", allow_pickle=False)
    result = fingerprint_quality_summary(
        pool.pool_key.to_numpy(), keys, vectors, batch_rows=batch_rows
    )
    provenance = artifacts / f"pool_fp_provenance{pool_tag}.json"
    if provenance.exists():
        with provenance.open("r", encoding="utf-8") as stream:
            manifest = json.load(stream)
        proof = verify_fingerprint_provenance(pool, keys, vectors, manifest)
        result["provenance"] = {
            "status": proof["status"], "schema": proof["schema"],
        }
    else:
        result["provenance"] = {
            "status": "UNVERIFIED_LEGACY_FINGERPRINTS",
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--batch-rows", type=int, default=65536)
    args = parser.parse_args()
    try:
        report = audit(batch_rows=args.batch_rows)
        rendered = json.dumps(report, indent=2) + "\n"
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x", encoding="utf-8") as file:
                file.write(rendered)
        print(rendered, end="")
        if report["status"] != "PASS_NONZERO_FINGERPRINTS":
            raise SystemExit(2)
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
