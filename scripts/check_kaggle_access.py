"""Read-only Kaggle CASMI26 access preflight. No downloads or submissions.

Credential input: KAGGLE_API_TOKEN in the *local process environment*.
Never pass credentials as CLI arguments or write them to disk. This program
cannot accept competition rules on your behalf.

Examples:
  python scripts/check_kaggle_access.py
  python scripts/check_kaggle_access.py --check-remote

--check-remote invokes only 'kaggle competitions files -c ...'; no dataset
content is downloaded and no predictions are uploaded.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess

SLUG = "enveda-CASMI26-molecule-id-mass-spectra"


def access_preflight(*, env=None, remote=False, which=None, runner=None) -> dict:
    """Return sanitized status: never include credential values or raw CLI errors."""
    env = os.environ if env is None else env
    which = shutil.which if which is None else which
    runner = subprocess.run if runner is None else runner
    configured = bool(str(env.get("KAGGLE_API_TOKEN", "")).strip())
    available = bool(which("kaggle"))
    report = {
        "competition": SLUG,
        "token_in_environment": configured,
        "kaggle_cli_installed": available,
        "rules": "Must be accepted in the Kaggle account before private data access",
        "downloaded_files": 0,
        "submitted_predictions": 0,
    }
    if not remote:
        report["status"] = "LOCAL_CHECK_ONLY"
        return report
    if not configured:
        report["status"] = "MISSING_ENV_TOKEN"
        return report
    if not available:
        report["status"] = "KAGGLE_CLI_MISSING"
        return report
    try:
        result = runner(
            ["kaggle", "competitions", "files", "-c", SLUG],
            env=dict(env),
            capture_output=True, text=True, timeout=60, check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        report["status"] = "REMOTE_CHECK_FAILED"
        return report
    if result.returncode == 0:
        report["status"] = "FILES_LISTED"
    else:
        # Deliberately never print CLI stderr; Kaggle can reject a token,
        # ask for rules acceptance, hit a rate limit, or be offline.
        report["status"] = "REMOTE_ACCESS_NOT_CONFIRMED"
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check-remote", action="store_true",
        help="read-only Kaggle API file listing; must already accept rules",
    )
    args = ap.parse_args()
    report = access_preflight(remote=args.check_remote)
    print(json.dumps(report, indent=2))
    if args.check_remote and report["status"] != "FILES_LISTED":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
