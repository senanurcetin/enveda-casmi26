"""Credential-safe Kaggle preflight: local checks and simulated listing only."""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_kaggle_access import access_preflight, SLUG  # noqa: E402


class KaggleAccessPreflightTests(unittest.TestCase):
    def test_local_default_never_calls_remote(self):
        def forbidden(*args, **kwargs):
            self.fail("Remote call must be explicit")
        result = access_preflight(
            env={"KAGGLE_API_TOKEN": "synthetic-secret-not-real"},
            remote=False, which=lambda name: "/bin/kaggle", runner=forbidden,
        )
        self.assertEqual(result["status"], "LOCAL_CHECK_ONLY")
        self.assertNotIn("synthetic-secret-not-real", str(result))
        self.assertEqual(result["downloaded_files"], 0)
        self.assertEqual(result["submitted_predictions"], 0)

    def test_remote_listing_does_not_pass_token_on_command_line(self):
        seen = {}
        def runner(command, **kwargs):
            seen["command"], seen["env"] = command, kwargs["env"]
            return SimpleNamespace(returncode=0, stdout="train.parquet", stderr="")
        report = access_preflight(
            env={"KAGGLE_API_TOKEN": "synthetic-secret-not-real"},
            remote=True, which=lambda name: "/bin/kaggle", runner=runner,
        )
        self.assertEqual(report["status"], "FILES_LISTED")
        self.assertEqual(seen["command"], [
            "kaggle", "competitions", "files", "-c", SLUG,
        ])
        self.assertNotIn("synthetic-secret-not-real", str(report))
        self.assertNotIn("synthetic-secret-not-real", " ".join(seen["command"]))

    def test_missing_auth_or_cli_reports_sanitized_failure(self):
        self.assertEqual(access_preflight(
            env={}, remote=True, which=lambda x: "/bin/kaggle",
        )["status"], "MISSING_ENV_TOKEN")
        self.assertEqual(access_preflight(
            env={"KAGGLE_API_TOKEN": "synthetic"}, remote=True, which=lambda x: None,
        )["status"], "KAGGLE_CLI_MISSING")

    def test_remote_rejection_never_emits_stderr(self):
        def runner(*args, **kwargs):
            return SimpleNamespace(
                returncode=1, stdout="",
                stderr="Authentication rejected: synthetic-secret-not-real",
            )
        result = access_preflight(
            env={"KAGGLE_API_TOKEN": "synthetic-secret-not-real"},
            remote=True, which=lambda x: "/bin/kaggle", runner=runner,
        )
        self.assertEqual(result["status"], "REMOTE_ACCESS_NOT_CONFIRMED")
        self.assertNotIn("synthetic-secret-not-real", str(result))


if __name__ == "__main__":
    unittest.main()
