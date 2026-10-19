"""Synthetic CLI comparison guards for missing/incomplete input provenance.

No competition data or scientific performance results are generated.
"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from casmi.validate.source_manifest import (  # noqa: E402
    MANIFEST_SCHEMA, required_validation_inputs,
)
from compare_runs import main  # noqa: E402


def fixture(ranked):
    n = 18
    return {
        "tag": "synthetic", "n_queries": n,
        "query_ids": [f"Q{i}" for i in range(n)],
        "query_signatures": [format(i + 1, "064x") for i in range(n)],
        "query_signature_version": "casmi-query-spectra-v2",
        "classes": {
            cls: {"ranks": [1 if ranked else 0] * n}
            for cls in ("C1", "C2", "C3")
        },
    }


def signed(run):
    run["source_manifest_sha256"] = "e" * 64
    run["source_manifest_schema"] = MANIFEST_SCHEMA
    run["source_manifest_pool_tag"] = ""
    run["source_manifest_coverage"] = "required-validation-inputs-v1"
    run["source_manifest_required_inputs"] = list(
        required_validation_inputs(analog=False)
    )
    run["source_manifest_file_count"] = 5
    run["config"] = {"channels": ["d1"]}
    return run


def compare(a, b):
    with tempfile.TemporaryDirectory() as directory:
        one, two = Path(directory) / "A.json", Path(directory) / "B.json"
        one.write_text(json.dumps(a), encoding="utf-8")
        two.write_text(json.dumps(b), encoding="utf-8")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            decision = main(one, two, n_boot=100, seed=3)
        return decision, out.getvalue()


class ModelComparisonManifestTests(unittest.TestCase):
    def test_reject_unverified_source_even_with_artificial_perfect_effect(self):
        verdict, output = compare(fixture(False), fixture(True))
        self.assertEqual(verdict, "UNVERIFIED")
        self.assertIn("QUERY-PAIRED ONLY", output)
        self.assertIn("VERDICT: UNVERIFIED", output)

    def test_matching_complete_inventory_enables_paired_verdict(self):
        verdict, output = compare(signed(fixture(False)), signed(fixture(True)))
        self.assertEqual(verdict, "ACCEPT")
        self.assertIn("VERIFIED query identities", output)
        self.assertIn("VERDICT: ACCEPT", output)

    def test_legacy_manifest_hash_without_coverage_is_only_exploratory(self):
        a, b = fixture(False), fixture(True)
        for record in (a, b):
            record["source_manifest_sha256"] = "e" * 64
        verdict, output = compare(a, b)
        self.assertEqual(verdict, "UNVERIFIED")
        self.assertIn("VERDICT: UNVERIFIED", output)

    def test_changed_source_digests_fail_closed(self):
        a, b = signed(fixture(False)), signed(fixture(True))
        b["source_manifest_sha256"] = "f" * 64
        with self.assertRaisesRegex(ValueError, "snapshot mismatch"):
            compare(a, b)


if __name__ == "__main__":
    unittest.main()
