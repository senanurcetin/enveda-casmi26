"""End-to-end paired bootstrap behavior on entirely synthetic CASMI run JSON.

Synthetic examples are unit fixtures, NOT empirical competition results.
Requires numpy, but neither RDKit nor Kaggle data.
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/compare_runs.py"


def fixtures(count=20):
    ids = [f"m{i:02d}" for i in range(count)]
    signatures = [f"{i:064x}" for i in range(count)]
    baseline = {
        "tag": "control", "n_queries": count, "query_ids": ids,
        "query_signatures": signatures,
        "classes": {
            cls: {"ranks": [0] * count, "mrr": 0.0}
            for cls in ("C1", "C2", "C3")
        },
    }
    experiment = {
        "tag": "candidate", "n_queries": count,
        "query_ids": list(reversed(ids)),
        "query_signatures": list(reversed(signatures)),
        "classes": {
            "C1": {"ranks": [1] * count, "mrr": 1.0},
            "C2": {"ranks": [1] * count, "mrr": 1.0},
            "C3": {"ranks": [0] * count, "mrr": 0.0},
        },
    }
    return baseline, experiment


def execute(a, b, extra=()):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        pa, pb = root / "a.json", root / "b.json"
        pa.write_text(json.dumps(a), encoding="utf-8")
        pb.write_text(json.dumps(b), encoding="utf-8")
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(pa), str(pb),
             "--bootstrap", "200", *extra],
            capture_output=True, text=True,
        )


class PairedBootstrapIntegrationTests(unittest.TestCase):
    def test_synthetic_signal_and_reverse_order(self):
        result = execute(*fixtures())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("VERIFIED paired identities", result.stdout)
        self.assertIn("VERDICT: ACCEPT", result.stdout)

    def test_same_count_changed_spectrum_rejected(self):
        baseline, candidate = fixtures()
        candidate["query_signatures"][1] = "f" * 64
        result = execute(baseline, candidate)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Spectrum/knockout mismatch", result.stderr)

    def test_legacy_override_never_produces_acceptance(self):
        baseline, candidate = fixtures()
        for report in (baseline, candidate):
            report.pop("query_ids")
            report.pop("query_signatures")
        result = execute(baseline, candidate)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("legacy runs", result.stderr)
        result = execute(baseline, candidate, ("--allow-unverified-legacy-order",))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("VERDICT: UNVERIFIED", result.stdout)
        self.assertNotIn("VERDICT: ACCEPT", result.stdout)


if __name__ == "__main__":
    unittest.main()
