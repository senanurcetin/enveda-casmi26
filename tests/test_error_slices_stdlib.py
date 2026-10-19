"""Source-aligned, label-free error-slice regressions on synthetic run JSON."""
import json
import subprocess
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.error_slices import (  # noqa: E402
    describe_query, diagnostic_slices,
)


def fixtures():
    ids = [f"q-{i}" for i in range(8)]
    context = [
        {
            "adduct": "[M+H]+" if i < 4 else "[M-H]-",
            "ionization_mode": "positive" if i < 4 else "negative",
            "spectra_bin": "1" if i < 4 else "4+",
            "precursor_mz_bin": "<300" if i < 4 else ">=600",
        }
        for i in range(8)
    ]
    a = {
        "tag": "A", "n_queries": 8, "query_ids": ids,
        "query_signatures": [f"{i+1:064x}" for i in range(8)],
        "query_signature_version": "casmi-query-spectra-v2",
        "query_context": context,
        "classes": {
            cls: {"ranks": [0] * 8} for cls in ("C1", "C2", "C3")
        },
    }
    b = deepcopy(a)
    b["tag"] = "B"
    order = list(reversed(range(8)))
    for field in ("query_ids", "query_signatures", "query_context"):
        b[field] = [a[field][i] for i in order]
    b["classes"]["C1"]["ranks"] = [1] * 8
    b["classes"]["C2"]["ranks"] = [1 if i < 4 else 0 for i in order]
    return a, b


class SliceAnalysisTests(unittest.TestCase):
    def test_input_only_context_handles_mixed_modes_and_precursor_bins(self):
        spectra = [
            {"adduct": "[M+H]+", "mode": "positive", "precursor_mz": 250.0},
            {"adduct": "[M+Na]+", "mode": "positive", "precursor_mz": 500.0},
        ]
        self.assertEqual(describe_query(spectra), {
            "adduct": "mixed", "ionization_mode": "positive",
            "spectra_bin": "2-3", "precursor_mz_bin": "300-600",
        })
        self.assertEqual(
            describe_query([dict(spectra[0], precursor_mz=float("nan"))])
            ["precursor_mz_bin"], "unknown"
        )
        with self.assertRaisesRegex(ValueError, "at least one"):
            describe_query([])

    def test_slices_align_by_id_and_suppress_small_subgroups(self):
        a, b = fixtures()
        report = diagnostic_slices(a, b, min_n=3)
        self.assertEqual(report["status"], "DESCRIPTIVE_ONLY")
        self.assertEqual(report["source_evidence"], "QUERY_PAIRED_ONLY_SOURCE_UNVERIFIED")
        self.assertEqual(report["n_molecules"], 8)
        self.assertEqual(len(report["slices"]["adduct"]), 2)
        self.assertTrue(all(item["status"] == "DESCRIPTIVE_ONLY"
                            for item in report["slices"]["adduct"]))
        by_adduct = {s["value"]: s for s in report["slices"]["adduct"]}
        self.assertEqual(
            by_adduct["[M+H]+"]["classes"]["C2"]["candidate_mrr"], 1.0
        )
        self.assertEqual(
            by_adduct["[M-H]-"]["classes"]["C2"]["candidate_mrr"], 0.0
        )
        self.assertFalse("verdict" in json.dumps(report).lower())
        suppressed = diagnostic_slices(a, b, min_n=5)
        for record in suppressed["slices"]["adduct"]:
            self.assertEqual(record["status"], "SUPPRESSED_SMALL_CELL")
            self.assertNotIn("classes", record)

    def test_fail_closed_on_inconsistent_metadata_and_query_population(self):
        a, b = fixtures()
        b["query_context"][0]["adduct"] = "unexpected"
        with self.assertRaisesRegex(ValueError, "metadata differ"):
            diagnostic_slices(a, b)
        a, b = fixtures()
        b["query_context"].pop()
        with self.assertRaisesRegex(ValueError, "query_context"):
            diagnostic_slices(a, b)
        a, b = fixtures()
        b["query_signatures"][0] = "f" * 64
        with self.assertRaisesRegex(ValueError, "Spectrum/knockout"):
            diagnostic_slices(a, b)

    def test_control_runs_are_not_model_diagnostics(self):
        a, b = fixtures()
        b["spectral_control"] = {"mode": "blank", "seed": 0}
        with self.assertRaisesRegex(ValueError, "diagnostic only"):
            diagnostic_slices(a, b)

    def test_minimum_threshold_must_be_integer(self):
        a, b = fixtures()
        for count in (0, 1, 2.5, True):
            with self.assertRaisesRegex(ValueError, "min_n"):
                diagnostic_slices(a, b, min_n=count)

    def test_cli_outputs_only_aggregate_counts_not_query_ids(self):
        a, b = fixtures()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            pa, pb = directory / "A.json", directory / "B.json"
            out = directory / "slices.json"
            pa.write_text(json.dumps(a), encoding="utf-8")
            pb.write_text(json.dumps(b), encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts" / "diagnose_validation_slices.py"),
                str(pa), str(pb), "--min-n", "3", "--output", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads(out.read_text()), json.loads(result.stdout))
            self.assertNotIn("q-7", result.stdout)


if __name__ == "__main__":
    unittest.main()
