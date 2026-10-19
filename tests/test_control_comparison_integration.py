"""Synthetic identity and paired-uncertainty tests for MS/MS controls."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from casmi.validate.control_alignment import align_control_runs  # noqa: E402
from casmi.validate.run_alignment import SIGNATURE_VERSION, align_runs  # noqa: E402
from compare_control_runs import diagnostic_summary  # noqa: E402


def baseline():
    ids = ["lib:A", "lib:B", "lib:C", "lib:D"]
    ranks = {"C1": [1, 2, 1, 0], "C2": [2, 1, 0, 4],
             "C3": [0, 0, 0, 0]}
    return {
        "tag": "reference",
        "n_queries": len(ids),
        "query_ids": ids,
        "query_signature_version": SIGNATURE_VERSION,
        "query_signatures": ["a"*64, "b"*64, "c"*64, "d"*64],
        "spectral_control": {"mode": "none", "seed": 0},
        "query_sampling": {
            "library": "enveda-np-examples", "max_spectra": 5,
            "max_queries": 4, "seed": 0,
            "spectral_control": "none", "control_seed": 0,
        },
        "config": {"channels": ["d1", "d2"], "analog": {"reference_policy": "legacy"}},
        "classes": {name: {"ranks": r} for name, r in ranks.items()},
    }


def control(ref, mode="blank"):
    order = [2, 0, 3, 1]
    report = deepcopy(ref)
    report["tag"] = "control"
    report["query_ids"] = [ref["query_ids"][i] for i in order]
    report["query_signatures"] = ["1"*64, "2"*64, "3"*64, "4"*64]
    report["spectral_control"] = {
        "mode": mode, "seed": 17,
        "source_query_signatures": [ref["query_signatures"][i] for i in order],
    }
    report["query_sampling"]["spectral_control"] = mode
    report["query_sampling"]["control_seed"] = 17
    for c in ("C1", "C2", "C3"):
        report["classes"][c]["ranks"] = [
            ref["classes"][c]["ranks"][i] for i in order
        ]
    report["classes"]["C1"]["ranks"] = [0] * 4
    return report


class PairedNegativeControlTests(unittest.TestCase):
    def test_reordered_control_ids_realign_by_actual_source_fingerprints(self):
        a = baseline()
        b = control(a)
        paired = align_control_runs(a, b)
        self.assertEqual(paired["n_queries"], 4)
        self.assertEqual(paired["changed_queries"], 4)
        self.assertEqual(paired["reference_ranks"]["C1"], [1, 2, 1, 0])
        self.assertEqual(paired["treatment_ranks"]["C1"], [0, 0, 0, 0])
        self.assertTrue(paired["verified_source_pairing"])

    def test_report_produces_no_model_verdict_even_with_strong_effect(self):
        a = baseline()
        b = control(a)
        report = diagnostic_summary(a, b, n_boot=200, seed=11)
        self.assertEqual(report["status"], "DIAGNOSTIC_ONLY")
        self.assertEqual(report["mode"], "blank")
        self.assertLess(report["classes"]["C1"]["control_minus_reference"], 0)
        self.assertAlmostEqual(report["classes"]["C2"]["control_minus_reference"], 0)
        self.assertEqual(report["mixture"]["n_molecules"], 4)
        self.assertNotIn("verdict", report["mixture"])
        self.assertEqual(len(report["mixture"]["scenarios"]), 12)
        with self.assertRaisesRegex(ValueError, "diagnostic only"):
            align_runs(a, b)

    def test_tampered_original_spectra_or_exclusions_fail_closed(self):
        a = baseline()
        b = control(a)
        b["spectral_control"]["source_query_signatures"][2] = "f"*64
        with self.assertRaisesRegex(ValueError, "source spectra"):
            align_control_runs(a, b)

    def test_changed_query_population_or_model_fails_closed(self):
        a = baseline()
        for change, error in [
            (lambda x: x["query_ids"].__setitem__(0, "NEW"), "population"),
            (lambda x: x["config"].update({"top_n": 12}), "configurations"),
            (lambda x: x["query_sampling"].update({"max_spectra": 1}), "sampling"),
            (lambda x: x["query_sampling"].update({"control_seed": 99}), "contradicts"),
            (lambda x: x.update({"query_signature_version": "v1"}), "recipe"),
        ]:
            b = control(a)
            change(b)
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValueError, error):
                    align_control_runs(a, b)

    def test_reject_missing_or_forged_manifest(self):
        a = baseline()
        for modify in [
            lambda b: b["spectral_control"].pop("source_query_signatures"),
            lambda b: b["spectral_control"].update(source_query_signatures=["x"]),
            lambda b: b.update(query_signatures=["g"*64]*4),
            lambda b: b["spectral_control"].update(seed=-1),
            lambda b: b.update(spectral_control={"mode": "none", "seed": 0}),
        ]:
            b = control(a)
            modify(b)
            with self.assertRaises(ValueError):
                align_control_runs(a, b)

    def test_control_without_observable_peak_perturbation_rejected(self):
        a = baseline()
        b = control(a)
        b["query_signatures"] = b["spectral_control"]["source_query_signatures"].copy()
        with self.assertRaisesRegex(ValueError, "No actual spectral perturbation"):
            align_control_runs(a, b)

    def test_diagnostic_bootstrap_validation(self):
        a = baseline()
        b = control(a, mode="intensity_shuffle")
        self.assertEqual(
            diagnostic_summary(a, b, n_boot=200)["mode"], "intensity_shuffle"
        )
        with self.assertRaisesRegex(ValueError, "bootstrap draws"):
            diagnostic_summary(a, b, n_boot=99)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            diagnostic_summary(a, b, seed=-1)


if __name__ == "__main__":
    unittest.main()
