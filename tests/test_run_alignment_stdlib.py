"""Guard against unpaired molecule identities, spectrum changes and legacy CI misuse."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.run_alignment import align_runs, query_signature  # noqa: E402


def run(ids, signatures=None):
    example = {
        "tag": "baseline", "n_queries": len(ids), "query_ids": list(ids),
        "query_signatures": signatures or [str(i) * 64 for i in range(1, len(ids) + 1)],
        "classes": {cls: {"ranks": [1] + [0] * (len(ids) - 1)}
                    for cls in ("C1", "C2", "C3")},
    }
    return example


class RunAlignmentTests(unittest.TestCase):
    def test_permutation_realigns_candidate_by_query_id(self):
        a = run(["a", "b"])
        a["classes"]["C1"]["ranks"] = [1, 0]
        b = run(["b", "a"], ["2" * 64, "1" * 64])
        b["classes"]["C1"]["ranks"] = [2, 1]
        result = align_runs(a, b)
        self.assertTrue(result["verified"])
        self.assertEqual(result["a"]["C1"], [1, 0])
        self.assertEqual(result["b"]["C1"], [1, 2])

    def test_same_count_different_query_ids_rejected(self):
        a, b = run(["a", "b"]), run(["a", "c"])
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            align_runs(a, b)

    def test_same_ids_changed_spectra_rejected(self):
        a, b = run(["a", "b"]), run(["a", "b"])
        b["query_signatures"][1] = "f" * 64
        with self.assertRaisesRegex(ValueError, "Spectrum/knockout mismatch"):
            align_runs(a, b)

    def test_signature_algorithm_versions_must_agree(self):
        a, b = run(["a"]), run(["a"])
        a["query_signature_version"] = "casmi-query-spectra-v1"
        b["query_signature_version"] = "casmi-query-spectra-v2"
        with self.assertRaisesRegex(ValueError, "signature algorithm versions"):
            align_runs(a, b)

    def test_duplicate_query_ids_rejected(self):
        a = run(["a", "a"])
        with self.assertRaisesRegex(ValueError, "duplicate query IDs"):
            align_runs(a, run(["a", "b"]))

    def test_missing_or_corrupt_manifest_rejected(self):
        a, b = run(["a"]), run(["a"])
        del b["query_signatures"]
        with self.assertRaisesRegex(ValueError, "legacy runs"):
            align_runs(a, b)
        b["query_signatures"] = ["not-a-sha256"]
        with self.assertRaisesRegex(ValueError, "malformed"):
            align_runs(a, b)

    def test_legacy_requires_explicit_override_and_cannot_verify(self):
        a, b = run(["a"]), run(["a"])
        for item in (a, b):
            item.pop("query_ids")
            item.pop("query_signatures")
        with self.assertRaisesRegex(ValueError, "legacy runs"):
            align_runs(a, b)
        self.assertFalse(
            align_runs(a, b, allow_unverified_legacy_order=True)["verified"]
        )

    def test_unequal_population_and_rank_shape_rejected(self):
        with self.assertRaisesRegex(ValueError, "same positive number"):
            align_runs(run(["a"]), run(["a", "b"]))
        a, b = run(["a", "b"]), run(["a", "b"])
        b["classes"]["C2"]["ranks"] = [1]
        with self.assertRaisesRegex(ValueError, "per-query ranks"):
            align_runs(a, b)

    def test_reject_bad_rank_and_inconsistent_summary(self):
        a, b = run(["a"]), run(["a"])
        b["classes"]["C1"]["ranks"] = [26]
        with self.assertRaisesRegex(ValueError, "0 to 25"):
            align_runs(a, b)
        b["classes"]["C1"]["ranks"] = [1]
        b["classes"]["C1"]["mrr"] = 0.2
        with self.assertRaisesRegex(ValueError, "disagrees"):
            align_runs(a, b)

    def test_query_spectral_fingerprint_is_deterministic_and_sensitive(self):
        q = SimpleNamespace(
            qid="enveda:test1", pool_key="mol", spectrum_rows=[3, 5],
            spectra=[dict(mz=[125.3, 230.2], it=[0.5, 1.0],
                          precursor_mz=270.1, adduct="[M+H]+",
                          ce_ev=float("nan"), mode="positive")],
        )
        fingerprint = query_signature(q)
        self.assertEqual(fingerprint, query_signature(q))
        self.assertEqual(len(fingerprint), 64)
        changed = deepcopy(q)
        changed.spectra[0]["it"][1] = 0.8
        self.assertNotEqual(fingerprint, query_signature(changed))
        changed = deepcopy(q)
        changed.spectrum_rows = [3, 6]
        self.assertNotEqual(fingerprint, query_signature(changed))


if __name__ == "__main__":
    unittest.main()
