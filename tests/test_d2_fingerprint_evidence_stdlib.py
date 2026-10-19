"""D2 source evidence acceptance gates; entirely synthetic standard-library tests."""
import unittest
from copy import deepcopy

from casmi.validate.run_alignment import align_runs
from casmi.validate.source_manifest import MANIFEST_SCHEMA, required_validation_inputs


def run(channels=("d1", "d2"), *, fingerprint=True):
    report = {
        "tag": "synthetic",
        "n_queries": 2,
        "query_ids": ["a", "b"],
        "query_signatures": ["a" * 64, "b" * 64],
        "query_signature_version": "casmi-query-spectra-v2",
        "classes": {cls: {"ranks": [1, 0]} for cls in ("C1", "C2", "C3")},
        "config": {"channels": list(channels)},
        "source_manifest_sha256": "1" * 64,
        "source_manifest_schema": MANIFEST_SCHEMA,
        "source_manifest_pool_tag": "",
        "source_manifest_coverage": "required-validation-inputs-v1",
        "source_manifest_required_inputs": list(
            required_validation_inputs(analog="d2" in channels)
        ),
        "source_manifest_file_count": 7,
        "fingerprint_provenance": None,
        "analog_index_provenance": None,
    }
    if "d2" in channels:
        report["config"]["analog"] = {"reference_policy": "legacy"}
        report["analog_index_provenance"] = {
            "schema": "casmi-d2-fresh-index-v1",
            "status": "REBUILT_WITHOUT_CACHE",
            "index_sha256": "4" * 64,
            "n_references": 50,
            "n_nonzero": 1000,
            "reference_policy": "legacy",
            "cache_used": False,
        }
    if "d2" in channels and fingerprint:
        report["fingerprint_provenance"] = {
            "status": "POOL_FINGERPRINTS_VERIFIED",
            "schema": "casmi-morgan-fingerprint-provenance-v1",
            "n_structures": 2,
            "pool_rows_sha256": "2" * 64,
            "fingerprint_bytes_sha256": "3" * 64,
        }
    return report


class D2EvidenceTests(unittest.TestCase):
    def test_d2_requires_fp_evidence_for_accept_reject(self):
        a, b = run(), run()
        self.assertTrue(align_runs(a, b)["decision_eligible"])
        b["fingerprint_provenance"] = None
        self.assertTrue(align_runs(a, b)["verified"])
        self.assertFalse(align_runs(a, b)["decision_eligible"])
        a["fingerprint_provenance"] = None
        self.assertFalse(align_runs(a, b)["source_verified"])

    def test_tampered_vectors_and_malformed_claim_rejected(self):
        a, b = run(), run()
        for mutation, expected in [
            (lambda x: x["fingerprint_provenance"].update(
                fingerprint_bytes_sha256="f" * 64), "fingerprint source mismatch"),
            (lambda x: x["fingerprint_provenance"].update(
                pool_rows_sha256="f" * 64), "fingerprint source mismatch"),
            (lambda x: x["fingerprint_provenance"].update(
                n_structures=9), "fingerprint source mismatch"),
            (lambda x: x["fingerprint_provenance"].update(
                schema="wrong"), "schema or status"),
            (lambda x: x["fingerprint_provenance"].update(
                status="not-verified"), "schema or status"),
            (lambda x: x["fingerprint_provenance"].update(
                fingerprint_bytes_sha256="bad"), "fingerprint provenance"),
            (lambda x: x["fingerprint_provenance"].update(
                n_structures=True), "row count"),
        ]:
            modified = deepcopy(b)
            mutation(modified)
            with self.subTest(error=expected):
                with self.assertRaisesRegex(ValueError, expected):
                    align_runs(a, modified)

    def test_d1_d2_ablation_can_share_full_manifest(self):
        a, b = run(("d1",), fingerprint=False), run()
        self.assertTrue(align_runs(a, b)["decision_eligible"])
        self.assertTrue(align_runs(b, a)["decision_eligible"])

    def test_d1_comparison_does_not_need_fp_sidecar(self):
        a, b = run(("d1",), fingerprint=False), run(("d1",), fingerprint=False)
        self.assertTrue(align_runs(a, b)["decision_eligible"])

    def test_cached_legacy_index_is_exploratory_even_with_verified_fingerprints(self):
        a, b = run(), run()
        a["analog_index_provenance"] = None
        self.assertTrue(align_runs(a, b)["verified"])
        self.assertFalse(align_runs(a, b)["decision_eligible"])
        a = run()
        b["analog_index_provenance"] = None
        self.assertFalse(align_runs(a, b)["decision_eligible"])

    def test_cache_claims_and_changed_reference_policy_fail_closed(self):
        a, b = run(), run()
        for mutation, expected in (
            (lambda p: p.update(cache_used=True), "bypass"),
            (lambda p: p.update(index_sha256="bad"), "content digest"),
            (lambda p: p.update(n_references=0), "reference count"),
            (lambda p: p.update(n_nonzero=-1), "nonzero count"),
            (lambda p: p.update(reference_policy="source_balanced"), "contradicts"),
            (lambda p: p.update(status="LOADED_FROM_CACHE"), "schema or status"),
        ):
            changed = deepcopy(b)
            mutation(changed["analog_index_provenance"])
            with self.subTest(expected=expected):
                with self.assertRaisesRegex(ValueError, expected):
                    align_runs(a, changed)

    def test_separately_fresh_rebuilt_d2_policy_variants_are_comparable(self):
        a, b = run(), run()
        b["config"]["analog"]["reference_policy"] = "source_balanced"
        b["analog_index_provenance"]["reference_policy"] = "source_balanced"
        b["analog_index_provenance"]["index_sha256"] = "5" * 64
        self.assertTrue(align_runs(a, b)["decision_eligible"])

    def test_source_snapshot_is_still_required(self):
        a, b = run(), run()
        b["source_manifest_sha256"] = "9" * 64
        with self.assertRaisesRegex(ValueError, "Source artifact snapshot mismatch"):
            align_runs(a, b)
        b = run()
        b["source_manifest_sha256"] = None
        with self.assertRaisesRegex(ValueError, "lacks a manifest"):
            align_runs(a, b)


if __name__ == "__main__":
    unittest.main()
