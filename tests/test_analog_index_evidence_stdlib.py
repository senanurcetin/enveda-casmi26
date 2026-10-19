"""Standalone D2 index provenance contract tests; stdlib-only, synthetic."""
import unittest

from casmi.validate.analog_index_evidence import verified_analog_index_evidence


def run():
    return {
        "config": {"channels": ["d1", "d2"],
                   "analog": {"reference_policy": "legacy"}},
        "analog_index_provenance": {
            "schema": "casmi-d2-fresh-index-v1",
            "status": "REBUILT_WITHOUT_CACHE",
            "index_sha256": "b" * 64,
            "n_references": 5,
            "n_nonzero": 78,
            "reference_policy": "legacy",
            "cache_used": False,
        },
    }


class D2FreshIndexEvidenceTests(unittest.TestCase):
    def test_fresh_index_is_accepted_and_missing_is_exploratory(self):
        good = run()
        self.assertTrue(verified_analog_index_evidence(good))
        good.pop("analog_index_provenance")
        self.assertFalse(verified_analog_index_evidence(good))

    def test_cached_or_unverifiable_claims_are_rejected(self):
        mutations = [
            ("cache_used", True, "bypass"),
            ("schema", "unknown", "schema"),
            ("status", "LOADED_FROM_CACHE", "status"),
            ("index_sha256", "bad", "digest"),
            ("n_references", 0, "reference count"),
            ("n_references", True, "reference count"),
            ("n_nonzero", -1, "nonzero count"),
            ("reference_policy", "source_balanced", "contradicts"),
        ]
        for field, value, pattern in mutations:
            with self.subTest(field=field, value=value):
                bad = run()
                bad["analog_index_provenance"][field] = value
                with self.assertRaisesRegex(ValueError, pattern):
                    verified_analog_index_evidence(bad)

    def test_model_policy_is_not_inferred_from_proof(self):
        bad = run()
        bad["config"]["analog"].pop("reference_policy")
        with self.assertRaisesRegex(ValueError, "reference-selection policy"):
            verified_analog_index_evidence(bad)
        bad = run()
        bad["analog_index_provenance"] = "legacy"
        with self.assertRaisesRegex(ValueError, "Malformed"):
            verified_analog_index_evidence(bad)


if __name__ == "__main__":
    unittest.main()
