"""Algebraic coverage ceilings must preserve the complete test denominator."""
import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path

from casmi.validate.coverage_budget import coverage_budget
from scripts.build_coverage_report import derive

ROOT = Path(__file__).resolve().parents[1]


class CoverageBudgetTests(unittest.TestCase):
    def test_known_rank_list_decomposes_without_dropping_misses(self):
        v = coverage_budget({"n": 4, "mrr": .385, "top1": .25, "top10": .5, "recall_at_25": .75}, 3)
        self.assertEqual(v["truth_candidates_absent"], 1)
        self.assertEqual(v["top25_hits"], 3)
        self.assertEqual(v["covered_outside_top25"], 0)
        self.assertAlmostEqual(v["conditional_mrr_given_coverage"], 1.54 / 3)
        self.assertAlmostEqual(v["coverage_loss"] + v["within_coverage_ranking_loss"] + v["observed_mrr"], 1)

    def test_zero_and_perfect_coverage_have_explicit_boundaries(self):
        zero = coverage_budget({"n": 4, "mrr": 0, "top1": 0, "top10": 0, "recall_at_25": 0}, 0)
        self.assertIsNone(zero["conditional_mrr_given_coverage"])
        self.assertIsNone(zero["conditional_hit25_given_coverage"])
        perfect = coverage_budget({"n": 4, "mrr": 1, "top1": 1, "top10": 1, "recall_at_25": 1}, 4)
        self.assertEqual(perfect["coverage_loss"] + perfect["within_coverage_ranking_loss"], 0)

    def test_impossible_molecule_counts_and_mrr_are_rejected(self):
        base = {"n": 4, "mrr": .385, "top1": .25, "top10": .5, "recall_at_25": .75}
        for patch, covered in [({"n": 0}, 0), ({}, 2), ({}, True), ({"top1": .1}, 3),
                               ({"mrr": .9}, 3), ({"mrr": .01}, 3), ({"top10": 0}, 3),
                               ({"mrr": float("nan")}, 3), ({"mrr": True}, 3)]:
            with self.subTest(patch=patch, covered=covered), self.assertRaises(ValueError):
                coverage_budget({**base, **patch}, covered)

    def test_committed_decomposition_is_bound_and_contains_all_variants(self):
        raw = (ROOT / "reports/advanced_study.json").read_text(encoding="utf-8").encode()
        source = json.loads(raw)
        publication = json.loads((ROOT / "reports/advanced_publication.json").read_text(encoding="utf-8"))
        derived = derive(source, publication, hashlib.sha256(raw).hexdigest())
        committed = json.loads((ROOT / "reports/candidate_loss_budget.json").read_text(encoding="utf-8"))
        self.assertEqual(derived, committed)
        self.assertEqual(set(derived["models"]), set(source["models"]))
        c3 = derived["models"]["full"]["C3"]
        self.assertEqual((c3["truth_candidates_absent"], c3["top25_hits"], c3["top1_hits"], c3["covered_outside_top25"]), (224, 24, 6, 2))
        for classes in derived["models"].values():
            for values in classes.values():
                self.assertAlmostEqual(sum(values[k] for k in ("observed_mrr", "within_coverage_ranking_loss", "coverage_loss")), 1)
        changed = deepcopy(publication)
        changed["aggregate_artifact_sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            derive(source, changed, hashlib.sha256(raw).hexdigest())


if __name__ == "__main__":
    unittest.main()
