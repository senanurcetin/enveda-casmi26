"""Rank-slot parsing never silently optimizes an invalid submission."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.ranked_slots import split_ranked_guesses  # noqa: E402


class RankedSlotsTests(unittest.TestCase):
    def test_blank_interior_candidate_costs_its_position(self):
        self.assertEqual(split_ranked_guesses("C;;CCO"), ["C", "", "CCO"])
        self.assertEqual(split_ranked_guesses(" ;CCO"), ["", "CCO"])

    def test_blank_trailing_candidate_still_counts(self):
        self.assertEqual(split_ranked_guesses("C;"), ["C", ""])
        self.assertEqual(split_ranked_guesses(""), [""])

    def test_whitespace_is_cleaned_without_reordering(self):
        self.assertEqual(split_ranked_guesses(" C ; CCO ; O "),
                         ["C", "CCO", "O"])

    def test_non_string_rejected(self):
        with self.assertRaises(TypeError):
            split_ranked_guesses(None)


if __name__ == "__main__":
    unittest.main()
