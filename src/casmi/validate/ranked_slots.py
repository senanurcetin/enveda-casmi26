"""Position-preserving parsing for CASMI ranked SMILES predictions.

A semicolon defines the next ranked position even when the candidate is
blank or invalid. Never filter an empty field when computing MRR.
"""
from __future__ import annotations


def split_ranked_guesses(text: str) -> list[str]:
    """Keep every guess slot, including interior/leading/trailing blanks."""
    if not isinstance(text, str):
        raise TypeError("Ranked SMILES must be a string")
    return [part.strip() for part in text.split(";")]
