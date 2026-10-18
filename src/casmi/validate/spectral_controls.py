"""Reproducible negative-control perturbations of query MS/MS peaks.

These transformations are *diagnostics*, never submission-model features:
- blank: suppress all fragment evidence while preserving precursor mass,
  adduct, ionization mode and collision energy; this exposes non-spectral
  candidate-pool/mass priors.
- intensity_shuffle: independently permute fragment intensities within each
  spectrum. m/z and precursor stay fixed, while intensity/mz relationships
  are intentionally broken. Identical intensities cannot be changed.

Original query spectra are never mutated. For molecule-level reproducibility
permutations derive from (seed, query ID, spectrum index), not list position.
"""
from __future__ import annotations

import hashlib

import numpy as np

MODES = ("none", "blank", "intensity_shuffle")


def perturb_spectra(spectra, *, mode: str = "none",
                    seed: int = 0, query_id: str = "") -> list[dict]:
    """Make an independent copy of each spectrum, perturbing only peak arrays."""
    if mode not in MODES:
        raise ValueError(f"Unsupported spectral control: {mode}")
    if type(seed) is not int or seed < 0:
        raise ValueError("Control seed must be a nonnegative integer")
    if not isinstance(query_id, str) or not query_id:
        raise ValueError("query_id must identify one molecule")
    controlled = []
    for index, spectrum in enumerate(spectra):
        original_mz = np.array(spectrum["mz"], copy=True)
        original_it = np.array(spectrum["it"], copy=True)
        if original_mz.ndim != 1 or original_it.ndim != 1 or len(original_mz) != len(original_it):
            raise ValueError("Fragment m/z and intensity arrays must have matching 1D shapes")
        item = dict(spectrum)
        if mode == "blank":
            item["mz"] = np.empty(0, dtype=original_mz.dtype)
            item["it"] = np.empty(0, dtype=original_it.dtype)
        elif mode == "intensity_shuffle" and len(original_it) > 1:
            digest = hashlib.sha256(
                f"casmi-negative-control-v1:{seed}:{query_id}:{index}".encode()
            ).digest()
            rng = np.random.default_rng(int.from_bytes(digest[:8], "big"))
            permutation = rng.permutation(len(original_it))
            if np.array_equal(permutation, np.arange(len(original_it))):
                permutation = np.roll(permutation, 1)
            item["mz"] = original_mz
            item["it"] = original_it[permutation]
        else:
            item["mz"], item["it"] = original_mz, original_it
        controlled.append(item)
    return controlled
