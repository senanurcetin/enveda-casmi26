"""Monoisotopic masses, formula parsing and adduct arithmetic.

Everything here is plain Python/NumPy (no RDKit) so it can run over millions of rows.
The adduct parser handles every string in train (`[2M+Na]+`, `[M-2H2O+H]+`, `[M+CH2O2-H]-`, ...),
not only the ten test adducts, because library spectra with other adducts are still evidence.
"""
import re
from functools import lru_cache

import numpy as np

ELECTRON = 0.000548579909
PROTON = 1.007276466621

MONOISOTOPIC = {
    "H": 1.00782503207, "D": 2.0141017778, "C": 12.0, "N": 14.0030740048, "O": 15.99491461956,
    "P": 30.97376163, "S": 31.97207100, "F": 18.99840322, "Cl": 34.96885268, "Br": 78.9183371,
    "I": 126.904473, "Si": 27.9769265325, "B": 11.0093054, "Se": 79.9165213, "Na": 22.9897692809,
    "K": 38.96370668, "Li": 7.01600455, "Mg": 23.9850417, "Ca": 39.96259098, "Fe": 55.9349375,
    "Cu": 62.9295975, "Zn": 63.9291422, "As": 74.9215965, "Al": 26.98153863, "Co": 58.933195,
    "Ni": 57.9353429, "Mn": 54.9380451, "Cr": 51.9405075, "Sn": 119.9021947, "Hg": 201.970643,
    "Pt": 194.9647911, "Ag": 106.905097, "Au": 196.9665687, "Gd": 157.9241039, "Ge": 73.9211778,
    "Ti": 47.9479463, "V": 50.9439595, "Mo": 97.9054082, "W": 183.9509312, "Te": 129.9062244,
    "Ru": 101.9043493, "Rh": 102.905504, "Pd": 105.903486, "Sb": 120.9038157, "Bi": 208.9803987,
    "Ba": 137.9052472, "Sr": 87.9056121, "Tl": 204.9744275, "Pb": 207.9766521, "Cs": 132.905451933,
    "Rb": 84.911789738,
}

_ELEMENT_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


@lru_cache(maxsize=200_000)
def formula_counts(formula: str) -> dict:
    """'C15H10O5' -> {'C': 15, 'H': 10, 'O': 5}. Ignores charge suffixes such as '+' or '-'."""
    counts = {}
    for element, n in _ELEMENT_RE.findall(formula.rstrip("+-")):
        counts[element] = counts.get(element, 0) + (int(n) if n else 1)
    return counts


@lru_cache(maxsize=200_000)
def formula_mass(formula: str) -> float:
    """Monoisotopic neutral mass of a molecular formula; NaN if it contains an element we do not
    tabulate (organometallics such as Ir or Os, far outside this competition's chemistry)."""
    counts = formula_counts(formula)
    if any(e not in MONOISOTOPIC for e in counts):
        return float("nan")
    return sum(MONOISOTOPIC[e] * n for e, n in counts.items())


# ----------------------------------------------------------------------------- adducts
_ADDUCT_RE = re.compile(r"^\[(\d*)M((?:[+-]\d*[A-Za-z0-9]+)*)\](\d*)([+-])$")
_PART_RE = re.compile(r"([+-])(\d*)([A-Za-z][A-Za-z0-9]*)")
# Common adduct shorthands that are not literal formulas
_ALIASES = {"ACN": "C2H3N", "FA": "CH2O2", "Hac": "C2H4O2", "HAc": "C2H4O2", "TFA": "C2HF3O2", "IsoProp": "C3H8O"}

TEST_ADDUCTS = ("[M+H]+", "[M+NH4]+", "[M-H2O+H]+", "[M-2H2O+H]+", "[M+Na]+", "[M+K]+",
                "[M-H]-", "[M-H2O-H]-", "[M+CH2O2-H]-", "[M+Cl]-")


@lru_cache(maxsize=4096)
def parse_adduct(adduct: str):
    """Return (n_molecules, mass_shift, charge) with ion_mz = (n*M + mass_shift) / |charge|.

    `mass_shift` already includes the electron mass for the charge state. Returns None when the
    string is not a recognisable adduct (the caller decides whether to drop the spectrum).
    """
    if not isinstance(adduct, str):
        return None
    m = _ADDUCT_RE.match(adduct.strip())
    if not m:
        return None
    n_mol = int(m.group(1)) if m.group(1) else 1
    z = int(m.group(3)) if m.group(3) else 1
    charge = z if m.group(4) == "+" else -z
    shift = 0.0
    for sign, mult, part in _PART_RE.findall(m.group(2)):
        mass = formula_mass(_ALIASES.get(part, part))
        if np.isnan(mass):
            return None
        shift += (1 if sign == "+" else -1) * (int(mult) if mult else 1) * mass
    shift -= charge * ELECTRON          # a cation has lost electrons, an anion gained them
    return n_mol, shift, charge


def ion_mz(neutral_mass: float, adduct: str) -> float:
    """Ion m/z for a neutral mass and adduct; NaN when the adduct string is not recognised."""
    parsed = parse_adduct(adduct)
    if parsed is None:
        return np.nan
    n, shift, charge = parsed
    return (n * neutral_mass + shift) / abs(charge)


def neutral_mass(precursor_mz, adduct: str):
    """Invert `ion_mz`. Works on scalars or NumPy arrays of precursor m/z for one adduct."""
    parsed = parse_adduct(adduct)
    if parsed is None:
        return np.nan if np.ndim(precursor_mz) == 0 else np.full(np.shape(precursor_mz), np.nan)
    n, shift, charge = parsed
    return (np.asarray(precursor_mz, dtype=float) * abs(charge) - shift) / n


def ppm(measured, theoretical):
    return (np.asarray(measured, dtype=float) - theoretical) / theoretical * 1e6
