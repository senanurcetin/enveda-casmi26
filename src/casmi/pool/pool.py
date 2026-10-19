"""A3 — candidate pool: unique structures (metric key) sorted by neutral monoisotopic mass.

Sources are tracked as a bitmask so the ranker can learn source priors (e.g. that structures known
only from the drug-like enveda-180 library are unlikely answers for a natural-product test set)
and so validation can knock out a structure from the pool (Class-3 simulation).
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

SOURCE_BITS = {"train": 1, "enveda180_only": 2, "coconut": 4, "chebi": 8, "lipidmaps": 16, "pubchem": 32, "pubchem_lit": 64}


@dataclass
class Window:
    idx: np.ndarray          # pool row indices
    ppm: np.ndarray          # signed ppm of each candidate relative to the query mass


class CandidatePool:
    def __init__(self, table: pd.DataFrame):
        """table: pool_key, smiles, mass, sources (bitmask) [+ any extra per-structure columns]."""
        needed = {"pool_key", "smiles", "mass"}
        if not needed.issubset(table.columns):
            raise ValueError(f"Candidate pool missing columns: {sorted(needed - set(table.columns))}")
        # Missing candidates can be dropped as before, but present impossible
        # masses must not silently produce NaN/inf ppm values or divide by 0.
        t = table.dropna(subset=["pool_key", "mass"]).copy()
        try:
            mass = pd.to_numeric(t["mass"], errors="raise").to_numpy(dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError("Candidate pool contains non-numeric masses") from exc
        if np.any(~np.isfinite(mass)) or np.any(mass <= 0.0):
            raise ValueError("Candidate pool masses must be finite and strictly positive")
        t["mass"] = mass
        t = t.sort_values("mass", kind="stable")
        t = t.drop_duplicates("pool_key", keep="first").reset_index(drop=True)
        self.table = t
        self.mass = t["mass"].to_numpy(np.float64)
        self.key = t["pool_key"].to_numpy()
        self.smiles = t["smiles"].to_numpy()
        self.row_of = dict(zip(self.key, range(len(t))))
        # generated (Class-3 hedge) structures live outside the table: key -> SMILES / packed fingerprint
        self.extra_smiles, self.extra_fp = {}, {}

    def smiles_of(self, key):
        r = self.row_of.get(key)
        return self.smiles[r] if r is not None else self.extra_smiles[key]

    def register_generated(self, keys, smiles):
        """Atomically register chemically valid generated structures.

        Generated candidate IDs are 'G' + standard InChIKey14, as produced
        by BaselineEngine._generated. Invalid records cannot contribute
        all-zero fingerprints or leave a half-written candidate overlay.
        """
        from rdkit import Chem
        from rdkit.Chem import rdFingerprintGenerator

        keys = list(keys)
        smiles = list(smiles)
        if len(keys) != len(smiles):
            raise ValueError("Generated key and SMILES counts differ")
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        staged = {}
        for i, (key, smi) in enumerate(zip(keys, smiles)):
            if key in self.extra_smiles or key in staged:
                continue
            if not isinstance(key, str) or not isinstance(smi, str):
                raise ValueError(f"Malformed generated candidate at row {i}")
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None or mol.GetNumHeavyAtoms() == 0:
                    raise ValueError
                expected = "G" + Chem.MolToInchiKey(mol)[:14]
                if key != expected:
                    raise ValueError
                bits = gen.GetFingerprintAsNumPy(mol).astype(np.uint8)
                if bits.size != 2048 or not np.any(bits):
                    raise ValueError
            except (ValueError, TypeError, RuntimeError):
                # Never log the source SMILES or derived molecular identity.
                raise ValueError(f"Invalid generated candidate identity at row {i}") from None
            staged[key] = (smi, np.packbits(bits))
        for key, (smi, fp) in staged.items():
            self.extra_smiles[key] = smi
            self.extra_fp[key] = fp

    def packed_fps(self, keys, fps):
        """Packed fingerprints for pool and generated keys alike."""
        out = np.empty((len(keys), fps.shape[1]), fps.dtype)
        for i, k in enumerate(keys):
            r = self.row_of.get(k)
            out[i] = fps[r] if r is not None else self.extra_fp[k]
        return out

    def __len__(self):
        return len(self.mass)

    def window(self, query_mass, tol_ppm=10.0, exclude=None) -> Window:
        """All candidates with |ppm(candidate vs query)| <= tol_ppm, minus excluded keys."""
        try:
            if isinstance(query_mass, (bool, np.bool_)) or not np.isscalar(query_mass):
                raise ValueError
            mass_value = float(query_mass)
            if not np.isfinite(mass_value) or mass_value <= 0:
                raise ValueError
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("Query mass must be a finite positive scalar") from exc
        try:
            if isinstance(tol_ppm, (bool, np.bool_)) or not np.isscalar(tol_ppm):
                raise ValueError
            tolerance = float(tol_ppm)
            if not np.isfinite(tolerance) or tolerance < 0:
                raise ValueError
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("Mass tolerance ppm must be a finite nonnegative scalar") from exc
        d = mass_value * tolerance * 1e-6
        lo = np.searchsorted(self.mass, mass_value - d, side="left")
        hi = np.searchsorted(self.mass, mass_value + d, side="right")
        idx = np.arange(lo, hi)
        if exclude is not None:
            excluded = list(exclude)
            if excluded:
                idx = idx[~np.isin(self.key[idx], excluded)]
        return Window(idx, (mass_value - self.mass[idx]) / self.mass[idx] * 1e6)
