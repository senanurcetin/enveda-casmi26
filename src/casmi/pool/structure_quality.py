"""Explicit whole-candidate quarantine before immutable pool publication.

This verifies parseability and nonzero Morgan features, not chemical validity,
provider identity accuracy or biological plausibility. No row details are logged.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np


def _invalid_reasons(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    reasons = []
    for value in smiles:
        try:
            mol = Chem.MolFromSmiles(value) if isinstance(value, str) else None
            if mol is None:
                reasons.append("unparseable")
            elif mol.GetNumHeavyAtoms() == 0:
                reasons.append("no_heavy_atoms")
            elif generator.GetFingerprint(mol).GetNumOnBits() == 0:
                reasons.append("zero_morgan")
            else:
                reasons.append("")
        except (ValueError, TypeError, RuntimeError):
            reasons.append("fingerprint_error")
    return reasons


def quarantine_invalid_candidates(pool, *, workers=1):
    """Return retained rows and aggregate accounting; never silently zero-fill."""
    if not pool.pool_key.is_unique:
        raise ValueError("Structure quarantine requires one candidate per identity")
    values = pool.smiles.tolist()
    chunks = [values[start:start + 5000] for start in range(0, len(values), 5000)]
    if workers == 1:
        nested = [_invalid_reasons(chunk) for chunk in chunks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            nested = list(executor.map(_invalid_reasons, chunks))
    reasons = np.asarray([reason for chunk in nested for reason in chunk], dtype=object)
    kept = pool.loc[reasons == ""].copy().reset_index(drop=True)
    if kept.empty:
        raise ValueError("No candidates remain after structure quarantine")
    return kept, {
        "policy": "quarantine_invalid_selected_candidate_identity",
        "input_candidates": len(pool), "retained_candidates": len(kept),
        "quarantined_candidates": int(np.count_nonzero(reasons != "")),
        "reasons": {name: int(np.count_nonzero(reasons == name))
                    for name in ("unparseable", "no_heavy_atoms", "zero_morgan", "fingerprint_error")},
        "scope": "Parseable, heavy-atom-containing, nonzero Morgan structures; not provider identity authentication",
    }
