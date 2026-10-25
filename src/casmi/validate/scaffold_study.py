"""Predeclared internal chemical-family partitions and joint paired uncertainty."""
from __future__ import annotations

import hashlib
import json

import numpy as np


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def family_key(smiles):
    """Generic Murcko ring topology; generic whole connectivity for acyclics.

    Atomic/bond types are intentionally removed to make exclusion conservative.
    This is a computational partition, not a chemical similarity guarantee.
    """
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    try:
        mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
        if mol is None or not mol.GetNumHeavyAtoms():
            return None
        parts = Chem.GetMolFrags(mol, asMols=True)
        mol = max(parts, key=lambda part: (part.GetNumHeavyAtoms(), Chem.MolToSmiles(part)))
        core = MurckoScaffold.GetScaffoldForMol(mol)
        prefix = "ring:" if core.GetNumHeavyAtoms() else "acyclic:"
        core = core if core.GetNumHeavyAtoms() else mol
        return prefix + Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(core), canonical=True)
    except (ValueError, RuntimeError, TypeError):
        return None


def partition(records, sizes, salt):
    """One representative per family; union metric identities across families.

    Input: key/family/metric plus query metadata. No model scores are accepted.
    The union prevents metric-equivalent labels crossing an unusual core change.
    """
    parents = {}
    def find(value):
        parents.setdefault(value, value)
        while parents[value] != value:
            parents[value] = parents[parents[value]]
            value = parents[value]
        return value
    first = {}
    for row in records:
        if not row.get("family") or not row.get("metric"):
            raise ValueError("Every partition row requires a valid family and metric identity")
        family = find(row["family"])
        previous = first.setdefault(row["metric"], family)
        a, b = find(previous), find(family)
        if a != b:
            parents[max(a, b)] = min(a, b)
    groups = {}
    for row in records:
        groups.setdefault(find(row["family"]), []).append(row)
    ordered = sorted(groups, key=lambda family: hashlib.sha256((salt + "|" + family).encode()).hexdigest())
    if set(sizes) != {"train", "dev", "test"} or any(type(n) is not int or n < 1 for n in sizes.values()):
        raise ValueError("Positive train/dev/test counts required")
    if len(ordered) < sum(sizes.values()):
        raise ValueError("Insufficient independent family groups for declared population")
    selected, offset = {}, 0
    for split in ("train", "dev", "test"):
        rows = []
        for family in ordered[offset:offset + sizes[split]]:
            representative = min(groups[family], key=lambda row: row["key"])
            rows.append({**representative, "group": family,
                         "fold": int(hashlib.sha256((salt + "|fold|" + family).encode()).hexdigest(), 16) % 2})
        selected[split] = rows
        offset += sizes[split]
    for identity in ("family", "group", "metric", "key"):
        sets = {name: {row[identity] for row in values} for name, values in selected.items()}
        if any(sets[a] & sets[b] for a, b in (("train", "dev"), ("train", "test"), ("dev", "test"))):
            raise ValueError("Partition identity overlap")
    return selected


def joint_intervals(effects, *, n_boot=4000, seed=0):
    """Joint max-absolute-deviation paired bootstrap across a fixed outcome family."""
    names = list(effects)
    values = np.column_stack([effects[name] for name in names])
    if values.ndim != 2 or not len(values) or not np.isfinite(values).all() or n_boot < 100:
        raise ValueError("Finite paired effect vectors and at least 100 draws required")
    mean = values.mean(axis=0)
    rng = np.random.default_rng(seed)
    maximum = np.empty(n_boot)
    for start in range(0, n_boot, 100):
        ids = rng.integers(len(values), size=(min(100, n_boot - start), len(values)))
        maximum[start:start + len(ids)] = np.max(np.abs(values[ids].mean(axis=1) - mean), axis=1)
    radius = float(np.quantile(maximum, .95))
    return {"method": "molecule-paired max-deviation bootstrap over all prespecified secondary contrasts",
            "n_molecules": len(values), "n_bootstrap": n_boot, "seed": seed,
            "family_size": len(names), "joint_radius": radius,
            "outcomes": {name: {"delta_mrr": float(value), "ci_low": float(value - radius),
                                "ci_high": float(value + radius)} for name, value in zip(names, mean)},
            "limitations": "Approximate fixed-family uncertainty; no correction for domain shift or historical hypothesis selection"}
