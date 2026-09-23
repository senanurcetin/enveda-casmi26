"""Lazy, persistent map from pool key (standard InChIKey14) to metric key (tautomer-canonical).

Measured: the two keys agree for only 88.5% of COCONUT structures (natural products carry many
tautomerisable groups), so scoring on the standard key would be wrong ~1 time in 9. Computing the
metric key costs 10–300 ms per structure, so we compute it only where it matters — query truths,
the head of each ranked list, and output deduplication — and cache it on disk across runs.
"""
from pathlib import Path

import pandas as pd

from casmi.chem.canon import canonical_ik14


class MetricKeys:
    def __init__(self, cache_path: Path):
        self.path = Path(cache_path)
        self.map = {}
        if self.path.exists():
            t = pd.read_parquet(self.path)
            self.map = dict(zip(t.pool_key, t.metric_key))
        self._dirty = 0

    def get(self, pool_key: str, smiles: str):
        key = self.map.get(pool_key)
        if key is None and pool_key not in self.map:
            key = canonical_ik14(smiles)
            self.map[pool_key] = key
            self._dirty += 1
        return key

    def save(self):
        if self._dirty:
            pd.DataFrame({"pool_key": list(self.map), "metric_key": list(self.map.values())}).to_parquet(self.path, index=False)
            self._dirty = 0

    def dedupe_ranked(self, pool_keys, smiles_of, top_n=25):
        """Map a ranked list of pool keys to metric keys, drop unparseable and duplicate structures,
        keep the first `top_n`. Returns [(pool_key, metric_key)] in rank order."""
        out, seen = [], set()
        for pk in pool_keys:
            mk = self.get(pk, smiles_of(pk))
            if mk is None or mk in seen:
                continue
            seen.add(mk)
            out.append((pk, mk))
            if len(out) == top_n:
                break
        return out
