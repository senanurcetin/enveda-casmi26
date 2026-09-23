"""Structure-level knockout validation.

A query is one structure with all its spectra from the query library (e.g. enveda-np-examples),
mimicking a test molecule. Each query is scored under three simulated novelty classes:

  C1  library minus the query's own spectra        pool keeps the structure
  C2  library minus EVERY spectrum of the structure, in every library   pool keeps it
  C3  as C2, and the structure is removed from the pool                 (retrieval must score 0)

Knockout is by metric key (tautomer-canonical InChIKey14), never by library name — holding a
library out by name leaks, because the same structures appear in other libraries.
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

CLASSES = ("C1", "C2", "C3")


@dataclass
class Query:
    qid: str
    pool_key: str                 # the structure's pool key (standard InChIKey14)
    spectra: list                 # list of dicts: mz, it, precursor_mz, adduct, ce_ev, mode
    spectrum_rows: np.ndarray     # store rows of the query's own spectra (always excluded)


@dataclass
class Knockout:
    """What the engine must not see for one (query, class)."""
    excluded_rows: np.ndarray                 # store rows hidden from library search
    excluded_keys: set = field(default_factory=set)   # structures hidden from the pool


def build_queries(store, key_of_row, library="enveda-np-examples", max_queries=None, seed=0):
    """Group the query library's spectra by structure into molecule-like queries."""
    meta = store.meta
    rows = np.flatnonzero((meta.ingest_lib == library).to_numpy() & meta.test_adduct.to_numpy())
    keys = key_of_row[rows]
    queries = []
    for key in pd.unique(keys):
        if key is None or (isinstance(key, float) and np.isnan(key)):
            continue
        r = rows[keys == key]
        spectra = []
        for i in r:
            mz, it = store.peaks(i)
            m = meta.iloc[i]
            spectra.append(dict(mz=mz, it=it, precursor_mz=float(m.precursor_mz), adduct=str(m.adduct),
                                ce_ev=float(m.ce_ev), mode=str(m.ionization_mode)))
        queries.append(Query(qid=f"{library}:{key}", pool_key=key, spectra=spectra, spectrum_rows=r))
    if max_queries and len(queries) > max_queries:
        rng = np.random.default_rng(seed)
        queries = [queries[i] for i in sorted(rng.choice(len(queries), max_queries, replace=False))]
    return queries


def knockout_for(query: Query, cls: str, rows_of_key, equivalent_keys) -> Knockout:
    """`equivalent_keys`: every pool key whose metric key equals the truth's (tautomer variants
    stored under a different standard key). Knocking out only the query's own key would leak."""
    if cls == "C1":
        return Knockout(query.spectrum_rows)
    rows = [rows_of_key.get(k, np.empty(0, np.int64)) for k in equivalent_keys]
    excluded = np.union1d(np.concatenate(rows) if rows else np.empty(0, np.int64), query.spectrum_rows)
    return Knockout(excluded, set(equivalent_keys) if cls == "C3" else set())


def mixture(mrr_by_class: dict, f1: float, f3: float) -> float:
    """Expected leaderboard score for class shares f1, f2 = 1 - f1 - f3, f3."""
    f2 = 1.0 - f1 - f3
    return f1 * mrr_by_class["C1"] + f2 * mrr_by_class["C2"] + f3 * mrr_by_class["C3"]


def mixture_grid(mrr_by_class, f1_values=(0.15, 0.20, 0.25, 0.30), f3_values=(0.2, 0.3, 0.4)):
    rows = [(f1, f3, mixture(mrr_by_class, f1, f3)) for f1 in f1_values for f3 in f3_values if f1 + f3 < 1]
    return pd.DataFrame(rows, columns=["f1", "f3", "score"])
