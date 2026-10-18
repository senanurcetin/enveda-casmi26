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
    spectrum_rows: np.ndarray     # selected source spectra used as query input
    source_library_rows: np.ndarray | None = None  # ALL same-key source-library spectra for C1 exclusion


@dataclass
class Knockout:
    """What the engine must not see for one (query, class)."""
    excluded_rows: np.ndarray                 # store rows hidden from library search
    excluded_keys: set = field(default_factory=set)   # structures hidden from the pool


def build_queries(store, key_of_row, library="enveda-np-examples", keys_subset=None,
                  max_spectra=None, max_queries=None, seed=0):
    """Group one library's test-adduct spectra by structure into molecule-like queries.
    `keys_subset` restricts the structures; `max_spectra` subsamples each query (test median is 3).
    A subsampled query still excludes ALL of that structure's reference rows in its own
    source library. Otherwise an unselected duplicate spectrum leaks into C1 scoring."""
    meta = store.meta
    rng = np.random.default_rng(seed)
    rows = np.flatnonzero((meta.ingest_lib == library).to_numpy() & meta.test_adduct.to_numpy()
                          & (meta.n_peaks.to_numpy() >= 3))
    keys = key_of_row[rows]
    if keys_subset is not None:
        m = np.isin(keys, list(keys_subset))
        rows, keys = rows[m], keys[m]
    order = np.argsort(keys.astype(str), kind="stable")
    rows, keys = rows[order], keys[order]
    starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1]) + 1] if len(keys) else np.empty(0, np.int64)
    ends = np.r_[starts[1:], len(keys)] if len(keys) else np.empty(0, np.int64)
    groups = list(zip(starts, ends))
    if max_queries and len(groups) > max_queries:
        groups = [groups[i] for i in sorted(rng.choice(len(groups), max_queries, replace=False))]
    # Collect all source-library rows for the selected query keys, not only
    # test-adduct/peak-filtered inputs. Reference indices have their own
    # filtering rules and cannot safely be inferred from sampled inputs.
    wanted = {str(keys[s_]) for s_, _ in groups}
    source_rows = np.flatnonzero((meta.ingest_lib == library).to_numpy())
    reference_rows = {key: [] for key in wanted}
    for index in source_rows:
        key = str(key_of_row[index])
        if key in reference_rows:
            reference_rows[key].append(int(index))

    queries = []
    for s_, e_ in groups:
        r = rows[s_:e_]
        if max_spectra and len(r) > max_spectra:
            r = np.sort(rng.choice(r, max_spectra, replace=False))
        spectra = []
        for i in r:
            mz, it = store.peaks(i)
            m = meta.iloc[i]
            spectra.append(dict(mz=mz, it=it, precursor_mz=float(m.precursor_mz), adduct=str(m.adduct),
                                ce_ev=float(m.ce_ev), mode=str(m.ionization_mode)))
        own_rows = np.asarray(reference_rows[str(keys[s_])], dtype=np.int64)
        queries.append(Query(
            qid=f"{library}:{keys[s_]}", pool_key=keys[s_],
            spectra=spectra, spectrum_rows=r,
            source_library_rows=own_rows,
        ))
    return queries


def knockout_for(query: Query, cls: str, rows_of_key, equivalent_keys) -> Knockout:
    """`equivalent_keys`: every pool key whose metric key equals the truth's (tautomer variants
    stored under a different standard key). Knocking out only the query's own key would leak."""
    if cls not in CLASSES:
        raise ValueError(f"Unknown knockout class: {cls}")
    query_rows = np.asarray(query.spectrum_rows, dtype=np.int64)
    own_rows = (query_rows if query.source_library_rows is None
                else np.union1d(query_rows, query.source_library_rows))
    if cls == "C1":
        return Knockout(own_rows)
    rows = [rows_of_key.get(k, np.empty(0, np.int64)) for k in equivalent_keys]
    excluded = np.union1d(
        np.concatenate(rows) if rows else np.empty(0, np.int64),
        own_rows,
    )
    return Knockout(excluded, set(equivalent_keys) if cls == "C3" else set())



def include_equivalent_source_rows(query: Query, equivalent_keys, source_rows, key_of_row):
    """Extend the query's C1 exclusions with its same-metric-key source refs.

    `source_rows` must contain only indices from the query's own library.
    `key_of_row` maps every store index to its *standard* structure key.
    Every equivalent-key variant must be removed, not just the sampled peaks.
    This mutates only the optional reference exclusion field, not spectra.
    """
    source_rows = np.asarray(source_rows, dtype=np.int64)
    key_of_row = np.asarray(key_of_row)
    if np.any(source_rows < 0) or np.any(source_rows >= len(key_of_row)):
        raise ValueError("Source-library reference row outside store")
    if not equivalent_keys:
        raise ValueError("Expected at least the query's own structure key")
    matching = source_rows[np.isin(key_of_row[source_rows], list(equivalent_keys))]
    base = query.spectrum_rows if query.source_library_rows is None else query.source_library_rows
    query.source_library_rows = np.union1d(base, matching).astype(np.int64, copy=False)
    return query.source_library_rows


def equivalent_keys(query_key, query_mass, pool, keys, window_ppm=30.0, truth_smiles=None):
    """Metric truth and every pool key in the mass window that is the same structure under the
    metric (tautomer variants stored under another standard key). The skeleton key is an exact
    necessary condition, so only skeleton matches pay for tautomer canonicalisation.
    `truth_smiles` is required when the query structure itself is not in this pool."""
    from casmi.chem.canon import skeleton_key
    if truth_smiles is None:
        truth_smiles = pool.smiles_of(query_key)
    truth, skel = keys.get(query_key, truth_smiles), skeleton_key(truth_smiles)
    same = {query_key} if query_key in pool.row_of else set()
    for k in pool.key[pool.window(query_mass, window_ppm).idx]:
        if k != query_key and skeleton_key(pool.smiles_of(k)) == skel and keys.get(k, pool.smiles_of(k)) == truth:
            same.add(k)
    return truth, same


def mixture(mrr_by_class: dict, f1: float, f3: float) -> float:
    """Expected leaderboard score for class shares f1, f2 = 1 - f1 - f3, f3."""
    f2 = 1.0 - f1 - f3
    return f1 * mrr_by_class["C1"] + f2 * mrr_by_class["C2"] + f3 * mrr_by_class["C3"]


def mixture_grid(mrr_by_class, f1_values=(0.15, 0.20, 0.25, 0.30), f3_values=(0.2, 0.3, 0.4)):
    rows = [(f1, f3, mixture(mrr_by_class, f1, f3)) for f1 in f1_values for f3 in f3_values if f1 + f3 < 1]
    return pd.DataFrame(rows, columns=["f1", "f3", "score"])
