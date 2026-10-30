"""Outcome-blind follow-on populations, preserving previously opened test families."""
from __future__ import annotations

import hashlib
import numpy as np
import pandas as pd

from .scaffold_study import family_key, partition, digest
from casmi.chem.canon import canonical_ik14


def prepare_internal(meta, truth_pool, previous_protocol, previous_audit, protocol, *, protected=None, return_protected=False, return_candidates=False):
    keys = meta.train_ik14.astype(str).to_numpy()
    labels = dict(zip(truth_pool.key, truth_pool.smiles))
    masses = dict(zip(truth_pool.key, truth_pool.mass))
    eligible = (meta.ingest_lib.isin(previous_protocol["query_libraries"]).to_numpy()
                & meta.test_adduct.to_numpy() & (meta.n_peaks.to_numpy() >= 5)
                & (np.abs(meta.ppm.to_numpy()) <= 20))
    rows = np.flatnonzero(eligible)
    choices = pd.DataFrame({"row": rows, "key": keys[rows], "library": meta.ingest_lib.astype(str).to_numpy()[rows]})
    choices = choices[choices.key.map(masses).between(*previous_protocol["query_mass_range"])]
    counts = choices.groupby(["key", "library"]).size().reset_index(name="count")
    best = counts.sort_values(["count", "key", "library"], ascending=[False, True, True]).drop_duplicates("key")
    groups = {key: group.row.to_numpy() for key, group in choices.merge(best[["key", "library"]], on=["key", "library"]).groupby("key")}
    historical_vin = {family_key(labels[key]) for key in pd.unique(keys[(meta.ingest_lib == "enveda-np-examples").to_numpy()]) if key in labels}
    if protected is not None:
        if protected["previous_partitions"] != {name: previous_audit[name] for name in ("train", "dev", "test")}:
            raise ValueError("Pre-purge historical partition evidence differs")
        historical_vin.update(protected["historical_vin"])
    representatives = {}
    for row in best.itertuples():
        family = family_key(labels[row.key])
        if family and family not in historical_vin and (family not in representatives or row.key < representatives[family]["key"]):
            representatives[family] = {"key": row.key, "family": family, "library": row.library}
    records = []
    for value in representatives.values():
        metric = canonical_ik14(labels[value["key"]])
        if metric:
            records.append({**value, "metric": metric})
    if protected is None:
        previous = partition(records, previous_protocol["split_counts"], previous_protocol["split_salt"])
        for name, values in previous.items():
            for field, private_values in (("identity_sha256", [v["key"] for v in values]),
                                          ("family_sha256", [v["group"] for v in values]),
                                          ("metric_sha256", [v["metric"] for v in values])):
                if digest(private_values) != previous_audit[name][field]:
                    raise ValueError("Historical partitions could not be reconstructed exactly")
        old_values = [v for values in previous.values() for v in values]
        old_metrics = {v["metric"] for v in old_values}
        old_families = {v[field] for v in old_values for field in ("family", "group")}
        protected = {"old_metrics": sorted(old_metrics), "old_families": sorted(old_families),
                     "historical_vin": sorted(v for v in historical_vin if v),
                     "previous_partitions": {name: previous_audit[name] for name in ("train", "dev", "test")}}
    else:
        old_metrics, old_families = set(protected["old_metrics"]), set(protected["old_families"])
    fresh_records = [v for v in records if v["metric"] not in old_metrics and v["family"] not in old_families]
    if return_candidates:
        # No partition selection, candidate retrieval or score inspection.
        return fresh_records, groups, labels, protected
    splits = partition(fresh_records, protocol["internal_counts"], protocol["split_salt"])
    audit = {}
    for name, values in splits.items():
        audit[name] = {"n_molecules": len(values), "identity_sha256": digest([r["key"] for r in values]),
                       "family_sha256": digest([r["group"] for r in values]), "metric_sha256": digest([r["metric"] for r in values]),
                       "libraries": pd.Series([r["library"] for r in values]).value_counts().to_dict()}
    selected = [v for values in splits.values() for v in values]
    metric_overlap = len({v["metric"] for v in selected} & old_metrics)
    family_overlap = len({v[field] for v in selected for field in ("family", "group")} & old_families)
    cross_split_overlap = sum(len({v["key"] for v in splits[a]} & {v["key"] for v in splits[b]})
                              for a, b in (("train", "dev"), ("train", "test"), ("dev", "test")))
    if metric_overlap or family_overlap or cross_split_overlap:
        raise ValueError("New population overlaps a protected historical or current partition")
    audit.update(previous_population_reconstructed=True, previous_selected_metric_overlap=metric_overlap,
                 previous_selected_family_overlap=family_overlap, cross_split_identity_overlaps=cross_split_overlap)
    result = (splits, groups, labels, audit)
    return (*result, protected) if return_protected else result


def select_external(novel, splits, salt, *, maximum=150, minimum=30):
    blocked = {v[field] for name in ("train", "dev") for v in splits[name] for field in ("family", "group")}
    excluded_metrics = set(novel.loc[novel.family.isin(blocked), "metric_key"])
    eligible = novel[~novel.metric_key.isin(excluded_metrics)]
    ordered = sorted(set(eligible.metric_key), key=lambda metric: hashlib.sha256((salt + "|external|" + metric).encode()).hexdigest())[:maximum]
    if len(ordered) < minimum:
        raise ValueError("Insufficient independent external molecules after supervised-family exclusion")
    selected = []
    for metric in ordered:
        spectra = eligible[eligible.metric_key == metric].sort_values("accession").head(5)
        first = spectra.iloc[0]
        selected.append({"metric": metric, "key": first.standard_ik14, "smiles": first.smiles,
                         "family": first.family, "records": spectra})
    overlap = len({v["family"] for v in selected} & blocked)
    if overlap:
        raise ValueError("External population shares supervised training/development families")
    audit = {"n_molecules": len(selected), "metric_sha256": digest(ordered),
             "family_sha256": digest([v["family"] for v in selected]),
             "accessions_sha256": digest([str(a) for v in selected for a in v["records"].accession]),
             "supervised_family_overlap": overlap, "excluded_supervised_family_molecules": len(excluded_metrics)}
    return selected, audit
