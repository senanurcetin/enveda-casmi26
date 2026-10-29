"""Deterministic external-family selection and whole-reference training exclusion."""
from __future__ import annotations

import hashlib


def select_external_families(source, salt, *, maximum=150, minimum=30):
    """Choose source chemistry before candidate or model inspection.

    Every family annotation of a chosen metric is protected, including variants
    outside the five spectra eventually used as query input.
    """
    if type(maximum) is not int or type(minimum) is not int or not 1 <= minimum <= maximum:
        raise ValueError("Positive bounded external population required")
    required = {"family", "metric_key", "accession", "standard_ik14"}
    if not required <= set(source) or source[list(required)].isna().any().any() or source.accession.duplicated().any():
        raise ValueError("Invalid external source identity contract")
    family_of = source.groupby("metric_key").family.agg(lambda values: set(values)).to_dict()
    chosen, blocked_families = [], set()
    ordered = sorted(set(source.family), key=lambda family: hashlib.sha256((salt + "|external-family|" + family).encode()).hexdigest())
    for family in ordered:
        if family in blocked_families:
            continue
        for metric in sorted(set(source.loc[source.family == family, "metric_key"])):
            if not family_of[metric] & blocked_families:
                chosen.append(metric); blocked_families.update(family_of[metric]); break
        if len(chosen) == maximum:
            break
    if len(chosen) < minimum:
        raise ValueError("Insufficient independently grouped external molecules")
    return source[source.metric_key.isin(chosen)].sort_values("accession").reset_index(drop=True), sorted(blocked_families)


def blocked_training_identities(rows, selected_source, selected_families, *, return_audit=False):
    """Inspect full labels; a blocked label blocks every row of its declared key."""
    from casmi.chem.canon import canonical_ik14, skeleton_key
    from casmi.validate.scaffold_study import family_key

    external_metrics = set(selected_source.metric_key)
    external_standard = set(selected_source.standard_ik14)
    external_skeletons = {skeleton_key(s) for s in selected_source.smiles}
    selected_families = set(selected_families)
    blocked, cache = set(external_standard), {}
    unassigned_labels, unassigned_keys = set(), set()
    for key, smiles in rows:
        if not isinstance(key, str) or not key:
            raise ValueError("Missing original training identity")
        if smiles not in cache:
            family = family_key(smiles)
            if family is None:
                # Undefined grouping cannot establish family separation. Remove
                # every row of this declared identity rather than retain it.
                unassigned_labels.add(smiles)
                cache[smiles] = True
            else:
                by_family = family in selected_families
                skeleton = skeleton_key(smiles)
                if skeleton is None:
                    raise ValueError("Unparseable original training connectivity")
                by_metric = False
                if skeleton in external_skeletons:
                    metric = canonical_ik14(smiles)
                    if metric is None:
                        raise ValueError("Unverifiable matching original metric")
                    by_metric = metric in external_metrics
                cache[smiles] = by_family or by_metric
        if smiles in unassigned_labels:
            unassigned_keys.add(key)
        if key in external_standard or cache[smiles]:
            blocked.add(key)
    audit = {"unassignable_family_labels_removed": len(unassigned_labels),
             "unassignable_family_training_identities_removed": len(unassigned_keys)}
    return (blocked, audit) if return_audit else blocked
