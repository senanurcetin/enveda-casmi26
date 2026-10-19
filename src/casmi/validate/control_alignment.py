"""Identity-verified, diagnostic-only pairing of spectral negative controls.

This module deliberately does not accept legacy/unverified runs or produce
model-improvement ACCEPT/REJECT judgments. It verifies that the reference and
perturbed reports use the same molecule IDs, underlying original spectra,
novelty exclusion scope, model configuration and query sampling.

Separate because align_runs must continue to reject altered MS/MS evidence.
No numerical/chemistry dependency is required.
"""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy

from casmi.validate.run_alignment import CLASSES, SIGNATURE_VERSION, _check_ranks


def _signatures(run, ids, key):
    values = run.get(key)
    if not isinstance(values, list) or len(values) != len(ids):
        raise ValueError(f"{key}: one spectrum SHA-256 required per query")
    if any(not isinstance(v, str) or len(v) != 64
           or any(c not in "0123456789abcdef" for c in v) for v in values):
        raise ValueError(f"{key}: invalid SHA-256 fingerprint")
    return dict(zip(ids, values))


def _identities(run, title):
    count = run.get("n_queries")
    ids = run.get("query_ids")
    if type(count) is not int or count < 1 or not isinstance(ids, list) or len(ids) != count:
        raise ValueError(f"{title}: invalid query count or identity manifest")
    if any(not isinstance(item, str) or not item for item in ids) or len(set(ids)) != count:
        raise ValueError(f"{title}: duplicate or missing query identities")
    version = run.get("query_signature_version")
    if version != SIGNATURE_VERSION:
        raise ValueError(f"{title}: unsupported query signature recipe")
    ranks = _check_ranks(run, count)
    return ids, _signatures(run, ids, "query_signatures"), ranks


def align_control_runs(reference: Mapping, treatment: Mapping) -> dict:
    """Align diagnostic ranks by query ID, verify source fingerprints.

    `reference` must be an unperturbed validation report. `treatment` must
    be a negative-control run with original signatures captured *before*
    treatment. A changed model/config/sampling invalidates the diagnostic.
    """
    if not isinstance(reference, Mapping) or not isinstance(treatment, Mapping):
        raise ValueError("Both reports must be mappings")
    rc, tc = reference.get("spectral_control"), treatment.get("spectral_control")
    if not isinstance(rc, Mapping) or rc.get("mode") != "none":
        raise ValueError("Reference must be a signed unperturbed control=none run")
    if not isinstance(tc, Mapping) or tc.get("mode") not in ("blank", "intensity_shuffle"):
        raise ValueError("Treatment must be a signed negative-control run")
    if type(tc.get("seed")) is not int or tc["seed"] < 0:
        raise ValueError("Treatment seed must be a nonnegative integer")
    a_ids, a_signatures, a_ranks = _identities(reference, "Reference")
    b_ids, b_signatures, b_ranks = _identities(treatment, "Treatment")
    if len(a_ids) != len(b_ids) or set(a_ids) != set(b_ids):
        raise ValueError("Control population/query ID mismatch")
    source_a = reference.get("source_manifest_sha256")
    source_b = treatment.get("source_manifest_sha256")
    if source_a is not None or source_b is not None:
        if (not isinstance(source_a, str) or len(source_a) != 64
                or any(ch not in "0123456789abcdef" for ch in source_a)
                or source_a != source_b):
            raise ValueError("Negative control source artifact snapshots differ")
    if reference.get("config") != treatment.get("config"):
        raise ValueError("Control cannot be compared across model configurations")
    sampling_a, sampling_b = reference.get("query_sampling"), treatment.get("query_sampling")
    if not isinstance(sampling_a, Mapping) or not isinstance(sampling_b, Mapping):
        raise ValueError("Missing query sampling provenance")
    keep_a = {k: v for k, v in sampling_a.items()
              if k not in ("spectral_control", "control_seed")}
    keep_b = {k: v for k, v in sampling_b.items()
              if k not in ("spectral_control", "control_seed")}
    if keep_a != keep_b:
        raise ValueError("Control query sampling differs")
    if sampling_a.get("spectral_control") != "none":
        raise ValueError("Reference sampling not marked as unperturbed")
    if (sampling_b.get("spectral_control") != tc["mode"]
            or sampling_b.get("control_seed") != tc["seed"]):
        raise ValueError("Treatment mode or seed contradicts query sampling")

    source_signatures = _signatures(
        {"query_signatures": tc.get("source_query_signatures",
                                  treatment.get("spectral_control", {}).get(
                                      "source_query_signatures"))},
        b_ids, "query_signatures",
    )
    if source_signatures != a_signatures:
        raise ValueError("Control source spectra or reference-exclusion scope differs")
    if not any(b_signatures[qid] != source_signatures[qid] for qid in a_ids):
        raise ValueError("No actual spectral perturbation detected")
    index_b = {qid: i for i, qid in enumerate(b_ids)}
    return {
        "n_queries": len(a_ids),
        "control_mode": tc["mode"],
        "control_seed": tc["seed"],
        "verified_source_pairing": True,
        "changed_queries": sum(
            b_signatures[qid] != source_signatures[qid] for qid in a_ids
        ),
        "reference_ranks": deepcopy(a_ranks),
        "treatment_ranks": {
            cls: [b_ranks[cls][index_b[qid]] for qid in a_ids] for cls in CLASSES
        },
    }
