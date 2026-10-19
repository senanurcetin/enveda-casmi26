"""Fail-closed identity checks for paired CASMI validation comparisons.

Paired confidence intervals require exactly the same query molecules and the
same source spectra, not merely the same number of queries. Unlike rank
summaries, this module needs no numerical or chemistry dependencies.
"""
from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence

CLASSES = ("C1", "C2", "C3")
SIGNATURE_VERSION = "casmi-query-spectra-v2"


def query_signature(query) -> str:
    """Fingerprint *actual* source spectra and exclusion rows for one query.

    Signatures are validation provenance, not a tamper-proof dataset checksum.
    This prevents pairing runs that resample spectra, change source-library
    exclusion scope, or change store row order. Version the signature recipe.
    """
    payload = {
        "version": SIGNATURE_VERSION,
        "qid": str(query.qid),
        "pool_key": str(query.pool_key),
        "selected_query_rows": sorted(int(i) for i in query.spectrum_rows),
        "source_library_excluded_rows": (
            sorted(int(i) for i in query.source_library_rows)
            if getattr(query, "source_library_rows", None) is not None else None
        ),
        "spectra": [
            {
                "mz": [float(x).hex() for x in spec["mz"]],
                "intensity": [float(x).hex() for x in spec["it"]],
                "precursor_mz": float(spec["precursor_mz"]).hex(),
                "adduct": str(spec["adduct"]),
                "ce_ev": float(spec["ce_ev"]).hex(),
                "mode": str(spec["mode"]),
            }
            for spec in query.spectra
        ],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _check_ranks(run: Mapping, n: int) -> dict[str, list[int]]:
    classes = run.get("classes")
    if not isinstance(classes, dict):
        raise ValueError("Missing classes in validation run")
    result = {}
    for cls in CLASSES:
        block = classes.get(cls)
        if not isinstance(block, dict):
            raise ValueError(f"Missing class {cls}")
        ranks = block.get("ranks")
        if not isinstance(ranks, list) or len(ranks) != n:
            raise ValueError(f"{cls}: expected {n} per-query ranks")
        for rank in ranks:
            if type(rank) is not int or not 0 <= rank <= 25:
                raise ValueError(f"{cls}: rank must be an integer from 0 to 25")
        computed = sum(1.0 / rank for rank in ranks if rank) / n
        recorded = block.get("mrr")
        if recorded is not None and (
            isinstance(recorded, bool) or not isinstance(recorded, (float, int))
            or not math.isfinite(recorded)
            or abs(computed - recorded) > 0.000005
        ):
            raise ValueError(f"{cls}: MRR summary disagrees with per-query ranks")
        result[cls] = ranks
    return result



def _source_snapshot_verified(a: Mapping, b: Mapping) -> bool:
    """Check local input inventory identity and complete known input coverage.

    Historical reports without inventory hashes remain exploratory. Merely
    naming an identical SHA-256 does not prove that a subset manifest covers
    the store and active candidate pool used by either run.
    """
    from casmi.validate.source_manifest import (
        MANIFEST_SCHEMA, required_validation_inputs,
    )

    hashes = (a.get("source_manifest_sha256"), b.get("source_manifest_sha256"))
    if hashes == (None, None):
        return False
    if None in hashes:
        raise ValueError("Source artifact snapshot mismatch: one run lacks a manifest")
    if any(not isinstance(s, str) or len(s) != 64 or
           any(ch not in "0123456789abcdef" for ch in s) for s in hashes):
        raise ValueError("Malformed source artifact snapshot SHA-256")
    if hashes[0] != hashes[1]:
        raise ValueError("Source artifact snapshot mismatch: input files differ")

    for run in (a, b):
        # Old signed reports did not yet validate required-file coverage.
        # They can be inspected, but never yield conclusive model verdicts.
        if run.get("source_manifest_coverage") != "required-validation-inputs-v1":
            return False
        if run.get("source_manifest_schema") != MANIFEST_SCHEMA:
            raise ValueError("Invalid source artifact manifest schema")
        pool_tag = run.get("source_manifest_pool_tag")
        channels = run.get("config", {}).get("channels") if isinstance(run.get("config"), Mapping) else None
        if not isinstance(channels, list) or "d1" not in channels:
            raise ValueError("Missing validation channel provenance")
        expected = required_validation_inputs(
            pool_tag=pool_tag, analog=("d2" in channels),
        )
        provided = run.get("source_manifest_required_inputs")
        if not isinstance(provided, list) or provided != list(expected):
            raise ValueError("Incomplete or contradictory source artifact inventory")
        count = run.get("source_manifest_file_count")
        if type(count) is not int or count < len(expected):
            raise ValueError("Invalid source artifact file count")
    from casmi.validate.fingerprint_evidence import verified_d2_fingerprint_evidence

    # An unchanged pool/fingerprint FILE inventory does not establish that
    # D2 vectors belong to the pool's exact SMILES, masses and row ordering.
    # Historical D2 runs without a verified sidecar remain exploratory.
    return verified_d2_fingerprint_evidence(a, b)


def align_runs(a: Mapping, b: Mapping, *, allow_unverified_legacy_order: bool = False) -> dict:
    """Validate evidence and align B to A by query ID.

    New-format reports must have unique IDs, per-query signatures, and the
    same ID->spectrum-signature mapping. Older runs can only be compared with
    an explicit unsafe override, and should not justify acceptance claims.
    """
    # Paired improvement inference requires identical *actual* MS/MS inputs.
    # Negative controls intentionally change spectra, so they must never
    # receive model-improvement ACCEPT/REJECT judgments from this comparator.
    for run in (a, b):
        control = run.get("spectral_control")
        if control is not None:
            if not isinstance(control, Mapping):
                raise ValueError("Malformed spectral-control provenance")
            if control.get("mode", "none") != "none":
                raise ValueError(
                    "Negative-control reports are diagnostic only; cannot "
                    "be compared as model-improvement runs"
                )

    n_a, n_b = a.get("n_queries"), b.get("n_queries")
    if type(n_a) is not int or type(n_b) is not int or n_a < 1 or n_a != n_b:
        raise ValueError("Runs must contain the same positive number of queries")
    n = n_a
    ar, br = _check_ranks(a, n), _check_ranks(b, n)
    a_ids, b_ids = a.get("query_ids"), b.get("query_ids")
    a_sig, b_sig = a.get("query_signatures"), b.get("query_signatures")
    has_manifest = all(isinstance(x, list) for x in (a_ids, b_ids, a_sig, b_sig))
    if not has_manifest:
        if not allow_unverified_legacy_order:
            raise ValueError(
                "Cannot validate paired query identities/spectra in legacy runs. "
                "Regenerate reports with query_ids + query_signatures; "
                "use --allow-unverified-legacy-order only for exploratory comparison."
            )
        return {"a": ar, "b": br, "n": n, "verified": False,
                "source_verified": False, "decision_eligible": False}

    for run_ids, run_sigs, name in ((a_ids, a_sig, "A"), (b_ids, b_sig, "B")):
        if len(run_ids) != n or len(run_sigs) != n:
            raise ValueError(f"Run {name}: manifest length differs from n_queries")
        if any(not isinstance(i, str) or not i for i in run_ids):
            raise ValueError(f"Run {name}: empty or invalid query ID")
        if len(set(run_ids)) != n:
            raise ValueError(f"Run {name}: duplicate query IDs")
        if any(not isinstance(s, str) or len(s) != 64 or
               any(c not in "0123456789abcdef" for c in s) for s in run_sigs):
            raise ValueError(f"Run {name}: malformed query signatures")

    # A changed signature recipe makes otherwise comparable digests
    # uninterpretable; old records without version tags remain supported.
    version_a = a.get("query_signature_version")
    version_b = b.get("query_signature_version")
    if version_a != version_b:
        raise ValueError("Query signature algorithm versions differ")
    if version_a is not None and (not isinstance(version_a, str) or not version_a):
        raise ValueError("Malformed query signature version")

    if set(a_ids) != set(b_ids):
        raise ValueError("Query identity mismatch: validation populations differ")
    a_by_id = dict(zip(a_ids, a_sig))
    b_by_id = dict(zip(b_ids, b_sig))
    if any(a_by_id[mid] != b_by_id[mid] for mid in a_ids):
        raise ValueError("Spectrum/knockout mismatch: query inputs differ")
    source_verified = _source_snapshot_verified(a, b)
    b_positions = {mid: pos for pos, mid in enumerate(b_ids)}
    return {
        "a": ar,
        "b": {cls: [br[cls][b_positions[mid]] for mid in a_ids] for cls in CLASSES},
        "n": n,
        "verified": True,
        "source_verified": source_verified,
        "decision_eligible": source_verified,
    }
