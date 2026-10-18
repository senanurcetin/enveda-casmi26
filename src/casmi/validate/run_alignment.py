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


def align_runs(a: Mapping, b: Mapping, *, allow_unverified_legacy_order: bool = False) -> dict:
    """Validate evidence and align B to A by query ID.

    New-format reports must have unique IDs, per-query signatures, and the
    same ID->spectrum-signature mapping. Older runs can only be compared with
    an explicit unsafe override, and should not justify acceptance claims.
    """
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
        return {"a": ar, "b": br, "n": n, "verified": False}

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
    b_positions = {mid: pos for pos, mid in enumerate(b_ids)}
    return {
        "a": ar,
        "b": {cls: [br[cls][b_positions[mid]] for mid in a_ids] for cls in CLASSES},
        "n": n,
        "verified": True,
    }
