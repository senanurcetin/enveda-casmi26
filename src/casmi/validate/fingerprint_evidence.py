"""Standard-library D2 fingerprint evidence gate for scientific run comparisons.

A verified local input inventory alone does not prove that an analog model
used Morgan fingerprints computed for the *same SMILES and masses*. Reports
without a verified D2 sidecar remain descriptive/exploratory; they must not
support a model ACCEPT/REJECT decision.

This validates report evidence, not the source data itself. The upstream
run_validation step must cryptographically verify the actual artifacts.
"""
from __future__ import annotations

from collections.abc import Mapping

FINGERPRINT_SCHEMA = "casmi-morgan-fingerprint-provenance-v1"
FINGERPRINT_STATUS = "POOL_FINGERPRINTS_VERIFIED"


def _valid_sha(value) -> bool:
    return (isinstance(value, str) and len(value) == 64
            and all(char in "0123456789abcdef" for char in value))


def verified_d2_fingerprint_evidence(a: Mapping, b: Mapping) -> bool:
    """Return whether all D2 runs use verified, mutually compatible vectors.

    D1-only records do not consume D2 fingerprints. A D1/D2 ablation is
    permitted if the D2 run has verified provenance and BOTH runs share a
    complete encompassing input manifest as enforced by run_alignment.
    """
    from casmi.validate.analog_index_evidence import verified_analog_index_evidence

    proofs = []
    missing = False
    for run in (a, b):
        cfg = run.get("config")
        channels = cfg.get("channels") if isinstance(cfg, Mapping) else None
        if not isinstance(channels, list):
            raise ValueError("Missing model channel provenance")
        if "d2" not in channels:
            continue
        proof = run.get("fingerprint_provenance")
        if proof is None:
            missing = True  # historical report: descriptive only
            continue
        if not isinstance(proof, Mapping):
            raise ValueError("Malformed D2 fingerprint provenance")
        if (proof.get("status") != FINGERPRINT_STATUS
                or proof.get("schema") != FINGERPRINT_SCHEMA):
            raise ValueError("Unverified D2 fingerprint provenance schema or status")
        count = proof.get("n_structures")
        if type(count) is not int or count < 1:
            raise ValueError("Invalid D2 fingerprint row count")
        for key in ("pool_rows_sha256", "fingerprint_bytes_sha256"):
            if not _valid_sha(proof.get(key)):
                raise ValueError(f"Invalid D2 fingerprint provenance: {key}")
        proofs.append(proof)
        # A valid fingerprint sidecar cannot certify a cached D2 reference
        # matrix. Unverified historical/cached runs remain exploratory.
        if not verified_analog_index_evidence(run):
            missing = True

    if missing:
        return False
    if len(proofs) == 2 and any(
        proofs[0][key] != proofs[1][key] for key in
        ("n_structures", "pool_rows_sha256", "fingerprint_bytes_sha256")
    ):
        raise ValueError("D2 fingerprint source mismatch between validation runs")
    return True
