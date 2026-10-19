"""Decision-grade D2 index evidence for signed CASMI validation reports.

Legacy cached analog index files were created without source provenance.
A run without a *freshly* rebuilt reference matrix remains exploratory even
when the source manifest and candidate Morgan bits have been verified.
"""
from __future__ import annotations

from collections.abc import Mapping

INDEX_SCHEMA = "casmi-d2-fresh-index-v1"
INDEX_STATUS = "REBUILT_WITHOUT_CACHE"


def verified_analog_index_evidence(run: Mapping) -> bool:
    """Validate reported cache-free index provenance, not external data truth."""
    proof = run.get("analog_index_provenance")
    if proof is None:
        return False
    if not isinstance(proof, Mapping):
        raise ValueError("Malformed D2 index provenance")
    if proof.get("schema") != INDEX_SCHEMA or proof.get("status") != INDEX_STATUS:
        raise ValueError("Unverified D2 index provenance schema or status")
    if proof.get("cache_used") is not False:
        raise ValueError("D2 index provenance must explicitly bypass the cache")
    digest = proof.get("index_sha256")
    if not isinstance(digest, str) or len(digest) != 64 or any(
        char not in "0123456789abcdef" for char in digest
    ):
        raise ValueError("Invalid D2 index content digest")
    rows = proof.get("n_references")
    nnz = proof.get("n_nonzero")
    if type(rows) is not int or rows < 1:
        raise ValueError("Invalid D2 index reference count")
    if type(nnz) is not int or nnz < 0:
        raise ValueError("Invalid D2 index nonzero count")
    cfg = run.get("config")
    analog = cfg.get("analog") if isinstance(cfg, Mapping) else None
    policy = analog.get("reference_policy") if isinstance(analog, Mapping) else None
    if not isinstance(policy, str) or policy not in ("legacy", "source_balanced"):
        raise ValueError("Missing or unsupported D2 reference-selection policy")
    if proof.get("reference_policy") != policy:
        raise ValueError("D2 index policy contradicts model configuration")
    return True
