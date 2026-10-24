"""Aggregate-only, exact-content public-example overlap audit.

Raw identifiers and peak arrays stay in memory and are never exported. Equality
of stored numeric arrays is distinct from proof of acquisition provenance.
"""
from __future__ import annotations

from collections import defaultdict

from .spectrum_overlap import coarse_spectrum_bucket, full_spectrum_match


def _spectrum(row):
    return (row.get("precursor_mz"), row.get("ms2_mzs"),
            row.get("ms2_normalized_intensities"))


def _identifier(value):
    if value is None or not str(value).strip():
        raise ValueError("Public example identifiers must be present")
    return str(value)


class PublicOverlapAudit:
    def __init__(self, public_rows):
        self.targets = list(public_rows)
        if not self.targets:
            raise ValueError("Public example input must not be empty")
        self.buckets = defaultdict(list)
        self.molecules = defaultdict(set)
        seen = set()
        for index, row in enumerate(self.targets):
            sid = _identifier(row.get("spectrum_id"))
            molecule = _identifier(row.get("molecule_id"))
            if sid in seen:
                raise ValueError("Duplicate public spectrum identifiers")
            seen.add(sid)
            spectrum = _spectrum(row)
            if not full_spectrum_match(*spectrum, *spectrum):
                raise ValueError("Public spectrum must have finite aligned nonempty arrays")
            self.buckets[coarse_spectrum_bucket(spectrum[0], spectrum[1])].append(index)
            self.molecules[molecule].add(index)
        self.precursors = {float(row["precursor_mz"]) for row in self.targets}
        self.scanned = 0
        self.candidate_rows = 0
        self.coarse_pairs = 0
        self.exact_pairs = 0
        self.false_pairs = 0
        self.matched = set()
        self.metadata_matched = set()
        self.library_matches = defaultdict(set)
        self.metadata_disagreement_pairs = 0

    def scan_batch(self, candidate_rows, *, scanned_count=None):
        """Count all scanned rows, considering only exact-precursor candidates.

        A caller can filter Arrow batches before converting peak lists to Python.
        `scanned_count` must include filtered-out rows; finish checks full coverage.
        """
        rows = list(candidate_rows)
        count = len(rows) if scanned_count is None else scanned_count
        if isinstance(count, bool) or not isinstance(count, int) or count < len(rows):
            raise ValueError("Scanned row count must cover every considered row")
        self.scanned += count
        for row in rows:
            try:
                precursor = float(row.get("precursor_mz"))
            except (TypeError, ValueError):
                continue
            if precursor not in self.precursors:
                continue
            self.candidate_rows += 1
            spectrum = _spectrum(row)
            bucket = coarse_spectrum_bucket(spectrum[0], spectrum[1])
            for index in self.buckets.get(bucket, ()):
                self.coarse_pairs += 1
                target = self.targets[index]
                if not full_spectrum_match(*spectrum, *_spectrum(target)):
                    self.false_pairs += 1
                    continue
                library = row.get("ingest_lib")
                if not isinstance(library, str) or not library.strip():
                    raise ValueError("Matching training spectra need a source-library label")
                self.exact_pairs += 1
                self.matched.add(index)
                self.library_matches[library].add(index)
                aligned = all(target.get(key) is not None and row.get(key) == target[key]
                              for key in ("adduct", "ionization_mode"))
                if aligned:
                    self.metadata_matched.add(index)
                else:
                    self.metadata_disagreement_pairs += 1

    def finish(self, *, expected_train_rows):
        if (isinstance(expected_train_rows, bool) or not isinstance(expected_train_rows, int)
                or expected_train_rows <= 0 or self.scanned != expected_train_rows):
            raise ValueError("A complete scan of the declared training row count is required")
        return {
            "schema_version": 1,
            "status": "COMPLETE_CONTENT_AUDIT",
            "matching_method": "exact stored precursor and full ordered m/z/intensity numeric arrays; no rounding or tolerance",
            "metadata_fields": ["adduct", "ionization_mode"],
            "train_rows_scanned": self.scanned,
            "public_spectra": len(self.targets),
            "public_molecules": len(self.molecules),
            "exact_precursor_candidate_rows": self.candidate_rows,
            "coarse_candidate_pairs_after_precursor_filter": self.coarse_pairs,
            "coarse_false_match_pairs": self.false_pairs,
            "exact_matching_pairs": self.exact_pairs,
            "content_matched_spectra": len(self.matched),
            "metadata_consistent_matched_spectra": len(self.metadata_matched),
            "metadata_disagreement_pairs": self.metadata_disagreement_pairs,
            "fully_content_matched_molecules": sum(indices <= self.matched for indices in self.molecules.values()),
            "fully_metadata_consistent_molecules": sum(indices <= self.metadata_matched for indices in self.molecules.values()),
            "content_matches_by_library": {key: len(value) for key, value in sorted(self.library_matches.items())},
            "limitations": [
                "Source-library match counts overlap and must not be summed as unique public spectra.",
                "Stored numerical equality does not independently prove common acquisition identity.",
                "This audit concerns mounted public examples, not the hidden competition evaluation set.",
                "No model evaluation, training or competition submission was performed.",
            ],
        }
