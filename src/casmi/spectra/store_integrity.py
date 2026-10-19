"""Fast fail-closed structural checks for an existing CASMI MS/MS peak store.

The store format consists of contiguous float32 m/z/intensity arrays, int64
offsets, and a separate metadata Parquet file. This module reads the entire
offset vector (small relative to spectra), but NEVER scans peak values.

Layout verification confirms file lengths, offsets, row counts and per-row
peak counts. It does not establish chemical correctness, non-finite peak
content, upstream Kaggle provenance, or source-library independence.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class StoreLayout:
    """Checked store offsets with metadata-only count diagnostics."""
    offsets: np.ndarray
    spectrum_count: int
    peak_count: int
    empty_spectra: int
    max_peaks_per_spectrum: int

    def summary(self) -> dict:
        return {
            "status": "STRUCTURE_VERIFIED",
            "spectrum_count": self.spectrum_count,
            "peak_count": self.peak_count,
            "empty_spectra": self.empty_spectra,
            "max_peaks_per_spectrum": self.max_peaks_per_spectrum,
            "method": "offsets and file byte lengths; not a full peak-content scan",
        }


def inspect_store_layout(root, *, expected_spectra=None,
                         metadata_n_peaks=None) -> StoreLayout:
    """Validate raw files before mmap; raise ValueError on structural drift.

    Integers in offsets.i64 are native-endian signed int64 by store design,
    and mz.f32 / it.f32 are native-endian float32. Check byte alignment
    before trusting them; avoid Python assertions (which can be disabled).
    """
    root = Path(root)
    offsets_file = root / "offsets.i64"
    mz_file = root / "mz.f32"
    it_file = root / "it.f32"
    offset_bytes = offsets_file.stat().st_size
    mz_bytes = mz_file.stat().st_size
    it_bytes = it_file.stat().st_size
    if offset_bytes < 16 or offset_bytes % 8:
        raise ValueError("Invalid offsets.i64 byte length: at least two aligned offsets required")
    if mz_bytes % 4 or it_bytes % 4:
        raise ValueError("Spectrum peak files must have float32-aligned byte lengths")
    if mz_bytes != it_bytes:
        raise ValueError("Spectrum m/z and intensity files have different peak counts")
    offsets = np.fromfile(offsets_file, dtype=np.int64)
    n_spectra = len(offsets) - 1
    n_peaks = mz_bytes // 4
    if expected_spectra is not None:
        if type(expected_spectra) is not int or expected_spectra < 0:
            raise ValueError("Expected spectrum row count must be a nonnegative integer")
        if n_spectra != expected_spectra:
            raise ValueError("Spectrum metadata row count and offset count differ")
    if offsets[0] != 0:
        raise ValueError("First spectrum offset must equal zero")
    if np.any(offsets < 0) or np.any(offsets[1:] < offsets[:-1]):
        raise ValueError("Spectrum offsets must be nonnegative and nondecreasing")
    if int(offsets[-1]) != n_peaks:
        raise ValueError("Final spectrum offset does not equal flat peak count")
    counts = offsets[1:] - offsets[:-1]

    if metadata_n_peaks is not None:
        try:
            declared = np.asarray(metadata_n_peaks, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid metadata n_peaks values") from exc
        if declared.ndim != 1 or len(declared) != n_spectra:
            raise ValueError("Metadata n_peaks must align with spectrum rows")
        if (not np.all(np.isfinite(declared)) or np.any(declared < 0)
                or np.any(declared != np.floor(declared))):
            raise ValueError("Metadata n_peaks contains missing or invalid counts")
        if not np.array_equal(declared, counts):
            raise ValueError("Metadata n_peaks disagrees with spectrum offsets")

    return StoreLayout(
        offsets=offsets, spectrum_count=n_spectra, peak_count=n_peaks,
        empty_spectra=int(np.count_nonzero(counts == 0)),
        max_peaks_per_spectrum=int(np.max(counts, initial=0)),
    )

def scan_store_peak_content(root, *, chunk_peaks: int = 262144) -> dict:
    """Optionally scan every stored m/z and intensity using bounded memory.

    The standard store reader and inspect_store_layout deliberately avoid the
    multi-GB content scan. This separate diagnostic does not load Parquet,
    mutate peak files, or reveal individual spectra. Only aggregate counts
    leave the analysis. A "clean" result is about numeric peak admissibility,
    not chemical correctness or the origin of the source files.
    """
    if type(chunk_peaks) is not int or chunk_peaks <= 0:
        raise ValueError("chunk_peaks must be a positive integer")
    root = Path(root)
    layout = inspect_store_layout(root)
    nonfinite_mz = nonpositive_mz = 0
    nonfinite_intensity = nonpositive_intensity = 0
    invalid_peak_pairs = 0
    read = 0

    with (root / "mz.f32").open("rb") as mz_file, (root / "it.f32").open("rb") as it_file:
        remaining = layout.peak_count
        while remaining:
            count = min(chunk_peaks, remaining)
            mz = np.fromfile(mz_file, dtype=np.float32, count=count)
            it = np.fromfile(it_file, dtype=np.float32, count=count)
            if len(mz) != count or len(it) != count:
                raise ValueError("Peak file truncated during content scan")
            mz_finite = np.isfinite(mz)
            it_finite = np.isfinite(it)
            mz_nonpositive = mz_finite & (mz <= 0)
            it_nonpositive = it_finite & (it <= 0)
            nonfinite_mz += int(np.count_nonzero(~mz_finite))
            nonpositive_mz += int(np.count_nonzero(mz_nonpositive))
            nonfinite_intensity += int(np.count_nonzero(~it_finite))
            nonpositive_intensity += int(np.count_nonzero(it_nonpositive))
            invalid_peak_pairs += int(np.count_nonzero(
                ~(mz_finite & it_finite) | mz_nonpositive | it_nonpositive
            ))
            read += count
            remaining -= count

        if mz_file.read(1) or it_file.read(1):
            raise ValueError("Peak file grew during content scan")

    invalid = invalid_peak_pairs > 0
    return {
        "status": "INVALID_CLEANED_PEAKS" if invalid else "VALID_CLEANED_PEAKS",
        "total_peaks_scanned": read,
        "nonfinite_mz": nonfinite_mz,
        "nonpositive_mz": nonpositive_mz,
        "nonfinite_intensity": nonfinite_intensity,
        "nonpositive_intensity": nonpositive_intensity,
        "invalid_peak_pairs": invalid_peak_pairs,
        "chunk_peaks": chunk_peaks,
        "method": "full numeric content scan in bounded float32 chunks",
        "interpretation": (
            "Aggregate cleaned-peak numeric validity only; not chemistry, "
            "data-lineage, leakage or instrument accuracy"
        ),
    }
