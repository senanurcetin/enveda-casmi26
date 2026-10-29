"""Strict MassBank 2.x CC0 MS2 record ingestion for external research inputs."""
from __future__ import annotations

import hashlib
import math


class RecordRejected(ValueError):
    pass


def parse_record(text):
    tags, peaks, in_peaks = {}, [], False
    if not text.rstrip().endswith("//"):
        raise RecordRejected("unterminated_record")
    for line in text.splitlines():
        if line == "//":
            break
        if line.startswith("PK$PEAK:"):
            in_peaks = True
            continue
        if in_peaks:
            fields = line.split()
            if len(fields) != 3:
                raise RecordRejected("malformed_peak_row")
            try:
                mz, intensity, relative = map(float, fields)
            except ValueError as exc:
                raise RecordRejected("nonnumeric_peak") from exc
            if not all(math.isfinite(v) for v in (mz, intensity, relative)) or mz <= 0 or intensity <= 0 or relative < 0:
                raise RecordRejected("invalid_peak")
            peaks.append((mz, intensity))
            continue
        if ": " in line:
            name, value = line.split(": ", 1)
            if name in {"AC$MASS_SPECTROMETRY", "MS$FOCUSED_ION"}:
                key, separator, value = value.partition(" ")
                if not separator:
                    raise RecordRejected("malformed_subtag")
                name += ":" + key
            tags.setdefault(name, []).append(value.strip())

    def one(name):
        values = tags.get(name, [])
        if len(values) != 1 or not values[0]:
            raise RecordRejected("missing_or_duplicate_" + name)
        return values[0]

    if one("LICENSE") != "CC0":
        raise RecordRejected("non_cc0_license")
    if one("AC$MASS_SPECTROMETRY:MS_TYPE") != "MS2" or one("AC$MASS_SPECTROMETRY:IONIZATION") != "ESI":
        raise RecordRejected("unsupported_acquisition")
    mode = one("AC$MASS_SPECTROMETRY:ION_MODE").lower()
    if mode not in ("positive", "negative"):
        raise RecordRejected("unsupported_mode")
    try:
        precursor = float(one("MS$FOCUSED_ION:PRECURSOR_M/Z"))
        declared_peaks = int(one("PK$NUM_PEAK"))
    except ValueError as exc:
        raise RecordRejected("invalid_numeric_header") from exc
    if not math.isfinite(precursor) or precursor <= 0 or declared_peaks != len(peaks) or len(peaks) < 5:
        raise RecordRejected("invalid_precursor_or_peak_count")
    if any(a[0] >= b[0] for a, b in zip(peaks, peaks[1:])):
        raise RecordRejected("unordered_or_duplicate_mz")
    maximum = max(v for _, v in peaks)
    return {"accession": one("ACCESSION"), "authors": one("AUTHORS"), "license": "CC0",
            "record_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
            "smiles": one("CH$SMILES"), "precursor_mz": precursor,
            "adduct": one("MS$FOCUSED_ION:PRECURSOR_TYPE"), "ionization_mode": mode,
            "instrument": one("AC$INSTRUMENT_TYPE"), "ms2_mzs": [m for m, _ in peaks],
            "ms2_normalized_intensities": [100 * v / maximum for _, v in peaks]}
