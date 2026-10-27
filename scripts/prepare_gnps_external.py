"""Prepare pinned, directly contributed GNPS2 standard spectra under GNPS CC0 policy."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def parse_record(record, library, source_url):
    metadata, peaks = record["metadata"], record["peaks"]
    text = " ".join(str(v) for k, v in metadata.items() if k in
                    ("INSTRUMENT", "SOURCE_INSTRUMENT", "IONSOURCE", "ION_SOURCE", "COMMENT")).lower()
    if metadata.get("MSLEVEL") != "2" or any(word in text for word in
            ("predicted", "in silico", "in-silico", "suspect", "propagated", "synthetic spectrum")):
        raise ValueError("unsupported_acquisition")
    if "esi" not in text or str(metadata.get("LIBRARYQUALITY", "1")) != "1":
        raise ValueError("unverified_acquisition_or_quality")
    mode = str(metadata.get("IONMODE", "")).strip().lower()
    if mode not in ("positive", "negative"):
        raise ValueError("invalid_mode")
    precursor = float(metadata["PEPMASS"])
    if not math.isfinite(precursor) or precursor <= 0:
        raise ValueError("invalid_precursor")
    if len(peaks) < 5 or any(len(p) != 2 or any(not math.isfinite(float(v)) or float(v) <= 0 for v in p) for p in peaks):
        raise ValueError("invalid_peaks")
    mz, intensity = zip(*[(float(m), float(i)) for m, i in peaks])
    if any(a >= b for a, b in zip(mz, mz[1:])):
        raise ValueError("unordered_peaks")
    raw_hash = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"accession": "GNPS2-SNAPSHOT:" + library + ":" + raw_hash,
            "record_sha256": raw_hash, "source_url": source_url, "library": library,
            "authors": str(metadata.get("DATACOLLECTOR", metadata.get("PI", "Not specified in source metadata"))),
            "license": "CC0", "license_basis": "GNPS direct-contribution policy; not an explicit per-record license field",
            "source_metadata_json": json.dumps(metadata, sort_keys=True),
            "smiles": metadata["SMILES"], "adduct": metadata["ADDUCT"], "ionization_mode": mode,
            "instrument": str(metadata.get("INSTRUMENT", metadata.get("SOURCE_INSTRUMENT", "unknown"))),
            "precursor_mz": precursor, "ms2_mzs": list(mz),
            "ms2_normalized_intensities": [100 * v / max(intensity) for v in intensity]}


def prepare(source_dir, output, protocol_path=None):
    import pandas as pd
    import rdkit
    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt
    from casmi.chem.canon import canonical_ik14, skeleton_key
    from casmi.chem.masses import TEST_ADDUCTS, neutral_mass, parse_adduct
    from casmi.validate.scaffold_study import family_key

    if rdkit.__version__ != "2026.03.3":
        raise ValueError("Official-metric RDKit required")
    protocol = json.loads((protocol_path or ROOT / "protocols/external-source-v2.json").read_text(encoding="utf-8"))
    rows, counts = [], Counter()
    for name, entry in protocol["external_source"]["files"].items():
        raw = (source_dir / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != entry["sha256"] or len(raw) != entry["bytes"]:
            raise ValueError("Pinned GNPS input differs")
        for record in json.loads(raw):
            counts["source_records"] += 1
            try:
                row = parse_record(record, name[:-5], entry["url"])
                if row["adduct"] not in TEST_ADDUCTS:
                    raise ValueError("unsupported_adduct")
                if (parse_adduct(row["adduct"])[-1] > 0) != (row["ionization_mode"] == "positive"):
                    raise ValueError("mode_adduct_conflict")
                mol = Chem.MolFromSmiles(row["smiles"])
                if mol is None or not mol.GetNumHeavyAtoms() or len(Chem.GetMolFrags(mol)) != 1:
                    raise ValueError("invalid_structure")
                mass = ExactMolWt(mol)
                if not 150 <= mass <= 1170:
                    raise ValueError("outside_mass_range")
                if abs(neutral_mass(row["precursor_mz"], row["adduct"]) - mass) / mass * 1e6 > 20:
                    raise ValueError("precursor_mass_conflict")
                metric, family = canonical_ik14(row["smiles"]), family_key(row["smiles"])
                if not metric or not family:
                    raise ValueError("invalid_metric_or_family")
                row.update(standard_ik14=Chem.MolToInchiKey(mol)[:14], metric_key=metric,
                           family=family, skeleton=skeleton_key(row["smiles"]), mass=mass)
                rows.append(row)
            except (ValueError, KeyError, TypeError) as exc:
                reason = str(exc) if str(exc) in {"unsupported_acquisition", "unverified_acquisition_or_quality", "invalid_mode",
                    "invalid_precursor", "invalid_peaks", "unordered_peaks", "unsupported_adduct", "mode_adduct_conflict",
                    "invalid_structure", "outside_mass_range", "precursor_mass_conflict", "invalid_metric_or_family"} else "malformed_record"
                counts[reason] += 1
    if not rows:
        raise ValueError("No eligible source records")
    output.mkdir(parents=True, exist_ok=False)
    table = pd.DataFrame(rows).drop_duplicates("accession").sort_values("accession").reset_index(drop=True)
    path = output / "gnps_cc0.parquet"
    table.to_parquet(path, index=False, compression="zstd")
    counts.update(eligible_spectra=len(table), eligible_metric_molecules=table.metric_key.nunique())
    manifest = {"schema_version": 1, "status": "PINNED_CC0_EXTERNAL_INPUT_PREPARED", "source": protocol["external_source"],
                "counts": dict(counts), "parquet_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "parquet_bytes": path.stat().st_size, "rdkit": rdkit.__version__,
                "scope": "Source preparation only; CASMI overlap and external independence not yet established"}
    (output / "external_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output / "SOURCE_ATTRIBUTION.md").write_text("# GNPS2 direct-contribution source snapshot\n\n"
        "Sources: " + ", ".join(name[:-5] for name in protocol["external_source"]["files"]) +
        ", retrieved 27 October 2026 from https://library.gnps2.org/libraries. "
        "GNPS policy assigns CC0 to directly contributed reference spectra: https://ccms-ucsd.github.io/GNPSDocumentation/browselibraries/#license. "
        "Imported third-party, propagated, MassQL-derived and suspect libraries are excluded. Original metadata, source URL and raw-record hashes "
        "are retained. Derived snapshot IDs do not claim official accession status. Some source records omit author names; missing names are disclosed. "
        "This preparation does not establish CASMI independence. No competition source records or trained weights are included.\n", encoding="utf-8")
    expanded = "GNPS-LIBRARY.json" in protocol["external_source"]["files"]
    metadata = {"id": "senanuretin/casmi26-gnps-cc0-expanded" if expanded else "senanuretin/casmi26-gnps2-cc0-snapshot",
        "title": "CASMI26 GNPS CC0 expanded sources" if expanded else "CASMI26 GNPS2 CC0 source snapshot",
        "licenses": [{"name": "CC0-1.0"}], "subtitle": "Directly contributed standard spectra for source-independence auditing",
        "description": "Pinned directly contributed GNPS standard libraries listed in the source manifest under the GNPS CC0 policy. "
            "Record hashes, original metadata, source URL and available attribution retained. See SOURCE_ATTRIBUTION.md. "
            "This is source preparation, not evidence of external independence. No CASMI competition data or model weights included."}
    (output / "dataset-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--protocol", type=Path)
    args = parser.parse_args()
    prepare(args.source_dir, args.output, args.protocol)
