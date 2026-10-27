"""Prepare pinned, explicitly CC0 external MS2 inputs; no CASMI data is read."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.validate.massbank_records import parse_record, RecordRejected


def prepare(archive, output):
    import pandas as pd
    import rdkit
    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt
    from casmi.chem.canon import canonical_ik14, skeleton_key
    from casmi.chem.masses import TEST_ADDUCTS, neutral_mass, parse_adduct
    from casmi.validate.scaffold_study import family_key

    if rdkit.__version__ != "2026.03.3":
        raise ValueError("Pinned official-metric RDKit required")
    protocol = json.loads((ROOT / "protocols/external-c3-v1.json").read_text(encoding="utf-8"))
    with archive.open("rb") as stream:
        archive_hash = hashlib.file_digest(stream, "sha256").hexdigest()
    if archive_hash != protocol["external_source"]["archive_sha256"]:
        raise ValueError("Archive differs from registered official source")
    rows, counts, accession_hashes = [], Counter(), {}
    with zipfile.ZipFile(archive) as source:
        for info in source.infolist():
            if not info.filename.endswith(".txt"):
                continue
            counts["source_records"] += 1
            if info.file_size > 2_000_000:
                raise ValueError("Unbounded source record")
            text = source.read(info).decode("utf-8", errors="strict")
            if "\nLICENSE: CC0\n" not in text:
                counts["non_cc0_records"] += 1
                continue
            try:
                row = parse_record(text)
                if row["adduct"] not in TEST_ADDUCTS:
                    raise RecordRejected("unsupported_adduct")
                if (parse_adduct(row["adduct"])[-1] > 0) != (row["ionization_mode"] == "positive"):
                    raise RecordRejected("mode_adduct_conflict")
                mol = Chem.MolFromSmiles(row["smiles"])
                if mol is None or not mol.GetNumHeavyAtoms() or len(Chem.GetMolFrags(mol)) != 1:
                    raise RecordRejected("invalid_or_disconnected_structure")
                mass = ExactMolWt(mol)
                if not 150 <= mass <= 1170:
                    raise RecordRejected("outside_mass_range")
                measured_mass = neutral_mass(row["precursor_mz"], row["adduct"])
                if abs(measured_mass - mass) / mass * 1e6 > 20:
                    raise RecordRejected("precursor_mass_conflict")
                metric, family = canonical_ik14(row["smiles"]), family_key(row["smiles"])
                if not metric or not family:
                    raise RecordRejected("invalid_metric_or_family")
                if row["accession"] in accession_hashes:
                    if accession_hashes[row["accession"]] != row["record_sha256"]:
                        raise ValueError("Conflicting source accession")
                    counts["duplicate_source_record"] += 1
                    continue
                accession_hashes[row["accession"]] = row["record_sha256"]
                row.update(standard_ik14=Chem.MolToInchiKey(mol)[:14], metric_key=metric,
                           family=family, skeleton=skeleton_key(row["smiles"]), mass=mass,
                           source_path=info.filename.split("/", 1)[-1],
                           source_url=f"https://github.com/MassBank/MassBank-data/blob/{protocol['external_source']['commit']}/{info.filename.split('/', 1)[-1]}")
                rows.append(row)
            except RecordRejected as exc:
                counts[str(exc)] += 1
    if not rows:
        raise ValueError("No eligible CC0 records")
    output.mkdir(parents=True, exist_ok=False)
    table = pd.DataFrame(rows).sort_values("accession").reset_index(drop=True)
    path = output / "massbank_cc0.parquet"
    table.to_parquet(path, index=False, compression="zstd")
    with path.open("rb") as stream:
        table_hash = hashlib.file_digest(stream, "sha256").hexdigest()
    counts.update(eligible_spectra=len(table), eligible_metric_molecules=table.metric_key.nunique())
    manifest = {"schema_version": 1, "status": "PINNED_CC0_EXTERNAL_INPUT_PREPARED",
                "source": protocol["external_source"], "counts": dict(counts),
                "parquet_sha256": table_hash, "parquet_bytes": path.stat().st_size,
                "rdkit": rdkit.__version__,
                "scope": "Source preparation only; CASMI overlap and external independence not yet established"}
    (output / "external_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output / "SOURCE_ATTRIBUTION.md").write_text(
        "# MassBank CC0 external inputs\n\nOnly records explicitly licensed CC0 are retained. "
        "Each row preserves original authors, accession, license, source URL and raw-record SHA-256. "
        "No CASMI training/competition records are included.\n\n"
        "Source: https://github.com/MassBank/MassBank-data/tree/" + protocol["external_source"]["commit"] +
        "\n\nPreparation filters are registered in protocols/external-c3-v1.json in the CASMI research repository. "
        "This input preparation does not establish non-overlap with CASMI.\n", encoding="utf-8")
    metadata = {"id": "senanuretin/casmi26-massbank-cc0-202603", "title": "CASMI26 MassBank CC0 2026.03",
                "licenses": [{"name": "CC0-1.0"}],
                "subtitle": "Pinned CC0 MS2 source records for independent-overlap auditing",
                "description": "Selected explicitly CC0 records from the official MassBank 2026.03 release, commit "
                    + protocol["external_source"]["commit"] + ". Source attribution, accession, URL and record hashes retained per row. "
                    "No competition data, trained weights or CASMI-derived source rows are present. Preparation is not proof of external independence; "
                    "the CASMI research study performs a separate full-training overlap audit. See SOURCE_ATTRIBUTION.md and the input manifest."}
    (output / "dataset-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.archive, args.output)
