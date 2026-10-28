"""Create a new training snapshot with every selected external family removed.

Original inputs stay read-only. Selected public source rows and historical
exclusion labels remain in private scratch storage; export aggregates only.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]


def prepare():
    import pandas as pd
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from casmi.spectra.store import SpectrumStore
    from casmi.pool.pool import CandidatePool
    from casmi.validate.external_population import prepare_internal
    from casmi.validate.external_purge import select_external_families, blocked_training_identities
    from casmi.validate.scaffold_study import digest
    from casmi.validate.source_manifest import sha256_file
    from scripts.run_external_audit import audit

    protocol_path = ROOT / "protocols/purged-external-v1.json"
    protocol = json.loads(protocol_path.read_text())
    original = Path(os.environ["CASMI_ORIGINAL_TRAIN"])
    source_path = Path(os.environ["CASMI_ORIGINAL_EXTERNAL"])
    source_manifest = json.loads(Path(os.environ["CASMI_ORIGINAL_EXTERNAL_MANIFEST"]).read_text())
    previous = json.loads((ROOT / "reports/advanced_study.json").read_text())
    expected_train = previous["raw_inputs"]["train.parquet"]
    original_hash = sha256_file(original)
    if original_hash != expected_train["sha256"] or original.stat().st_size != expected_train["bytes"]:
        raise ValueError("Original training snapshot differs from protected historical input")
    if (source_manifest["source"] != protocol["external_source"]
            or sha256_file(source_path) != source_manifest["parquet_sha256"]):
        raise ValueError("Pinned external source or source policy differs")
    artifacts = Path(os.environ["CASMI_ARTIFACTS"])
    original_store = SpectrumStore(artifacts / "store")
    truth_pool = CandidatePool(pd.read_parquet(artifacts / "structures.parquet"))
    previous_protocol = json.loads((ROOT / "protocols/advanced-internal-v1.json").read_text())
    *_, protected = prepare_internal(original_store.meta, truth_pool, previous_protocol,
        previous["partitions"], protocol, return_protected=True)
    protected["original_training_sha256"] = original_hash
    print("Prepared verified original historical population exclusions before source filtering", flush=True)
    source = pd.read_parquet(source_path)
    selected, families = select_external_families(source, protocol["split_salt"])
    print(f"Prepared {selected.metric_key.nunique()} external molecules before fitting or candidate inspection", flush=True)
    unique_labels = pd.read_parquet(original, columns=["inchikey14", "normalized_smiles"]).drop_duplicates()
    blocked = blocked_training_identities(unique_labels.itertuples(index=False, name=None), selected, families)
    actual_blocked = set(unique_labels.inchikey14) & blocked
    output = Path(os.environ["CASMI_PURGED_DATA"])
    output.mkdir(parents=True, exist_ok=False)
    protected_path = output / "protected_population.json"
    protected_path.write_text(json.dumps(protected, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    selected_path = output / "selected_external.parquet"
    selected.to_parquet(selected_path, index=False, compression="zstd")
    counts = dict(source_manifest["counts"])
    counts.update(selection_excluded_eligible_spectra=len(source) - len(selected),
                  eligible_spectra=len(selected), eligible_metric_molecules=int(selected.metric_key.nunique()))
    derived_manifest = {**source_manifest, "counts": counts, "parquet_sha256": sha256_file(selected_path),
                        "parquet_bytes": selected_path.stat().st_size,
                        "scope": "Outcome-blind source-family selection; independence is relative to a separately purged training snapshot"}
    manifest_path = output / "selected_external_manifest.json"
    manifest_path.write_text(json.dumps(derived_manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    train = pq.ParquetFile(original)
    retained, removed, seen = 0, 0, 0
    clean_path = output / "train.parquet"
    blocked_array = pa.array(sorted(blocked), type=pa.string())
    with pq.ParquetWriter(clean_path, train.schema_arrow, compression="zstd") as writer:
        for index, batch in enumerate(train.iter_batches(batch_size=16384)):
            mask = pc.invert(pc.is_in(batch.column("inchikey14"), value_set=blocked_array))
            clean = batch.filter(mask)
            if clean.num_rows:
                writer.write_table(pa.Table.from_batches([clean]))
            retained += clean.num_rows; removed += batch.num_rows - clean.num_rows; seen += batch.num_rows
            if index % 25 == 0:
                print(f"Prepared purged reference rows {seen}/{train.metadata.num_rows}", flush=True)
    if retained < 100000 or retained + removed != train.metadata.num_rows or removed == 0:
        raise ValueError("Purged training row accounting or minimum size failed")
    independence, novel = audit(clean_path, selected_path, manifest_path, ROOT / "protocols/external-source-v3.json")
    identities = independence["identity_audit"]
    if (independence["status"] != "EXTERNAL_INDEPENDENCE_FEASIBLE"
            or identities["excluded_external_molecules_union"] != 0 or len(novel) != len(selected)):
        raise ValueError("Selected external source retains training identity/content overlap after purge")
    remaining = unique_labels[~unique_labels.inchikey14.isin(blocked)]
    from casmi.validate.scaffold_study import family_key
    retained_family_overlap = sum(family_key(s) in set(families) for s in set(remaining.normalized_smiles))
    if retained_family_overlap:
        raise ValueError("Selected external family remains in filtered training")
    report = {"schema_version": 1, "status": "EXTERNAL_FAMILIES_PURGED_AND_INDEPENDENCE_VERIFIED",
        "protocol_sha256": sha256_file(protocol_path), "original_training": expected_train,
        "original_source_preparation": source_manifest, "selected_source_preparation": derived_manifest,
        "original_training_rows": train.metadata.num_rows, "retained_training_rows": retained, "removed_training_rows": removed,
        "removed_training_identities": len(actual_blocked), "selected_external_molecules": int(selected.metric_key.nunique()),
        "selected_external_families": len(families), "selected_external_spectra_available": len(selected),
        "selected_metric_sha256": digest(sorted(set(selected.metric_key))), "selected_family_sha256": digest(families),
        "selected_accessions_sha256": digest(sorted(selected.accession)), "retained_external_family_overlap": retained_family_overlap,
        "protected_population_sha256": sha256_file(protected_path), "historical_population_reconstructed_before_purge": True,
        "post_purge_independence": independence, "model_performance_inspected": False,
        "scope": "External-family holdout relative to newly purged training/reference data; original snapshot novelty feasibility remains unsuccessful"}
    (output / "purge_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("Prepared independently audited family-purged training snapshot; no model outcomes inspected", flush=True)


if __name__ == "__main__": prepare()
