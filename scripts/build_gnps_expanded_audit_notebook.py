"""Broad community source feasibility under the unchanged external novelty gate."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.build_external_audit_notebook import build_configured

KERNEL = "senanuretin/casmi26-gnps-community-external-independence"
DATASETS = ["senanuretin/casmi26-gnps-cc0-expanded", "metric/rdkit-2026-3-3-wheel"]
OUTPUT = ROOT / "kaggle/audits/gnps-expanded-independence"
REGISTERED_COMMIT = "6be58a5"


def build():
    publication_path = ROOT / "reports/gnps_expanded_external_audit_publication.json"
    if publication_path.exists():
        from scripts.research_snapshot import snapshot_hashes
        from scripts.check_kaggle_publications import source_hash
        publication = json.loads(publication_path.read_text(encoding="utf-8"))
        report = json.loads((ROOT / "reports/gnps_expanded_external_audit.json").read_text(encoding="utf-8"))
        notebook = json.loads((OUTPUT / "notebook.ipynb").read_text(encoding="utf-8"))
        if source_hash(notebook) != publication["source_sha256"] or snapshot_hashes(OUTPUT / "notebook.ipynb") != report["source_code_sha256"]:
            raise ValueError("Completed expanded GNPS source snapshot was altered")
        print("Verified frozen completed expanded GNPS independence bundle")
        return
    build_configured(OUTPUT, KERNEL, DATASETS, ROOT / "protocols/external-source-v3.json",
                     "gnps_cc0.parquet", REGISTERED_COMMIT, "GNPS Community")


if __name__ == "__main__": build()
