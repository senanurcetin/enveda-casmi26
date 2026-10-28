"""Authenticate expanded GNPS feasibility without conflating publication and science."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.build_gnps_expanded_audit_notebook import KERNEL, DATASETS, REGISTERED_COMMIT
from scripts.check_external_audit_publication import authenticate, validate_result

OPTIONS = dict(protocol_path="protocols/external-source-v3.json", prepared_path="reports/gnps_expanded_source_preparation.json",
               notebook_path="kaggle/audits/gnps-expanded-independence/notebook.ipynb", external_filename="gnps_cc0.parquet",
               registered_commit=REGISTERED_COMMIT)


if __name__ == "__main__":
    authenticate(kernel=KERNEL, datasets=DATASETS, output_prefix="gnps_expanded_external_audit", **OPTIONS)
