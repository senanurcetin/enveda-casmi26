"""Build deterministic owned Kaggle notebook bundles without publishing or accessing data."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUNDLES = (
    ("shortcut", "build_shortcut_notebook.py", "casmi26-your-ranker-learns-how-you-built-the-pool",
     "CASMI26: Your Ranker Learns How You Built the Pool", []),
    ("pool", "build_pool_notebook.py", "casmi26-a-bigger-candidate-pool-lowered-my-score",
     "CASMI26: A Bigger Candidate Pool Lowered My Score", ["senanuretin/casmi26-structure-sources"]),
)


def main():
    for folder, builder, slug, title, datasets in BUNDLES:
        output = ROOT / "kaggle/publications" / folder
        output.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, str(ROOT / "notebooks" / builder), str(output / "notebook.ipynb")], check=True)
        metadata = {"id": "senanuretin/" + slug, "title": title, "code_file": "notebook.ipynb",
                    "language": "python", "kernel_type": "notebook", "is_private": False,
                    "enable_gpu": False, "enable_tpu": False, "enable_internet": False,
                    "dataset_sources": datasets, "kernel_sources": [],
                    "competition_sources": ["enveda-CASMI26-molecule-id-mass-spectra"]}
        (output / "kernel-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
