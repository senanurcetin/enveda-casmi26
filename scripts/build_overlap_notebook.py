"""Build the aggregate-only Kaggle overlap audit from tracked source; never publish."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell as C, new_markdown_cell as M, new_notebook

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "notebooks"))
from research_style import STYLE  # noqa: E402

OUTPUT = ROOT / "kaggle/audits/public-overlap"
KERNEL = "senanuretin/casmi26-exact-public-example-overlap-audit"
SOURCES = ["src/casmi/validate/spectrum_overlap.py", "src/casmi/validate/overlap_audit.py",
           "scripts/audit_public_overlap.py"]


def build():
    sources = {name: (ROOT / name).read_text(encoding="utf-8") for name in SOURCES}
    setup = """# Source files are bundled from the tracked repository; no internet/install step.
import json, sys
from pathlib import Path
bundle = Path('/kaggle/working/casmi_overlap_source')
sources = json.loads(SOURCE_JSON)
for name, text in sources.items():
    path = bundle / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8')
sys.path.insert(0, str(bundle / 'src'))
sys.path.insert(0, str(bundle / 'scripts'))
from audit_public_overlap import run_audit
print('Loaded content-audit source; no model or prediction code is invoked.')
""".replace("SOURCE_JSON", repr(json.dumps(sources, ensure_ascii=False)))
    run = """# There must be exactly one attached competition train/test pair.
pairs = [path for path in Path('/kaggle/input').rglob('train.parquet')
         if path.with_name('test.parquet').is_file()]
if len(pairs) != 1:
    raise RuntimeError('Expected exactly one attached competition train/test pair')
result = run_audit(pairs[0], pairs[0].with_name('test.parquet'), batch_size=8192,
                   progress=lambda count, total: print(f'Scanned {count:,}/{total:,} train rows', flush=True))
result['execution'] = {
    'platform': 'Kaggle',
    'competition_source': 'enveda-CASMI26-molecule-id-mass-spectra',
    'input_scope': 'mounted public train.parquet and test.parquet; not hidden evaluation',
}
with Path('/kaggle/working/public_overlap_audit.json').open('x', encoding='utf-8') as output:
    json.dump(result, output, indent=2, allow_nan=False)
    output.write('\\n')
print(json.dumps(result, indent=2, allow_nan=False))
"""
    plot = """labels = ['Public examples', 'Exact numeric content', '+ Adduct/mode agreement']
counts = [result['public_spectra'], result['content_matched_spectra'],
          result['metadata_consistent_matched_spectra']]
fig, ax = plt.subplots(figsize=(10, 4))
bars = ax.barh(labels[::-1], counts[::-1], color=['#8F9BFF', '#00E0FF', '#4EA8FF'])
ax.bar_label(bars, padding=5)
ax.set_xlim(0, max(counts) * 1.15)
ax.set_xlabel('Public-example spectra (unique target rows)')
ax.set_title('Executed exact-content overlap audit — not a model evaluation')
plt.tight_layout()
plt.show()
"""
    cells = [M("""<img src="https://raw.githubusercontent.com/senanurcetin/enveda-casmi26/main/site/assets/nexus-midnight.svg" alt="NEXUS official identity" width="360">

# CASMI26: Exact Public-Example Overlap Audit

**NEXUS / MOLECULAR IDENTIFICATION RESEARCH** · Senanur Cetin · Noncommercial research

The historical profile found coarse-key matches between public examples and training spectra. A rounded precursor, peak count and first peak cannot establish exact overlap. This notebook reruns that question on the attached public files, comparing **every stored fragment m/z and normalized intensity, in order, plus the unrounded precursor value**.

It separately checks adduct and ionization-mode agreement. Equality of stored numeric content does not independently authenticate acquisition identity. The mounted public example file is distinct from the competition's hidden evaluation set.

Only aggregate counts, input file SHA-256 hashes, source-code hashes and software versions are exported. Raw molecule/spectrum identifiers, SMILES, peak arrays, models and predictions are not published. No model training, performance scoring or competition submission occurs. CPU only; internet disabled.
"""), C(STYLE), M("## 1. Load the versioned audit code"), C(setup),
             M("## 2. Scan all training rows and verify full content\n\nArrow batches limit memory. Exact precursor equality is a necessary prefilter, followed by full-array comparison. The scanner fails if the visited row count does not match the Parquet footer. Input hashes describe the exact mounted files used by this run."),
             C(run), M("## 3. Read the result without overstating it"), C(plot),
             M("""A matched public example is unsuitable as an independent unseen-spectrum holdout. This finding says nothing by itself about hidden-test contamination, model ranking quality or prospective chemical identification.

Source-library match counts are **not additive**: one public spectrum can match multiple training records/libraries. A molecule is fully matched only if all its public spectra have matches. A matching content pair with conflicting acquisition metadata is counted separately.

An exact comparison is deliberately stricter than tolerance-based spectral similarity. A stored floating-point difference or a changed peak order produces a non-match even if the underlying acquisition was related. Inspect the JSON for exact/metadata distinctions and denominator coverage; do not substitute rounded-key matches for exact matches.

[Repository methodology and model card](https://github.com/senanurcetin/enveda-casmi26/blob/main/MODEL_CARD.md). Historical model scores remain unchanged by this data audit.
""")]
    nb = new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                            "language_info": {"name": "python"}})
    for index, cell in enumerate(nb.cells):
        cell.id = hashlib.sha256(f"{index}:{cell.cell_type}:{cell.source}".encode()).hexdigest()[:12]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, OUTPUT / "notebook.ipynb")
    meta = {"id": KERNEL, "title": "CASMI26: Exact Public-Example Overlap Audit", "code_file": "notebook.ipynb",
            "language": "python", "kernel_type": "notebook", "is_private": False,
            "enable_gpu": False, "enable_tpu": False, "enable_internet": False,
            "competition_sources": ["enveda-CASMI26-molecule-id-mass-spectra"],
            "dataset_sources": [], "kernel_sources": []}
    (OUTPUT / "kernel-metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print("Built aggregate-only overlap audit bundle")


if __name__ == "__main__":
    build()
