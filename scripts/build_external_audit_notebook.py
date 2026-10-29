"""Build an offline CPU source-independence audit before any new model selection."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "kaggle/audits/external-independence"
KERNEL = "senanuretin/casmi26-massbank-external-independence"
DATASETS = ["senanuretin/casmi26-massbank-cc0-202603", "metric/rdkit-2026-3-3-wheel"]


def build():
    publication_path = ROOT / "reports/external_audit_publication.json"
    if publication_path.exists():
        import sys
        sys.path.insert(0, str(ROOT))
        from scripts.research_snapshot import snapshot_hashes
        from scripts.check_kaggle_publications import source_hash
        publication = json.loads(publication_path.read_text(encoding="utf-8"))
        report = json.loads((ROOT / "reports/external_audit.json").read_text(encoding="utf-8"))
        notebook = json.loads((OUTPUT / "notebook.ipynb").read_text(encoding="utf-8"))
        if (source_hash(notebook) != publication["source_sha256"]
                or snapshot_hashes(OUTPUT / "notebook.ipynb") != report["source_code_sha256"]):
            raise ValueError("Completed external audit source snapshot was altered")
        print("Verified frozen completed external-independence bundle")
        return
    build_configured(OUTPUT, KERNEL, DATASETS, ROOT / "protocols/external-c3-v1.json", "massbank_cc0.parquet", "6164c5e", "MassBank")


def build_configured(output_dir, kernel, datasets, protocol_path, external_filename, registered_commit, source_label):
    import nbformat
    from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell
    paths = sorted((ROOT / "src").rglob("*.py")) + [ROOT / "scripts/run_external_audit.py"]
    sources = {p.relative_to(ROOT).as_posix(): p.read_text(encoding="utf-8") for p in paths}
    protocol = protocol_path.read_text(encoding="utf-8")
    code = '''import hashlib,json,os,subprocess,sys,re
from pathlib import Path
sources = json.loads(SOURCE_LITERAL)
protocol_text = PROTOCOL_LITERAL
work=Path('/tmp/casmi_external_audit')
work.mkdir(exist_ok=False)
for name,text in sources.items():
    path=work/name
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8')
protocol=work/'protocol.json'
protocol.write_text(protocol_text,encoding='utf-8')
def unique(pattern):
    paths=sorted(Path('/kaggle/input').rglob(pattern))
    if len(paths)!=1:raise ValueError('Required input missing or ambiguous')
    return paths[0]
wheel=unique(f'rdkit-2026.3.3-cp{sys.version_info.major}{sys.version_info.minor}-*.whl')
subprocess.run([sys.executable,'-m','pip','install','--no-deps','--force-reinstall',str(wheel)],check=True)
output=Path('/kaggle/working/external_audit.json')
command=[sys.executable,str(work/'scripts/run_external_audit.py'),'--train',str(unique('train.parquet')),
         '--external',str(unique(EXTERNAL_FILENAME)),'--manifest',str(unique('external_input_manifest.json')),
         '--protocol',str(protocol),'--output',str(output)]
env={**os.environ,'PYTHONPATH':str(work/'src')+os.pathsep+str(work),'PYTHONUNBUFFERED':'1'}
with (work/'audit-private.log').open('w',encoding='utf-8') as log:
    process=subprocess.Popen(command,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    for line in process.stdout:
        log.write(line)
        if line.startswith(('External audit scanned','External metadata/content')):print(line.strip(),flush=True)
    status=process.wait()
if status:
    text=(work/'audit-private.log').read_text(encoding='utf-8',errors='replace')
    errors=re.findall(r'^([A-Za-z.]*Error|AssertionError):',text,flags=re.M)
    Path('/kaggle/working/external_audit_failure.json').write_text(json.dumps({'status':'FAILED_CLOSED','exception_type':errors[-1] if errors else 'Unknown'}))
    raise RuntimeError('External source audit failed closed')
report=json.loads(output.read_text())
report['source_code_sha256']={n:hashlib.sha256(t.encode()).hexdigest() for n,t in sources.items()}
report['execution']={'platform':'Kaggle','internet_enabled':False,'gpu_enabled':False,'model_training':False,'submitted_predictions':False,'protocol_registered_commit':REGISTERED_COMMIT}
output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\\n')
print('Aggregate external feasibility evidence exported; no raw source records or model outcomes exported.')
'''.replace("SOURCE_LITERAL", repr(json.dumps(sources))).replace("PROTOCOL_LITERAL", repr(protocol)).replace("EXTERNAL_FILENAME", repr(external_filename)).replace("REGISTERED_COMMIT", repr(registered_commit))
    nb = new_notebook(cells=[new_markdown_cell('''<img src="https://raw.githubusercontent.com/senanurcetin/enveda-casmi26/main/site/assets/nexus-midnight.svg" width="360" alt="NEXUS">

# CASMI26: SOURCE_LABEL External Independence

Feasibility audit before model selection: pinned CC0 source records are compared with every CASMI training row for standard/metric identity and exact stored spectral content. The embedded protocol records the source, license basis and eligibility policy. This stage reports overlap and eligibility counts only. It trains no model, inspects no prediction scores and does not establish laboratory replication. Only aggregate counts and provenance hashes are exported.
'''.replace("SOURCE_LABEL", source_label)), new_code_cell(code)], metadata={"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                                    "language_info": {"name": "python"}})
    for i, cell in enumerate(nb.cells):
        cell.id = hashlib.sha256(f"{i}:{cell.source}".encode()).hexdigest()[:12]
    output_dir.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, output_dir / "notebook.ipynb")
    metadata = {"id": kernel, "title": f"CASMI26: {source_label} External Independence", "code_file": "notebook.ipynb",
                "language": "python", "kernel_type": "notebook", "is_private": False, "enable_gpu": False,
                "enable_tpu": False, "enable_internet": False, "competition_sources": ["enveda-CASMI26-molecule-id-mass-spectra"],
                "dataset_sources": datasets, "kernel_sources": []}
    (output_dir / "kernel-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print("Built aggregate-only external-independence audit")


if __name__ == "__main__":
    build()
