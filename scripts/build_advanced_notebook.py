"""Build a distinct offline fresh-training experiment; never overwrite prior evidence."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
KERNEL = "senanuretin/casmi26-fresh-training-family-holdout"
OUTPUT = ROOT / "kaggle/audits/advanced-internal"
DATASETS = ["senanuretin/casmi-coconut-b", "senanuretin/casmi-pubchemlite-b", "metric/rdkit-2026-3-3-wheel"]
SCRIPTS = ["build_structures.py", "build_store.py", "build_pool.py", "build_fingerprints.py", "run_advanced_study.py"]


def build():
    publication_path = ROOT / "reports/advanced_publication.json"
    if publication_path.exists():
        from scripts.research_snapshot import snapshot_hashes
        from scripts.check_kaggle_publications import source_hash
        publication = json.loads(publication_path.read_text(encoding="utf-8"))
        report = json.loads((ROOT / "reports/advanced_study.json").read_text(encoding="utf-8"))
        notebook = json.loads((OUTPUT / "notebook.ipynb").read_text(encoding="utf-8"))
        if (source_hash(notebook) != publication["source_sha256"]
                or snapshot_hashes(OUTPUT / "notebook.ipynb") != report["source_code_sha256"]):
            raise ValueError("Completed fresh-study source snapshot was altered")
        print("Verified frozen completed fresh-training bundle")
        return
    import nbformat
    from nbformat.v4 import new_code_cell as C, new_markdown_cell as M, new_notebook

    paths = sorted((ROOT / "src").rglob("*.py")) + [ROOT / "scripts" / name for name in SCRIPTS]
    sources = {p.relative_to(ROOT).as_posix(): p.read_text(encoding="utf-8") for p in paths}
    protocol_text = (ROOT / "protocols/advanced-internal-v1.json").read_text(encoding="utf-8")
    setup = """import hashlib, importlib.metadata, json, os, re, subprocess, sys, time
from pathlib import Path
sources = json.loads(SOURCES_LITERAL)
protocol_text = PROTOCOL_LITERAL
work = Path('/tmp/casmi_advanced_internal')
code = work / 'code'
artifacts = work / 'artifacts'
code.mkdir(parents=True, exist_ok=False)
artifacts.mkdir()
for name, source in sources.items():
    path = code / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(source, encoding='utf-8')
(code / 'protocols').mkdir()
(code / 'protocols/advanced-internal-v1.json').write_text(protocol_text,encoding='utf-8')
def unique(pattern):
    values = sorted(Path('/kaggle/input').rglob(pattern))
    if len(values) != 1:
        raise RuntimeError('Required input is missing or ambiguous: '+pattern)
    return values[0]
train = unique('train.parquet')
coconut = unique('coconut_csv_lite-09-2026.csv')
pcl = unique('PubChemLite_CCSbase_*.csv')
wheel = unique(f'rdkit-2026.3.3-cp{sys.version_info.major}{sys.version_info.minor}-*.whl')
subprocess.run([sys.executable,'-m','pip','install','--no-deps','--force-reinstall',str(wheel)],check=True)
env = {**os.environ,'CASMI_DATA':str(train.parent),'CASMI_ARTIFACTS':str(artifacts),
       'CASMI_COCONUT':str(coconut),'CASMI_PUBCHEMLITE':str(pcl),'CASMI_POOL_TAG':'_advanced_internal',
       'CASMI_QUARANTINE_MASS_CONFLICTS':'1','CASMI_QUARANTINE_INVALID_STRUCTURES':'1',
       'PYTHONPATH':str(code / 'src')+os.pathsep+str(code),'PYTHONUNBUFFERED':'1',
       'PYTHONHASHSEED':'20261025',
       'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','MKL_NUM_THREADS':'2',
       'CUBLAS_WORKSPACE_CONFIG':':4096:8','CASMI_ADVANCED_OUTPUT':'/kaggle/working/advanced_study.json'}
print('Fresh supervised training configured; no pretrained model dataset is attached.',flush=True)
""".replace("SOURCES_LITERAL", repr(json.dumps(sources, ensure_ascii=False))).replace("PROTOCOL_LITERAL", repr(protocol_text))
    execute = """def run_stage(name):
    print('Stage '+name,flush=True)
    log = work / (name+'.log')
    with log.open('w',encoding='utf-8') as output:
        process = subprocess.Popen([sys.executable,str(code/'scripts'/name)],env=env,
                                   stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        for line in process.stdout:
            output.write(line)
            output.flush()
            if re.match(r'^(Prepared |Trained fresh |Fitted prespecified |Sealed all |Completed sealed |row group [0-9])',line):
                print(line.strip(),flush=True)
        status = process.wait()
    if status:
        text = log.read_text(encoding='utf-8',errors='replace')
        errors = re.findall(r'^([A-Za-z.]*Error|AssertionError):',text,flags=re.M)
        locations = re.findall(r'File "'+re.escape(str(code))+r'/([A-Za-z0-9_/]+\\.py)", line ([0-9]+)',text)
        failure = {'status':'FAILED_CLOSED','stage':name,'exception_type':errors[-1] if errors else 'Unknown',
                   'last_project_location':list(locations[-1]) if locations else None}
        Path('/kaggle/working/advanced_failure.json').write_text(json.dumps(failure,indent=2))
        raise RuntimeError('Stage failed closed: '+name+'; '+failure['exception_type'])
for stage in ['build_structures.py','build_store.py','build_pool.py','build_fingerprints.py','run_advanced_study.py']:
    run_stage(stage)
output = Path('/kaggle/working/advanced_study.json')
report = json.loads(output.read_text())
def file_record(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):
            digest.update(block)
    return {'bytes':path.stat().st_size,'sha256':digest.hexdigest()}
report['raw_inputs'] = {'train.parquet':file_record(train),'coconut.csv':file_record(coconut),'pubchemlite.csv':file_record(pcl)}
report['source_code_sha256'] = {name:hashlib.sha256(source.encode()).hexdigest() for name,source in sources.items()}
report['execution'] = {'platform':'Kaggle','competition_source':'enveda-CASMI26-molecule-id-mass-spectra',
                       'dataset_sources':DATASET_LITERAL,'pretrained_models_attached':False,
                       'internet_enabled':False,'submitted_predictions':False,
                       'protocol_preregistered_commit':'4f8ed74','protocol_registered_before_any_new_test_results':True}
report['execution']['source_amendment_sha256'] = AMENDMENT_HASH
# Record exact installed versions without filesystem paths, URLs or environment secrets.
report['installed_packages'] = dict(sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions()
                                   if d.metadata['Name'] and re.fullmatch(r'[A-Za-z0-9_.+-]+',d.version)))
output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\\n')
print('Completed fresh-training family holdout; aggregate report sealed.',flush=True)
""".replace("DATASET_LITERAL", repr(DATASETS)).replace("AMENDMENT_HASH", repr(hashlib.sha256(
        (ROOT / "protocols/advanced-internal-v1-amendment.json").read_text(encoding="utf-8").encode()).hexdigest()))
    notebook = new_notebook(cells=[M("""<img src="https://raw.githubusercontent.com/senanurcetin/enveda-casmi26/main/site/assets/nexus-midnight.svg" width="360" alt="NEXUS research identity">

# CASMI26: Fresh Training and Chemical-Family Holdout

This predeclared internal benchmark fits every supervised weight from scratch. No archived fingerprint model or ranker is attached. Training, development and final test queries are separated by conservative generic Murcko/acyclic connectivity groups and metric identity. Two fresh fingerprint networks supply out-of-fold training features; every ranker and feature list is sealed before test feature generation.

The declared population is 600 training, 150 development and 250 test molecules. All feature ablations, blank-fragment and intensity-shuffle controls, the primary C2 contrast, joint secondary uncertainty and minimum error-slice size were fixed before test execution. Grouping applies to supervised fitting: the task-defined retrieval corpus remains available under C1/C2/C3 exclusions. This is not a newly acquired external dataset, laboratory replication or competition submission. Raw chemistry and model weights stay in scratch storage.
"""), C(setup), C(execute), M("""The report includes fresh training lineage, input/code/model hashes, partition manifests, the pre-test model seal, all favorable and adverse outcomes, paired uncertainty, candidate-coverage ceilings, descriptive source/adduct/instrument slices and exact installed package versions. Feature ablations share a fixed generated candidate universe; removing generation features is not a generator-removal experiment. No production service, billing, paid infrastructure or leaderboard submission is created.
""")], metadata={"kernelspec":{"name":"python3","display_name":"Python 3","language":"python"},"language_info":{"name":"python"}})
    for i, cell in enumerate(notebook.cells):
        cell.id = hashlib.sha256(f"{i}:{cell.cell_type}:{cell.source}".encode()).hexdigest()[:12]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, OUTPUT / "notebook.ipynb")
    meta = {"id":KERNEL,"title":"CASMI26: Fresh Training Family Holdout","code_file":"notebook.ipynb",
            "language":"python","kernel_type":"notebook","is_private":False,"enable_gpu":True,
            "enable_tpu":False,"enable_internet":False,"competition_sources":["enveda-CASMI26-molecule-id-mass-spectra"],
            "dataset_sources":DATASETS,"kernel_sources":[]}
    (OUTPUT / "kernel-metadata.json").write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print("Built distinct fresh-training family-holdout bundle")


if __name__ == "__main__":
    build()
