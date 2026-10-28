"""Prepare the separately registered family-purged external study on free CPU."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUTPUT = ROOT / "kaggle/audits/purged-external"
KERNEL = "senanuretin/casmi26-purged-family-external-validation"
DATASETS = ["senanuretin/casmi-coconut-b", "senanuretin/casmi-pubchemlite-b", "metric/rdkit-2026-3-3-wheel",
            "senanuretin/casmi26-gnps-cc0-expanded"]
REGISTERED_COMMIT = "1b7a6cc"


def build():
    publication_path = ROOT / "reports/purged_external_publication.json"
    if publication_path.exists():
        from scripts.research_snapshot import snapshot_hashes
        from scripts.check_kaggle_publications import source_hash
        publication = json.loads(publication_path.read_text(encoding="utf-8"))
        report = json.loads((ROOT / "reports/purged_external_study.json").read_text(encoding="utf-8"))
        notebook = json.loads((OUTPUT / "notebook.ipynb").read_text(encoding="utf-8"))
        if source_hash(notebook) != publication["source_sha256"] or snapshot_hashes(OUTPUT / "notebook.ipynb") != report["source_code_sha256"]:
            raise ValueError("Completed purged external study source was altered")
        print("Verified frozen completed family-purged external study")
        return
    original = json.loads((ROOT / "reports/gnps_expanded_external_audit.json").read_text(encoding="utf-8"))
    publication = json.loads((ROOT / "reports/gnps_expanded_external_audit_publication.json").read_text(encoding="utf-8"))
    if original["status"] != "EXTERNAL_INDEPENDENCE_NOT_ESTABLISHED" or publication["status"] != "pass":
        raise ValueError("Preserve authenticated unsuccessful original novelty feasibility first")
    import nbformat
    from nbformat.v4 import new_notebook, new_markdown_cell as M, new_code_cell as C
    script_names = ["build_structures.py", "build_store.py", "build_pool.py", "build_fingerprints.py", "run_external_audit.py",
        "run_advanced_study.py", "external_study_training.py", "run_external_c3_study.py", "prepare_purged_external.py"]
    paths = sorted((ROOT / "src").rglob("*.py")) + [ROOT / "scripts" / name for name in script_names]
    paths += [ROOT / "protocols" / name for name in ("purged-external-v1.json", "purged-external-v1-execution.json",
        "external-source-v3.json", "advanced-internal-v1.json")]
    paths += [ROOT / "reports/advanced_study.json", ROOT / "reports/gnps_expanded_external_audit.json"]
    sources = {p.relative_to(ROOT).as_posix(): p.read_text(encoding="utf-8") for p in paths}
    setup = '''import hashlib,importlib.metadata,json,os,re,subprocess,sys,time
from pathlib import Path
sources=json.loads(SOURCES_LITERAL)
work=Path('/tmp/casmi_purged_external');code=work/'code';original_artifacts=work/'original-artifacts'
code.mkdir(parents=True,exist_ok=False);original_artifacts.mkdir()
for name,text in sources.items():
    path=code/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
def unique(pattern):
    paths=sorted(Path('/kaggle/input').rglob(pattern))
    if len(paths)!=1:raise ValueError('Required declared input missing or ambiguous')
    return paths[0]
train=unique('train.parquet');coconut=unique('coconut_csv_lite-09-2026.csv');pcl=unique('PubChemLite_CCSbase_*.csv')
external=unique('gnps_cc0.parquet');external_manifest=unique('external_input_manifest.json')
wheel=unique(f'rdkit-2026.3.3-cp{sys.version_info.major}{sys.version_info.minor}-*.whl')
subprocess.run([sys.executable,'-m','pip','install','--no-deps','--force-reinstall',str(wheel)],check=True)
purged_data=work/'purged-data';artifacts=work/'filtered-artifacts'
env={**os.environ,'CASMI_DATA':str(train.parent),'CASMI_ARTIFACTS':str(original_artifacts),
    'CASMI_COCONUT':str(coconut),'CASMI_PUBCHEMLITE':str(pcl),'CASMI_POOL_TAG':'_purge_original',
    'CASMI_QUARANTINE_MASS_CONFLICTS':'1','CASMI_QUARANTINE_INVALID_STRUCTURES':'1',
    'CASMI_ORIGINAL_TRAIN':str(train),'CASMI_ORIGINAL_EXTERNAL':str(external),
    'CASMI_ORIGINAL_EXTERNAL_MANIFEST':str(external_manifest),'CASMI_PURGED_DATA':str(purged_data),
    'PYTHONPATH':str(code/'src')+os.pathsep+str(code),'PYTHONUNBUFFERED':'1','PYTHONHASHSEED':'20261028',
    'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','MKL_NUM_THREADS':'2','NUMBA_NUM_THREADS':'2'}
print('Original-snapshot novelty feasibility failed and remains preserved. Separate purged-training experiment begins.',flush=True)
'''.replace("SOURCES_LITERAL", repr(json.dumps(sources)))
    execute = '''def run_stage(name,label=None):
    print('Stage '+(label or name),flush=True)
    log=work/((label or name)+'.log')
    with log.open('w',encoding='utf-8') as output:
        process=subprocess.Popen([sys.executable,str(code/'scripts'/name)],env=env,stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,text=True,bufsize=1)
        for line in process.stdout:
            output.write(line);output.flush()
            if re.match(r'^(Prepared |Trained fresh |Fitted prespecified |Sealed all |Completed sealed |External audit scanned |row group [0-9])',line):
                print(line.strip(),flush=True)
        status=process.wait()
    if status:
        text=log.read_text(encoding='utf-8',errors='replace')
        errors=re.findall(r'^([A-Za-z.]*Error|AssertionError):',text,flags=re.M)
        locations=re.findall(r'File "'+re.escape(str(code))+r'/([A-Za-z0-9_/]+\\.py)", line ([0-9]+)',text)
        failure={'status':'FAILED_CLOSED','stage':label or name,'exception_type':errors[-1] if errors else 'Unknown',
                 'last_project_location':list(locations[-1]) if locations else None}
        Path('/kaggle/working/purged_external_failure.json').write_text(json.dumps(failure,indent=2))
        raise RuntimeError('Purged-training stage failed closed; sanitized failure evidence exported')
for name in ['build_structures.py','build_store.py']:
    run_stage(name,'original-'+name)
run_stage('prepare_purged_external.py')
artifacts.mkdir()
env.update(CASMI_DATA=str(purged_data),CASMI_ARTIFACTS=str(artifacts),CASMI_POOL_TAG='_purged_external',
    CASMI_EXTERNAL_INPUT=str(purged_data/'selected_external.parquet'),
    CASMI_EXTERNAL_MANIFEST=str(purged_data/'selected_external_manifest.json'),
    CASMI_PURGE_REPORT=str(purged_data/'purge_report.json'),
    CASMI_PROTECTED_POPULATION=str(purged_data/'protected_population.json'),
    CASMI_STUDY_PROTOCOL='protocols/purged-external-v1.json',
    CASMI_EXECUTION_PROTOCOL='protocols/purged-external-v1-execution.json',
    CASMI_EXTERNAL_PROTOCOL='protocols/external-source-v3.json',
    CASMI_EXTERNAL_STUDY_OUTPUT='/kaggle/working/purged_external_study.json')
for name in ['build_structures.py','build_store.py','build_pool.py','build_fingerprints.py','run_external_c3_study.py']:
    run_stage(name,'filtered-'+name)
output=Path('/kaggle/working/purged_external_study.json');report=json.loads(output.read_text())
def file_record(path):
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest}
report['raw_inputs']={'original_train.parquet':file_record(train),'purged_train.parquet':file_record(purged_data/'train.parquet'),
    'coconut.csv':file_record(coconut),'pubchemlite.csv':file_record(pcl),'gnps_cc0.parquet':file_record(external)}
report['source_code_sha256']={n:hashlib.sha256(t.encode()).hexdigest() for n,t in sources.items()}
report['execution']={'platform':'Kaggle','competition_source':'enveda-CASMI26-molecule-id-mass-spectra',
    'dataset_sources':DATASET_LITERAL,'pretrained_models_attached':False,'gpu_enabled':False,'internet_enabled':False,
    'submitted_predictions':False,'protocol_registered_commit':REGISTERED_COMMIT,'new_model_test_results_before_registration':False}
report['installed_packages']=dict(sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions()
    if d.metadata['Name'] and re.fullmatch(r'[A-Za-z0-9_.+-]+',d.version)))
output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\\n')
print('Completed separately registered purged-training external-family study; original novelty failures remain explicit.',flush=True)
'''.replace("DATASET_LITERAL", repr(DATASETS)).replace("REGISTERED_COMMIT", repr(REGISTERED_COMMIT))
    nb = new_notebook(cells=[M('''<img src="https://raw.githubusercontent.com/senanurcetin/enveda-casmi26/main/site/assets/nexus-midnight.svg" width="360" alt="NEXUS">

# CASMI26: Purged Family External Validation

The original-snapshot novelty study did not find the required independent population. This is a separately registered experiment: choose external generic families before predictions, remove their spectra and molecular identities from the entire training/reference corpus, rebuild every artifact, fit every supervised weight from scratch, then seal both candidate budgets and the development choice before test features. Historical test families remain protected. Generalization is relative to this purged training snapshot. No original-snapshot novelty, independent laboratory replication, paid compute, service deployment or competition submission is claimed.
'''), C(setup), C(execute)], metadata={"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"}, "language_info": {"name": "python"}})
    for index, cell in enumerate(nb.cells): cell.id = hashlib.sha256(f"{index}:{cell.source}".encode()).hexdigest()[:12]
    OUTPUT.mkdir(parents=True, exist_ok=True); nbformat.write(nb, OUTPUT / "notebook.ipynb")
    metadata = {"id": KERNEL, "title": "CASMI26: Purged Family External Validation", "code_file": "notebook.ipynb",
        "language": "python", "kernel_type": "notebook", "is_private": False, "enable_gpu": False, "enable_tpu": False,
        "enable_internet": False, "competition_sources": ["enveda-CASMI26-molecule-id-mass-spectra"],
        "dataset_sources": DATASETS, "kernel_sources": []}
    (OUTPUT / "kernel-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print("Built distinct CPU-only family-purged external study")


if __name__ == "__main__": build()
