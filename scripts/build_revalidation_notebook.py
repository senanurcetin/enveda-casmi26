"""Build a scratch-only, offline Kaggle full-source revalidation notebook."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUTPUT = ROOT / "kaggle/audits/corrected-revalidation"
KERNEL = "senanuretin/casmi26-corrected-models-and-negative-controls"
SCRIPTS = ["build_structures.py", "build_store.py", "build_pool.py", "build_fingerprints.py",
           "run_revalidation_study.py", "run_validation.py", "make_submission.py", "compare_control_runs.py"]
DATASETS = ["senanuretin/casmi-coconut-b", "senanuretin/casmi-pubchemlite-b",
            "senanuretin/casmi-models-b", "metric/rdkit-2026-3-3-wheel"]


def build():
    # An executed public bundle is a historical experiment, not a moving copy
    # of main. New experiments get a distinct notebook and source snapshot.
    if (ROOT / "reports/revalidation_publication.json").exists():
        from scripts.research_snapshot import snapshot_hashes
        from scripts.check_kaggle_publications import source_hash
        publication = json.loads((ROOT / "reports/revalidation_publication.json").read_text(encoding="utf-8"))
        result = json.loads((ROOT / "reports/revalidation_study.json").read_text(encoding="utf-8"))
        notebook = json.loads((OUTPUT / "notebook.ipynb").read_text(encoding="utf-8"))
        if (source_hash(notebook) != publication["source_sha256"]
                or snapshot_hashes(OUTPUT / "notebook.ipynb") != result["source_code_sha256"]):
            raise ValueError("Executed historical source snapshot was altered")
        print("Verified frozen executed revalidation bundle; new work requires a distinct study")
        return
    import nbformat
    from nbformat.v4 import new_code_cell as C, new_markdown_cell as M, new_notebook

    paths = sorted((ROOT / "src").rglob("*.py")) + [ROOT / "scripts" / name for name in SCRIPTS]
    sources = {path.relative_to(ROOT).as_posix(): path.read_text(encoding="utf-8") for path in paths}
    setup = """import hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path
sources = json.loads(SOURCE_JSON)
work = Path('/tmp/casmi_corrected_study')
code = work / 'code'
artifacts = work / 'artifacts'
code.mkdir(parents=True, exist_ok=False)
for name, source in sources.items():
    path = code / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(source, encoding='utf-8')

def unique(pattern):
    found = sorted(Path('/kaggle/input').rglob(pattern))
    if len(found) != 1:
        raise RuntimeError('Required input is missing or ambiguous: '+pattern)
    return found[0]

train = unique('train.parquet')
coconut = unique('coconut_csv_lite-09-2026.csv')
pcl = unique('PubChemLite_CCSbase_*.csv')
models = unique('fp_fold0.pt').parent
wheel = unique(f'rdkit-2026.3.3-cp{sys.version_info.major}{sys.version_info.minor}-*.whl')
subprocess.run([sys.executable, '-m', 'pip', 'install', '--no-deps', '--force-reinstall', str(wheel)], check=True)
artifacts.mkdir()
(artifacts / 'models').mkdir()
names = ['fp_fold0.pt','fp_fold1.pt'] + [f'ranker-v16-gen2_seed{seed}.txt' for seed in range(5)]
for name in names:
    shutil.copyfile(models / name, artifacts / 'models' / name)
env = {**os.environ, 'CASMI_ARTIFACTS':str(artifacts), 'CASMI_DATA':str(train.parent),
       'CASMI_COCONUT':str(coconut), 'CASMI_PUBCHEMLITE':str(pcl), 'CASMI_POOL_TAG':'_corrected_study',
       'CASMI_QUARANTINE_MASS_CONFLICTS':'1',
       'CASMI_QUARANTINE_INVALID_STRUCTURES':'1',
       'PYTHONPATH':str(code / 'src')+os.pathsep+str(code), 'PYTHONUNBUFFERED':'1',
       'OMP_NUM_THREADS':'2', 'OPENBLAS_NUM_THREADS':'2', 'MKL_NUM_THREADS':'2',
       'CUBLAS_WORKSPACE_CONFIG':':4096:8', 'CASMI_STUDY_OUTPUT':'/kaggle/working/revalidation_study.json'}
print('Fresh source rebuild configured; derived spectra/models stay under /tmp, not notebook outputs.')
""".replace("SOURCE_JSON", repr(json.dumps(sources, ensure_ascii=False)))
    execute = """import re
def run_stage(name):
    print('Starting '+name, flush=True)
    log = work / (name + '.log')
    with log.open('w', encoding='utf-8') as output:
        process = subprocess.Popen([sys.executable, str(code / 'scripts' / name)], env=env,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in process.stdout:
            output.write(line)
            output.flush()
            # Publish only stage counts/timing, never source rows or chemistry errors.
            if re.match(r'^(row group [0-9]|Prepared [0-9]|Completed [0-9]|[0-9,]+ fingerprints)', line):
                print(line.strip(), flush=True)
        result = process.wait()
    if result:
        text = log.read_text(encoding='utf-8', errors='replace')
        exceptions = re.findall(r'^(?:[A-Za-z.]*Error|AssertionError):.*$', text, flags=re.M)
        message = exceptions[-1] if exceptions else 'Stage exited without a completed aggregate'
        # These project guards produce aggregate-only messages. Other errors get only their type.
        safe = message if re.match(r'^(ValueError: (Conflicting neutral mass for [0-9]+ shared InChIKey14 groups|Malformed candidate source InChIKey14|Fingerprint provenance mismatch|Invalid generated candidate identity at row [0-9]+)|FileExistsError: Refusing)', message) else message.split(':',1)[0]
        Path('/kaggle/working/revalidation_failure.json').write_text(json.dumps({'status':'FAILED_CLOSED','stage':name,'safe_reason':safe},indent=2))
        raise RuntimeError('Revalidation stopped at '+name+'; '+safe)
for step in ['build_structures.py','build_store.py','build_pool.py','build_fingerprints.py','run_revalidation_study.py']:
    run_stage(step)
report_path = Path('/kaggle/working/revalidation_study.json')
report = json.loads(report_path.read_text())
def file_record(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):
            digest.update(block)
    return {'bytes':path.stat().st_size,'sha256':digest.hexdigest()}
report['raw_inputs'] = {'train.parquet':file_record(train),'coconut.csv':file_record(coconut),'pubchemlite.csv':file_record(pcl)}
report['source_code_sha256'] = {name:hashlib.sha256(source.encode()).hexdigest() for name,source in sources.items()}
report['execution'] = {'platform':'Kaggle','competition_source':'enveda-CASMI26-molecule-id-mass-spectra',
                       'dataset_sources':DATASET_JSON,'internet_enabled':False,'submitted_predictions':False}
report_path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\\n')
print('Completed corrected-protocol retrospective revalidation and negative controls.')
for model, block in report['models'].items():
    print(model, {mode:{cls:round(values['mrr'],5) for cls,values in classes.items()} for mode,classes in block['conditions'].items()})
""".replace("DATASET_JSON", repr(DATASETS))
    nb = new_notebook(cells=[M("""<img src="https://raw.githubusercontent.com/senanurcetin/enveda-casmi26/main/site/assets/nexus-midnight.svg" width="360" alt="NEXUS research identity">

# CASMI26: Corrected Models and Negative Controls

This is a retrospective study of the designated V-in natural-product population using rebuilt full-source artifacts and frozen original weights. D1, D1+D2 and the original five-seed v16 learned ranker are evaluated under corrected source-library/tautomer exclusions. The same molecules also receive blank-fragment and deterministic intensity-shuffle treatments.

No public-example test file is evaluated. No new model is trained and no predictions are submitted. Historical model selection and incomplete original weight-training lineage prevent interpreting this as a new untouched generalization benchmark. Reports contain aggregates and hashes only; source arrays, structures, identifiers, per-molecule ranks and trained weights remain in scratch storage.
"""), C(setup), C(execute), M("""The exported report separates absolute retrospective summaries, molecule-paired control diagnostics and model contrasts. Confidence intervals describe this selected population and the fixed novelty-mixture grid; they do not correct for historical model selection or instrument shift. Blank controls retain precursor, adduct, collision energy and candidate-mass evidence. Intensity shuffling retains fragment locations. No ACCEPT/REJECT promotion decision or new competition score is justified by this study alone.
""")], metadata={"kernelspec":{"name":"python3","display_name":"Python 3","language":"python"},"language_info":{"name":"python"}})
    for index, cell in enumerate(nb.cells):
        cell.id = hashlib.sha256(f"{index}:{cell.cell_type}:{cell.source}".encode()).hexdigest()[:12]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, OUTPUT / "notebook.ipynb")
    meta = {"id":KERNEL,"title":"CASMI26: Corrected Models and Negative Controls","code_file":"notebook.ipynb",
            "language":"python","kernel_type":"notebook","is_private":False,"enable_gpu":True,"enable_tpu":False,
            "enable_internet":False,"competition_sources":["enveda-CASMI26-molecule-id-mass-spectra"],
            "dataset_sources":DATASETS,"kernel_sources":[]}
    (OUTPUT / "kernel-metadata.json").write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print("Built scratch-only corrected revalidation notebook")


if __name__ == "__main__":
    build()
