"""Build a distinct CPU-only unopened-family experiment; never overwrite an attempt."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.build_purged_external_notebook import DATASETS
from scripts.check_kaggle_publications import source_hash
OUTPUT=ROOT/'kaggle/audits/unopened-coverage'
KERNEL='senanuretin/casmi26-unopened-family-coverage'
REGISTERED_COMMIT='3703ba1'

def build():
 history=ROOT/'reports/unopened_coverage_execution_history.json'
 if history.exists():
  record=json.loads(history.read_text())['attempts'][-1]
  if source_hash(json.loads((OUTPUT/'notebook.ipynb').read_text()))!=record['source_sha256']:raise ValueError('Registered attempt source changed')
  print('Verified registered unopened study; no regeneration');return
 import nbformat
 scripts=['build_structures.py','build_store.py','build_pool.py','build_fingerprints.py','run_advanced_study.py','run_unopened_coverage_study.py']
 paths=sorted((ROOT/'src').rglob('*.py'))+[ROOT/'scripts'/name for name in scripts]
 paths+=sorted((ROOT/'protocols').glob('*.json'))
 paths+=[ROOT/'reports'/name for name in ['advanced_study.json','purged_external_study.json','gnps_expanded_source_preparation.json']]
 sources={p.relative_to(ROOT).as_posix():p.read_text(encoding='utf-8') for p in paths}
 setup="""import hashlib,importlib.metadata,json,os,re,subprocess,sys
from pathlib import Path
sources=json.loads(SOURCE_LITERAL)
work=Path('/tmp/casmi_unopened_coverage');code=work/'code';artifacts=work/'artifacts'
code.mkdir(parents=True,exist_ok=False);artifacts.mkdir()
for name,text in sources.items():
 path=code/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
def unique(pattern):
 paths=sorted(Path('/kaggle/input').rglob(pattern))
 if len(paths)!=1:raise ValueError('Missing or ambiguous declared input')
 return paths[0]
train=unique('train.parquet');coconut=unique('coconut_csv_lite-09-2026.csv');pcl=unique('PubChemLite_CCSbase_*.csv');external=unique('gnps_cc0.parquet')
wheel=unique(f'rdkit-2026.3.3-cp{sys.version_info.major}{sys.version_info.minor}-*.whl')
subprocess.run([sys.executable,'-m','pip','install','--no-deps','--force-reinstall',str(wheel)],check=True)
def file_record(path):
 with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
 return {'bytes':path.stat().st_size,'sha256':digest}
raw={'train.parquet':file_record(train),'coconut.csv':file_record(coconut),'pubchemlite.csv':file_record(pcl),'gnps_cc0.parquet':file_record(external)}
old=json.loads((code/'reports/advanced_study.json').read_text())
for name in ('train.parquet','coconut.csv','pubchemlite.csv'):
 if raw[name]!=old['raw_inputs'][name]:raise ValueError('Declared source bytes differ from original protected input')
env={**os.environ,'CASMI_DATA':str(train.parent),'CASMI_ARTIFACTS':str(artifacts),'CASMI_POOL_TAG':'','CASMI_COCONUT':str(coconut),'CASMI_PUBCHEMLITE':str(pcl),'CASMI_HISTORY_EXTERNAL':str(external),'CASMI_QUARANTINE_MASS_CONFLICTS':'1','CASMI_QUARANTINE_INVALID_STRUCTURES':'1','CASMI_UNOPENED_OUTPUT':'/kaggle/working/unopened_coverage_study.json','PYTHONPATH':str(code/'src')+os.pathsep+str(code),'PYTHONUNBUFFERED':'1','PYTHONHASHSEED':'20261030','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','MKL_NUM_THREADS':'2','NUMBA_NUM_THREADS':'2'}
""".replace('SOURCE_LITERAL',repr(json.dumps(sources)))
 execute="""for name in ['build_structures.py','build_store.py','build_pool.py','build_fingerprints.py','run_unopened_coverage_study.py']:
 print('Stage '+name,flush=True)
 log=work/(name+'.log')
 with log.open('w',encoding='utf-8') as out:
  proc=subprocess.Popen([sys.executable,str(code/'scripts'/name)],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  for line in proc.stdout:
   out.write(line);out.flush()
   if re.match(r'^(Prepared |Sealed |Completed |row group [0-9])',line):print(line.strip(),flush=True)
  status=proc.wait()
 if status:
  text=log.read_text(encoding='utf-8',errors='replace')
  errors=re.findall(r'^([A-Za-z.]*Error|AssertionError):',text,flags=re.M)
  locations=re.findall(r'File "'+re.escape(str(code))+r'/([A-Za-z0-9_/]+\.py)", line ([0-9]+)',text)
  failure={'status':'FAILED_CLOSED','stage':name,'exception_type':errors[-1] if errors else 'Unknown','last_project_location':list(locations[-1]) if locations else None}
  Path('/kaggle/working/unopened_coverage_failure.json').write_text(json.dumps(failure,indent=2))
  raise RuntimeError('Unopened-family study failed closed; sanitized evidence exported')
output=Path(env['CASMI_UNOPENED_OUTPUT']);result=json.loads(output.read_text())
result['raw_inputs']=raw
result['source_code_sha256']={name:hashlib.sha256(text.encode()).hexdigest() for name,text in sources.items()}
result['execution']={'registered_commit':REGISTERED,'dataset_sources':DATASETS,'competition_source':'enveda-CASMI26-molecule-id-mass-spectra','gpu_enabled':False,'internet_enabled':False,'supervised_models_fitted':0,'submitted_predictions':False}
result['installed_packages']={name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','pyarrow','rdkit','numba']}
output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\\n')
print('Completed aggregate-only unopened study',flush=True)
""".replace('REGISTERED',repr(REGISTERED_COMMIT)).replace('DATASETS',repr(DATASETS))
 nb=nbformat.v4.new_notebook(cells=[nbformat.v4.new_markdown_cell('<img src="https://raw.githubusercontent.com/senanurcetin/enveda-casmi26/main/site/assets/nexus-midnight.svg" width="360" alt="NEXUS">\n\n# CASMI26: Unopened Family Coverage\n\nA separately preregistered internal query-family study. Reconstruct and protect every prior selected population before new queries; compare fixed C3 candidate budgets and C1 nonlearned scorers. No new supervised fitting, external/laboratory replication, paid services or competition submission.'),nbformat.v4.new_code_cell(setup),nbformat.v4.new_code_cell(execute)],metadata={'kernelspec':{'name':'python3','display_name':'Python 3','language':'python'},'language_info':{'name':'python'}})
 for index,cell in enumerate(nb.cells):cell.id=hashlib.sha256(f'{index}:{cell.source}'.encode()).hexdigest()[:12]
 OUTPUT.mkdir(parents=True,exist_ok=True);nbformat.write(nb,OUTPUT/'notebook.ipynb')
 metadata={'id':KERNEL,'title':'CASMI26: Unopened Family Coverage','code_file':'notebook.ipynb','language':'python','kernel_type':'notebook','is_private':False,'enable_gpu':False,'enable_tpu':False,'enable_internet':False,'competition_sources':['enveda-CASMI26-molecule-id-mass-spectra'],'dataset_sources':DATASETS,'kernel_sources':[]}
 (OUTPUT/'kernel-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
 print('Built new unopened-family CPU study')

if __name__=='__main__':build()
