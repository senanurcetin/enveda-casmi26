"""Authenticate the unopened-family comparison; original studies remain frozen."""
from pathlib import Path
from datetime import datetime,UTC
import hashlib,json,math,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from scripts.build_unopened_coverage_notebook import KERNEL,DATASETS,REGISTERED_COMMIT
from scripts.check_kaggle_publications import source_hash
from scripts.research_snapshot import snapshot_hashes
from casmi.validate.scaffold_study import digest

def require(value,message):
 if not value:raise ValueError(message)

def validate(result):
 protocol=json.loads((ROOT/'protocols/unopened-coverage-v1.json').read_text())
 require(result['status']=='COMPLETE_UNOPENED_INTERNAL_COVERAGE' and result['protocol']==protocol,'Registered result design mismatch')
 require(result['source_code_sha256']==snapshot_hashes(ROOT/'kaggle/audits/unopened-coverage/notebook.ipynb'),'Executed inventory mismatch')
 seal=result['seal'];require(digest(seal)==result['seal_sha256'],'Seal mismatch')
 require(seal['protocol_sha256']==hashlib.sha256((ROOT/'protocols/unopened-coverage-v1.json').read_text().encode()).hexdigest(),'Protocol binding differs')
 require(seal['supervised_models_fitted']==0 and seal['population']==result['population'] and seal['protection']==result['protection'],'Pre-test population binding differs')
 require(datetime.fromisoformat(seal['sealed_at'])<=datetime.fromisoformat(result['checked_at']),'Seal ordering differs')
 guard=result['protection']
 require(all(guard[k] is True for k in ['advanced_partitions_reconstructed','purged_partitions_reconstructed','external_families_reconstructed']),'Prior populations were not reconstructed')
 require(guard['selected_protected_family_overlap']==guard['selected_protected_metric_overlap']==0,'Opened query population overlap')
 require(guard['protected_family_count']>=1000 and guard['protected_metric_count']>=1550,'Incomplete prior population protection')
 for name,count in protocol['counts'].items():
  block=result['population'][name];require(block['n_molecules']==count,'Population denominator differs')
  require(all(isinstance(block[k],str) and len(block[k])==64 for k in ['identity_sha256','family_sha256','metric_sha256']),'Population hashes absent')
 n=protocol['counts']['test'];coverage=result['coverage']
 require(set(coverage)==set(protocol['policies']),'Missing fixed budget outcome')
 for name,cell in coverage.items():
  require(cell['n_molecules']==n and type(cell['covered']) is int and 0<=cell['covered']<=n,'Invalid coverage numerator')
  require(0<=cell['mean_candidates']<=cell['max_candidates'],'Invalid candidate cost')
 require(coverage['legacy']['covered']<=coverage['expanded']['covered']<=coverage['extended']['covered'],'Nested candidate coverage lost')
 primary=result['primary'];delta=(coverage['extended']['covered']-coverage['expanded']['covered'])/n
 require(primary['n_molecules']==n and primary['losses']==0 and primary['gains']==round(delta*n) and math.isclose(primary['delta_coverage'],delta,abs_tol=1e-12) and primary['bootstrap_draws']==4000,'Primary contrast accounting differs')
 require(-1<=primary['ci_low']<=primary['ci_high']<=1,'Invalid primary interval')
 require(set(result['C1'])=={'entropy','optimal_cosine'},'Missing scorer')
 for cell in result['C1'].values():
  require(cell['n']==n and all(math.isfinite(cell[k]) and 0<=cell[k]<=1 for k in ['mrr','ci_low','ci_high','top1','top10','recall_at_25']),'Invalid C1 summary')
  require(cell['top1']<=cell['top10']<=cell['recall_at_25'] and cell['top1']<=cell['mrr']<=cell['recall_at_25'],'Impossible C1 metric')
 joint=result['secondary'];require(set(joint['outcomes'])==set(protocol['secondary']) and joint['family_size']==2 and joint['n_molecules']==n and joint['n_bootstrap']==4000,'Secondary family differs')
 effects={'C3_coverage_expanded_minus_legacy':(coverage['expanded']['covered']-coverage['legacy']['covered'])/n,'C1_MRR_optimal_cosine_minus_entropy':result['C1']['optimal_cosine']['mrr']-result['C1']['entropy']['mrr']}
 for name,value in effects.items():
  cell=joint['outcomes'][name];require(math.isclose(cell['delta_mrr'],value,abs_tol=1e-12) and math.isclose(cell['ci_low'],value-joint['joint_radius'],abs_tol=1e-12) and math.isclose(cell['ci_high'],value+joint['joint_radius'],abs_tol=1e-12),'Secondary absolute/paired inconsistency')
 dev=seal['development'];require(dev['losses']==0 and dev['recommendation']==('extended' if dev['gains']>0 else 'expanded'),'Development choice differs')
 require(result['analog_provenance']['cache_used'] is False and result['fingerprint_provenance']['status']=='POOL_FINGERPRINTS_VERIFIED','Invalid reference provenance')
 require(digest(result['source_manifest'])==seal['source_manifest_sha256'],'Static source manifest seal differs')
 require(result['execution']=={'registered_commit':REGISTERED_COMMIT,'dataset_sources':DATASETS,'competition_source':'enveda-CASMI26-molecule-id-mass-spectra','gpu_enabled':False,'internet_enabled':False,'supervised_models_fitted':0,'submitted_predictions':False},'Execution attachments/policy differ')
 original=json.loads((ROOT/'reports/advanced_study.json').read_text())
 for name in ['train.parquet','coconut.csv','pubchemlite.csv']:require(result['raw_inputs'][name]==original['raw_inputs'][name],'Pinned input mismatch')
 source=json.loads((ROOT/'reports/gnps_expanded_source_preparation.json').read_text());require(result['raw_inputs']['gnps_cc0.parquet']=={'sha256':source['parquet_sha256'],'bytes':source['parquet_bytes']},'Historical external input mismatch')
 def privacy(value):
  if isinstance(value,dict):
   require(not set(value)&{'qid','smiles','query_ids','accessions','selected_rows','ranks','ms2_mzs'},'Private per-query data exported')
   for child in value.values():privacy(child)
  elif isinstance(value,list):
   for child in value:privacy(child)
 privacy(result)

def main():
 from kaggle.api.kaggle_api_extended import KaggleApi
 api=KaggleApi();api.authenticate()
 history=json.loads((ROOT/'reports/unopened_coverage_execution_history.json').read_text())
 with tempfile.TemporaryDirectory(prefix='casmi-unopened-') as temporary:
  folder=Path(temporary);api.kernels_pull(KERNEL,path=temporary,metadata=True,quiet=True)
  metadata=json.loads((folder/'kernel-metadata.json').read_text());remote=json.loads((folder/metadata['code_file']).read_text())
  require(source_hash(remote)==history['attempts'][-1]['source_sha256'],'Latest registered public source mismatch')
  state=api.kernels_status(KERNEL);state=state if isinstance(state,dict) else json.loads(str(state))
  require(state['status']=='COMPLETE' and not metadata['is_private'] and not metadata['enable_gpu'] and not metadata['enable_internet'] and not metadata['kernel_sources'] and sorted(metadata['dataset_sources'])==sorted(DATASETS) and metadata['competition_sources']==['enveda-CASMI26-molecule-id-mass-spectra'],'Completed public CPU attachments not authenticated')
  api.kernels_output(KERNEL,path=temporary,file_pattern=r'^unopened_coverage_study\.json$',quiet=True)
  artifact=folder/'unopened_coverage_study.json';require(artifact.stat().st_size<=2000000,'Unbounded aggregate')
  raw=artifact.read_bytes();validate(json.loads(raw))
 (ROOT/'reports/unopened_coverage_study.json').write_bytes(raw)
 record={'checked_at':datetime.now(UTC).isoformat(),'status':'pass','source_sha256':source_hash(remote),'aggregate_sha256':hashlib.sha256(raw).hexdigest(),'url':f'https://www.kaggle.com/code/{KERNEL}','scope':'Unopened internal query-family coverage and nonlearned C1 comparison; no laboratory replication or learned-model improvement claim'}
 (ROOT/'reports/unopened_coverage_publication.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))

if __name__=='__main__':main()
