"""Execute preregistered, unopened internal families; aggregate-only output."""
from __future__ import annotations
import hashlib,json,os,sys,time
from pathlib import Path
from dataclasses import replace,asdict
from datetime import datetime,UTC
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]

def main():
 import numpy as np
 import pandas as pd
 import rdkit
 from casmi.config import ARTIFACT_DIR,POOL_FILE,POOL_FP_FILE,POOL_FP_KEYS_FILE
 from casmi.spectra.store import SpectrumStore
 from casmi.pool.pool import CandidatePool
 from casmi.validate.external_population import prepare_internal
 from casmi.validate.external_purge import select_external_families,blocked_training_identities
 from casmi.validate.unopened_population import select_unopened
 from casmi.validate.scaffold_study import digest,joint_intervals
 from casmi.validate.source_manifest import create_manifest,sha256_file,required_validation_inputs
 from casmi.pool.fingerprint_provenance import verify_fingerprint_provenance
 from casmi.channels.library import LibraryIndex,LibraryConfig
 from casmi.channels.analog import AnalogIndex
 from casmi.chem.keys import MetricKeys
 from casmi.chem.canon import skeleton_key
 from casmi.engine import BaselineEngine,EngineConfig,molecule_mass
 from casmi.generate.expanded import preserve_budget
 from casmi.validate.knockout import Query,include_equivalent_source_rows,knockout_for
 from casmi.validate.metric import summarize
 from casmi.spectra.clean import clean_peaks_fast
 from scripts.run_advanced_study import SearchMemo
 from casmi.validate.knockout import equivalent_keys
 started=time.time()
 protocol_path=ROOT/'protocols/unopened-coverage-v1.json'
 protocol=json.loads(protocol_path.read_text());seed=protocol['seed']
 if rdkit.__version__ not in ('2026.03.3','2026.3.3'):raise ValueError('Pinned metric runtime required')
 old=json.loads((ROOT/'reports/advanced_study.json').read_text())
 purged=json.loads((ROOT/'reports/purged_external_study.json').read_text())
 old_proto=json.loads((ROOT/'protocols/advanced-internal-v1.json').read_text())
 purge_proto=json.loads((ROOT/'protocols/purged-external-v1.json').read_text())
 train=Path(os.environ['CASMI_DATA'])/'train.parquet'
 if sha256_file(train)!=old['raw_inputs']['train.parquet']['sha256']:raise ValueError('Original snapshot mismatch')
 store=SpectrumStore(ARTIFACT_DIR/'store');truth_pool=CandidatePool(pd.read_parquet(ARTIFACT_DIR/'structures.parquet'))
 records,groups,labels,protected=prepare_internal(store.meta,truth_pool,old_proto,old['partitions'],purge_proto,return_candidates=True)
 source=Path(os.environ['CASMI_HISTORY_EXTERNAL'])
 expected=json.loads((ROOT/'reports/gnps_expanded_source_preparation.json').read_text())
 if sha256_file(source)!=expected['parquet_sha256']:raise ValueError('Historical source mismatch')
 selected,external_families=select_external_families(pd.read_parquet(source),purge_proto['split_salt'])
 if digest(sorted(external_families))!=purged['purged_training']['selected_family_sha256']:raise ValueError('External family reconstruction differs')
 full_labels=pd.read_parquet(train,columns=['inchikey14','normalized_smiles']).drop_duplicates()
 blocked=blocked_training_identities(full_labels.itertuples(index=False,name=None),selected,external_families)
 filtered_meta=store.meta[~store.meta.train_ik14.isin(blocked)].reset_index(drop=True)
 filtered_truth=CandidatePool(truth_pool.table[~truth_pool.table.pool_key.isin(blocked)].reset_index(drop=True))
 old_splits,_,_,old_audit=prepare_internal(filtered_meta,filtered_truth,old_proto,old['partitions'],purge_proto,protected=protected)
 for name in ('train','dev','test'):
  if old_audit[name]!=purged['partitions'][name]:raise ValueError('Opened purged-study partitions could not be reconstructed')
 families=set(protected['old_families'])|set(protected['historical_vin'])|set(external_families)
 metrics=set(protected['old_metrics'])|set(selected.metric_key)
 for values in old_splits.values():
  for item in values:families.update((item['family'],item['group']));metrics.add(item['metric'])
 splits=select_unopened(records,families,metrics,protocol['counts'],protocol['split_salt'])
 population={name:{'n_molecules':len(items),'identity_sha256':digest([r['key'] for r in items]),'family_sha256':digest([r['group'] for r in items]),'metric_sha256':digest([r['metric'] for r in items])} for name,items in splits.items()}
 guard={'advanced_partitions_reconstructed':True,'purged_partitions_reconstructed':True,'external_families_reconstructed':True,'protected_family_count':len(families),'protected_metric_count':len(metrics),'selected_protected_family_overlap':0,'selected_protected_metric_overlap':0}
 if any(r['family'] in families or r['group'] in families or r['metric'] in metrics for items in splits.values() for r in items):raise ValueError('Opened population leakage')
 del full_labels,filtered_meta,filtered_truth,selected
 print('Prepared authenticated unopened populations before candidate inspection',flush=True)
 pool=CandidatePool(pd.read_parquet(ARTIFACT_DIR/POOL_FILE));fps=np.load(ARTIFACT_DIR/POOL_FP_FILE,mmap_mode='r',allow_pickle=False)
 fp_keys=pd.read_parquet(ARTIFACT_DIR/POOL_FP_KEYS_FILE).pool_key.to_numpy()
 fp_status=verify_fingerprint_provenance(pool.table,fp_keys,fps,json.loads((ARTIFACT_DIR/'pool_fp_provenance.json').read_text()))
 names=list(required_validation_inputs(pool_tag='',analog=True))+['structures.parquet','pool_source_quality.json']
 manifest=create_manifest(ARTIFACT_DIR,names)
 row_keys=store.meta.train_ik14.astype(str).to_numpy();cfg=EngineConfig(generate=True,gen_consensus=True)
 entropy=LibraryIndex(store,row_keys,cfg.library);cosine=LibraryIndex(store,row_keys,replace(cfg.library,similarity='optimal_cosine'))
 analog=AnalogIndex(store,row_keys,cfg.analog,cache=None);memo=SearchMemo(analog)
 engines={name:BaselineEngine(pool,entropy,replace(cfg,gen_max_analogs=b['seeds'],gen_cap=b['cap'],gen_expand=b['ordered_pairs']),memo,fps) for name,b in protocol['policies'].items()}
 keys=MetricKeys(ARTIFACT_DIR/'unopened_metric_keys.parquet')
 libraries={name:np.flatnonzero((store.meta.ingest_lib==name).to_numpy()) for name in old_proto['query_libraries']}
 outcomes={};details={}
 def covered(candidate_keys,truth,skel):
  return int(any(skeleton_key(pool.smiles_of(k))==skel and keys.get(k,pool.smiles_of(k))==truth for k in candidate_keys))
 def evaluate(items):
  result={name:[] for name in engines};sizes={name:[] for name in engines};ranks={'entropy':[],'optimal_cosine':[]}
  for position,item in enumerate(items):
   rows=groups[item['key']];rng=np.random.default_rng(int(hashlib.sha256(item['key'].encode()).hexdigest()[:16],16))
   chosen=np.sort(rng.choice(rows,5,replace=False)) if len(rows)>5 else rows
   spectra=[]
   for index in chosen:
    mz,intensity=store.peaks(index);row=store.meta.iloc[index]
    spectra.append({'mz':mz,'it':intensity,'precursor_mz':float(row.precursor_mz),'adduct':str(row.adduct),'mode':str(row.ionization_mode),'ce_ev':float(row.ce_ev)})
   query=Query(item['library']+':'+item['key'],item['key'],spectra,chosen)
   mass=molecule_mass(spectra)
   truth,same=equivalent_keys(item['key'],mass,truth_pool,keys,truth_smiles=labels[item['key']])
   candidate_truth,candidates=equivalent_keys(item['key'],mass,pool,keys,truth_smiles=labels[item['key']])
   if truth!=item['metric'] or candidate_truth!=truth:raise ValueError('Metric population mismatch')
   same.update(candidates);same.add(item['key']);include_equivalent_source_rows(query,same,libraries[item['library']],row_keys)
   for s in spectra:s['mz_c'],s['it_c']=clean_peaks_fast(s['mz'],s['it'],s['precursor_mz'],cfg.library.clean)
   c3=knockout_for(query,'C3',entropy.rows_of,same);memo.reset()
   search=memo.search(spectra,mass,c3.excluded_rows)
   ordinary=pool.key[pool.window(mass,cfg.tol_ppm,c3.excluded_keys).idx]
   candidate_sets={};generated_sets={}
   for name,engine in engines.items():
    frame=engine._generated(search,mass,c3.excluded_keys)
    previous=generated_sets['legacy'] if name=='expanded' else generated_sets['expanded'] if name=='extended' else []
    generated_sets[name]=preserve_budget(previous,list(frame.pool_key) if len(frame) else [],protocol['policies'][name]['cap'])
    current=set(ordinary)|set(generated_sets[name])
    candidate_sets[name]=current;result[name].append(covered(current,truth,skeleton_key(labels[item['key']])));sizes[name].append(len(current))
   c1=knockout_for(query,'C1',entropy.rows_of,same);window=pool.window(mass,cfg.tol_ppm)
   for name,index in [('entropy',entropy),('optimal_cosine',cosine)]:
    scored=[(index.score(spectra,k,c1.excluded_rows)[0],abs(float(ppm)),str(k)) for k,ppm in zip(pool.key[window.idx],window.ppm)]
    ordered=[k for _,_,k in sorted(scored,key=lambda v:(-v[0],v[1],v[2]))]
    ranked=keys.dedupe_ranked(ordered,pool.smiles_of,25)
    ranks[name].append(next((i for i,(_,metric) in enumerate(ranked,1) if metric==truth),0))
   if (position+1)%10==0:print(f'Prepared unopened query {position+1}/{len(items)}',flush=True)
  return result,sizes,ranks
 dev=evaluate(splits['dev']);gains=sum(a>b for a,b in zip(dev[0]['extended'],dev[0]['expanded']));losses=sum(a<b for a,b in zip(dev[0]['extended'],dev[0]['expanded']))
 seal={'sealed_at':datetime.now(UTC).isoformat(),'protocol_sha256':sha256_file(protocol_path),'population':population,'protection':guard,'source_manifest_sha256':digest(manifest),'development':{'gains':gains,'losses':losses,'recommendation':'extended' if gains>0 and losses==0 else 'expanded'},'supervised_models_fitted':0}
 print('Sealed unopened populations and development recommendation before test',flush=True)
 test=evaluate(splits['test']);n=len(splits['test']);effect=np.asarray(test[0]['extended'])-np.asarray(test[0]['expanded']);rng=np.random.default_rng(seed)
 boots=effect[rng.integers(n,size=(4000,n))].mean(axis=1)
 rr=lambda r:np.where(np.asarray(r)>0,1/np.maximum(np.asarray(r),1),0)
 secondary=joint_intervals({'C3_coverage_expanded_minus_legacy':np.asarray(test[0]['expanded'])-np.asarray(test[0]['legacy']),'C1_MRR_optimal_cosine_minus_entropy':rr(test[2]['optimal_cosine'])-rr(test[2]['entropy'])},n_boot=4000,seed=seed)
 if digest(create_manifest(ARTIFACT_DIR,names))!=seal['source_manifest_sha256']:raise ValueError('Static source changed during test')
 result={'schema_version':1,'status':'COMPLETE_UNOPENED_INTERNAL_COVERAGE','checked_at':datetime.now(UTC).isoformat(),'protocol':protocol,'seal':seal,'seal_sha256':digest(seal),'protection':guard,'population':population,'fingerprint_provenance':fp_status,'analog_provenance':analog.fresh_index_provenance(),'source_manifest':manifest,'primary':{'n_molecules':n,'delta_coverage':float(effect.mean()),'ci_low':float(np.quantile(boots,.025)),'ci_high':float(np.quantile(boots,.975)),'gains':int((effect>0).sum()),'losses':int((effect<0).sum()),'bootstrap_draws':4000},'secondary':secondary,'coverage':{name:{'covered':sum(values),'n_molecules':n,'mean_candidates':float(np.mean(test[1][name])),'max_candidates':max(test[1][name])} for name,values in test[0].items()},'C1':{name:asdict(summarize(np.asarray(values),n_boot=4000,seed=seed)) for name,values in test[2].items()},'private_details_sha256':digest({'development':dev,'test':test}),'elapsed_seconds':time.time()-started,'environment':{'python':sys.version.split()[0],'rdkit':rdkit.__version__},'limitations':['New internal query-family cohort; no external or laboratory replication.','No learned D3/NP-score ablation or state-of-the-art superiority estimate.','Candidate-coverage improvement does not establish achieved C3 ranking accuracy.','Historical sources and shared task-defined candidate knowledge remain.','Both positive and negative outcomes are retained.']}
 Path(os.environ['CASMI_UNOPENED_OUTPUT']).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print('Completed preregistered unopened-family coverage comparison',flush=True)

if __name__=='__main__':main()
