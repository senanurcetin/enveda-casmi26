"""Render only authenticated unopened-family aggregates, preserving null findings."""
from pathlib import Path
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.check_unopened_coverage_publication import validate

def build(tables_only=False):
 raw=(ROOT/'reports/unopened_coverage_study.json').read_bytes();result=json.loads(raw);validate(result)
 record=json.loads((ROOT/'reports/unopened_coverage_publication.json').read_text())
 if record['status']!='pass' or record['aggregate_sha256']!=hashlib.sha256(raw).hexdigest():raise ValueError('Authenticated aggregate required')
 primary=result['primary'];dev=result['seal']['development']
 lines=['![NEXUS](../site/assets/nexus-midnight.svg)','','# Unopened internal-family coverage and reference comparison','',
 'This separately registered study protects all previously selected query families and evaluates new internal families. It fits no supervised model and does not establish external or laboratory replication. Candidate truth coverage is not achieved C3 ranking accuracy.','',
 '[Protocol](../protocols/unopened-coverage-v1.json) · [Aggregate](unopened_coverage_study.json) · [Authenticated publication](unopened_coverage_publication.json)','',
 '![NEXUS unopened-family comparison](../docs/assets/unopened-coverage.png)','',
 '[Vector figure](../docs/assets/unopened-coverage.svg) · [PDF](../docs/assets/unopened-coverage.pdf)','',
 f"Development molecules: 100; test molecules: 300. Development recommendation: **{dev['recommendation']}**, with {dev['gains']} gains and {dev['losses']} losses versus expanded. Both larger budgets remain reported regardless of this choice.",'',
 '## C3 candidate coverage','','| Fixed policy | Covered / all test molecules | Mean candidates | Maximum candidates |','|---|---:|---:|---:|']
 for name,v in result['coverage'].items():lines.append(f"| {name} | {v['covered']}/{v['n_molecules']} | {v['mean_candidates']:.1f} | {v['max_candidates']} |")
 lines+=['',f"Primary extended-minus-expanded coverage change: **{primary['delta_coverage']:+.4f} [{primary['ci_low']:+.4f}, {primary['ci_high']:+.4f}]**, paired 95% percentile interval from 4,000 molecule resamples. Gains: {primary['gains']}; losses: {primary['losses']}. All 300 molecules remain in every denominator.",'',
 'The proposal caps apply to generated candidates; ordinary mass-compatible pool candidates are additional. The extended policy preserves earlier proposals before filling its cap. A primary interval touching zero does not establish a strictly positive population benefit.','',
 '## C1 reference scoring','','| Scorer | MRR@25 [95% CI] | Top 1 | Top 10 | Hit 25 |','|---|---:|---:|---:|---:|']
 for name,v in result['C1'].items():lines.append(f"| {name} | {v['mrr']:.4f} [{v['ci_low']:.4f}, {v['ci_high']:.4f}] | {v['top1']:.4f} | {v['top10']:.4f} | {v['recall_at_25']:.4f} |")
 lines+=['','C1 removes all same-library metric-equivalent truth acquisitions. Direct scorers need an eligible cross-library reference to supply truth spectral evidence; queries without one remain in the denominator. Raw-intensity optimal cosine and intensity-weighted entropy have different objectives. An absolute interval does not establish a between-scorer gain.','',
 '## Joint secondary comparisons','','| Comparison | Difference [joint 95% CI] |','|---|---:|']
 for name,v in result['secondary']['outcomes'].items():lines.append(f"| {name} | {v['delta_mrr']:+.4f} [{v['ci_low']:+.4f}, {v['ci_high']:+.4f}] |")
 lines+=['','The two registered secondary comparisons share a molecule-paired max-deviation radius. Coverage differences retain binary-coverage units even though the generic interval record uses the field name delta_mrr. No adverse result is suppressed.','',
 '## Limits','']+['- '+v for v in result['limitations']]+['','No learned D3/NP-score ablation was performed. These unresolved effects cannot be declared resolved by this study. Historical studies use different populations; direct score comparisons do not isolate causal improvement. Source annotations remain computationally checked, not independently authenticated laboratory measurements.','']
 (ROOT/'reports/unopened_coverage_review.md').write_text('\n'.join(lines),encoding='utf-8')
 if not tables_only:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  fig,axes=plt.subplots(1,2,figsize=(12,6),facecolor='#0B0F19');fig.suptitle('NEXUS / CASMI26 · Unopened internal families',color='#00E0FF',fontsize=17)
  for ax in axes:
   ax.set_facecolor('#0B0F19');ax.tick_params(colors='#B9C5CE')
   for spine in ax.spines.values():spine.set_color('#344456')
   ax.set_ylim(0,1);ax.grid(axis='y',alpha=.2)
  names=list(result['coverage']);values=[result['coverage'][n]['covered']/300 for n in names]
  axes[0].bar(names,values,color=['#B5A4FF','#00E0FF','#4EA8FF']);axes[0].set_title('C3 truth coverage / all 300 molecules',color='white')
  for i,(name,value) in enumerate(zip(names,values)):
   axes[0].text(i,value+.02,f"{result['coverage'][name]['covered']}/300",ha='center',color='#B9C5CE')
  names=list(result['C1']);values=[result['C1'][n]['mrr'] for n in names]
  lo=[v-result['C1'][n]['ci_low'] for n,v in zip(names,values)];hi=[result['C1'][n]['ci_high']-v for n,v in zip(names,values)]
  axes[1].bar(names,values,yerr=[lo,hi],color=['#00E0FF','#B5A4FF'],error_kw={'ecolor':'#B9C5CE','capsize':4});axes[1].set_title('C1 MRR@25 · bootstrap 95% intervals',color='white')
  fig.text(.03,.065,f"Primary coverage change: {primary['delta_coverage']:+.4f} [{primary['ci_low']:+.4f}, {primary['ci_high']:+.4f}]; no positive gain established when interval touches zero.",color='#B9C5CE',fontsize=9)
  fig.text(.03,.025,'Internal unopened-family cohort; no supervised fitting, external validation or laboratory replication.',color='#B9C5CE',fontsize=9)
  fig.tight_layout(rect=(0,.10,1,.9))
  for extension in ['png','svg','pdf']:fig.savefig(ROOT/f'docs/assets/unopened-coverage.{extension}',dpi=180,facecolor=fig.get_facecolor())
  plt.close(fig)
 print('Rendered all authenticated unopened-family outcomes')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--tables-only',action='store_true');build(p.parse_args().tables_only)
