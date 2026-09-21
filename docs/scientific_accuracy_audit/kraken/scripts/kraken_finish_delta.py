#!/usr/bin/env python3
"""Join both delta extrema and complete independent ensembles into outlier evidence."""
import datetime,json,math
from pathlib import Path
import pandas as pd
from kraken_analyze import ROOT,FIELDS,metrics,save,sha
OUT=ROOT/'delta_interpretation'
def value(r,variant,field):
 v=r.get(variant,{})
 if not v or 'error' in v:return None
 f='buried_volume' if field=='percent_buried_volume' else field
 x=v['orientations'][0][f] if variant.startswith('matched_sut') and f in ['buried_volume','near_vbur','far_vbur'] else v[f]
 return x*100/(4*math.pi*3.5**3/3) if field=='percent_buried_volume' else x

def main():
 OUT.mkdir(exist_ok=False)
 mf=json.loads((ROOT/'delta_outliers/manifest_frozen.json').read_text());assert sha(ROOT/'delta_outliers/reference_raw.jsonl')==mf['raw_sha256']
 frame=pd.read_csv(ROOT/'analysis/conformers.csv');all_top=json.loads((ROOT/'interpretation/all_top20_outlier_dossiers.json').read_text());full={v['id']:v for v in map(json.loads,(ROOT/'full_analytic_reference/reference_raw.jsonl').open())};refs={v['id']:v for v in map(json.loads,(ROOT/'morfeus_outliers/reference_raw.jsonl').open())};refs.update({v['id']:v for v in map(json.loads,(ROOT/'delta_outliers/reference_raw.jsonl').open())});mids=set(json.loads((ROOT/'delta_outliers/selection_frozen.json').read_text())['ligands'])
 all_aggregate=[];dossiers={};variants=['matched_sut_001','matched_sut_0001','primary_001','primary_0001']
 for name,rows in all_top.items():
  if not name.endswith('_delta'):continue
  records=[]
  for row in rows:
   mid=row['molecule_id'];f=row['descriptor'];g=frame[frame.molecule_id==mid];lo=g.loc[g[f].idxmin()];hi=g.loc[g[f].idxmax()];ext=[]
   for which,p in [('min',lo),('max',hi)]:
    key=p['id'];cid=int(p['conformer_id']);ext.append({'extremum':which,'id':key,'conformer_id':cid,'sut_value':float(p[f]),'input':f'primary/conformers/{cid}.body','input_sha256':sha(ROOT/'primary/conformers'/f'{cid}.body'),'independent_primary_analytic':full[key],'independent_convention_variants':refs.get(key)})
   comparisons=[]
   if not f.startswith(('sterimol','pyr')):
    for variant in variants:
     vals=[(k,value(refs[k],variant,f)) for k in g.id];good=[(k,v) for k,v in vals if v is not None]
     if len(good)!=len(g):comparisons.append({'variant':variant,'error':'Incomplete reference ensemble','count':len(good),'expected':len(g)});continue
     mn=min(good,key=lambda p:p[1]);mx=max(good,key=lambda p:p[1]);delta=mx[1]-mn[1];comparisons.append({'variant':variant,'N_conformers':len(good),'independent_minimum':mn[1],'independent_minimum_id':mn[0],'independent_maximum':mx[1],'independent_maximum_id':mx[0],'delta':delta,'published_delta':row['reference'],'native_delta':row['sut'],'residual_vs_published':delta-row['reference'],'native_minus_reference_delta':row['sut']-delta})
    row['independent_complete_ensemble']=comparisons
    primary=next((v for v in comparisons if v['variant']=='primary_0001' and 'delta' in v),None);matched=next((v for v in comparisons if v['variant']=='matched_sut_001' and 'delta' in v),None)
    row['likely_cause']=f"Native delta {row['sut']:.17g}; published {row['reference']:.17g}. Matched-default independent delta {matched['delta'] if matched else 'unavailable'}; primary-center/density independent delta {primary['delta'] if primary else 'unavailable'}. Four complete-ensemble convention variants quantify center/grid effects; remaining primary-minus-published residual is not assigned automatically to either implementation."
    row['result_classification']='UNCERTAIN for any remaining historical-reference mismatch; independently measured convention differences are quantified separately.'
    all_aggregate.extend({'metric':name,'molecule_id':mid,**v} for v in comparisons)
   else:
    row['result_classification']='Independent complete-ensemble delta available; residual historical-source cause UNCERTAIN where disagreement persists.'
   row.pop('input_id',None);row.pop('input_path',None);row['input_ids']=[v['id'] for v in ext];row['input_paths']=[v['input'] for v in ext];row['native_extremizers']=ext;row['selected_conformer_diagnostics']={v['extremum']:{'id':v['id'],'primary_analytic':v['independent_primary_analytic'],'convention_variants':v['independent_convention_variants']} for v in ext};row['limitation']='Both native extremizers are mapped explicitly; independent full ensembles select their own extrema. Original full-precision historical geometries/energies/dependency versions are unavailable; remaining published residuals do not establish which source is incorrect.';records.append(row)
  dossiers[name]=records
 save(OUT/'delta_outlier_dossiers.json',dossiers);save(OUT/'delta_variant_aggregates.json',all_aggregate)
 # Authoritative replacement for the derivative dossier; preserve prior analysis rather than erase it.
 previous=ROOT/'interpretation/all_top20_outlier_dossiers.json';target=OUT/'all_top20_outlier_dossiers.json';save(target,all_top)
 # Combined reference kernel/percent evidence across all distinct diagnostic cases.
 sut={v['id']:v for v in map(json.loads,(ROOT/'sut/stdout.jsonl').open()) if v['id'] in refs};kernel=[]
 for key,r in refs.items():
  s=sut[key]['buried_volume']
  if 'error' in s:continue
  for f in list(FIELDS)[:8]+['percent_buried_volume']:
   ref=value(r,'matched_sut_001',f)
   if ref is None:continue
   kernel.append({'id':key,'descriptor':f,'reference':ref,'sut':s[f],'residual':s[f]-ref})
 k=pd.DataFrame(kernel);k.to_csv(OUT/'expanded_matched_default.csv',index=False);summary={}
 for f,g in k.groupby('descriptor'):
  x=metrics(g.reference,g.sut);x['maximum_case']=g.loc[g.residual.abs().idxmax()].to_dict();summary[f]=x
 save(OUT/'expanded_matched_default_metrics.json',summary)
 save(OUT/'manifest.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'script_sha256':sha(Path(__file__)),'frozen_delta_manifest_sha256':sha(ROOT/'delta_outliers/manifest_frozen.json'),'original_dossier_sha256':sha(previous),'updated_dossier_sha256':sha(target),'delta_metric_count':len(dossiers),'delta_entries':sum(map(len,dossiers.values())),'all_delta_entries_have_two_inputs':all(len(v['input_ids'])==2 for rows in dossiers.values() for v in rows),'outputs':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()}})
 print('Completed',len(dossiers),'delta metrics',sum(map(len,dossiers.values())),'entries')
 for v in all_aggregate:
  if v['molecule_id']==369 and v['metric'] in ['buried_volume_delta','far_vbur_delta']:print(json.dumps(v))
 print('MATCHEDMAX',json.dumps({k:(v['N'],v['maximum_absolute_error'],v['maximum_case']['id']) for k,v in summary.items()}))
if __name__=='__main__':main()
