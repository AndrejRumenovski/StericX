#!/usr/bin/env python3
"""Compare frozen full-precision SUT measurements against fresh primary records."""
from __future__ import annotations
import collections,datetime,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
FIELDS={
 'buried_volume':('buried_volume','buried_volume','vbur_vbur','Å³'),
 'qvbur_min':('buried_volume','qvbur_min','vbur_qvbur_min','Å³'),
 'qvbur_max':('buried_volume','qvbur_max','vbur_qvbur_max','Å³'),
 'max_delta_qvbur':('buried_volume','max_delta_qvbur','vbur_max_delta_qvbur','Å³'),
 'ovbur_min':('buried_volume','ovbur_min','vbur_ovbur_min','Å³'),
 'ovbur_max':('buried_volume','ovbur_max','vbur_ovbur_max','Å³'),
 'near_vbur':('buried_volume','near_vbur','vbur_near_vbur','Å³'),
 'far_vbur':('buried_volume','far_vbur','vbur_far_vbur','Å³'),
 'sterimol_l':('sterimol_coordination','l','sterimol_L','Å'),
 'sterimol_b1':('sterimol_coordination','b1','sterimol_B1','Å'),
 'sterimol_b5':('sterimol_coordination','b5','sterimol_B5','Å'),
 'pyr_p':('pyramidalization','pyr_p','pyr_P','dimensionless'),
 'pyr_alpha':('pyramidalization','pyr_alpha','pyr_alpha','degrees'),
 'percent_buried_volume':('buried_volume','percent_buried_volume','vbur_vbur','percentage points; reference derived from absolute Vbur'),
}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def metrics(x,y):
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);d=y-x
 n=len(x);sst=float(np.sum((x-x.mean())**2));r2=1-float(d@d)/sst if sst>0 else None
 slope,intercept=np.polyfit(x,y,1) if n>1 and np.ptp(x)>0 else (None,None)
 pearson=float(np.corrcoef(x,y)[0,1]) if n>1 and np.ptp(x)>0 and np.ptp(y)>0 else None
 return {'N':n,'MAE':float(np.abs(d).mean()),'RMSE':float(np.sqrt(np.mean(d*d))),'maximum_absolute_error':float(np.abs(d).max()),'median_absolute_error':float(np.median(np.abs(d))),'mean_signed_error':float(d.mean()),'R2_one_to_one':r2,'slope':float(slope) if slope is not None else None,'intercept':float(intercept) if intercept is not None else None,'pearson_r':pearson}
def main():
 out=ROOT/'analysis';out.mkdir(exist_ok=False)
 frozen=json.loads((ROOT/'sut/manifest_frozen.json').read_text())
 assert frozen['returncode']==0 and frozen['frozen_before_reference_comparison']
 for name,m in frozen['raw_output_files'].items():assert sha(ROOT/'sut'/name)==m['sha256']
 inventory={r['id']:r for r in map(json.loads,(ROOT/'prepared_all/inventory.jsonl').open()) if 'id' in r}
 ligands=json.loads((ROOT/'primary/ligands_manifest.json').read_text())['ligands']
 errors=[];records=[];bond_diffs=[];seen=set()
 with (ROOT/'sut/stdout.jsonl').open() as stream:
  for line in stream:
   r=json.loads(line);key=r['id'];assert key not in seen;seen.add(key)
   inp=inventory[key];row={k:inp[k] for k in ['id','molecule_id','conformer_id','topology_category','P_H_count','donor','reference','atom_count','sanitization_flag']}
   detected=[v['index'] for v in r['bonded_neighbors']]
   if set(detected)!=set(inp['neighbors']):bond_diffs.append({'id':key,'topological':inp['neighbors'],'sut_geometric':detected})
   for operation in ['buried_volume','sterimol_coordination','pyramidalization']:
    if 'error' in r[operation]:errors.append({'id':key,'molecule_id':inp['molecule_id'],'conformer_id':inp['conformer_id'],'operation':operation,'error':r[operation]['error'],'independent_neighbors':inp['neighbors'],'sut_neighbors':detected})
   for column,(operation,field,_,_) in FIELDS.items():row[column]=r[operation].get(field)
   if isinstance(r['center'],list):row.update({f'center_{axis}':r['center'][i] for i,axis in enumerate('xyz')})
   records.append(row)
 assert seen==set(inventory),(len(seen),len(inventory))
 frame=pd.DataFrame(records);frame.to_csv(out/'conformers.csv',index=False)
 save(out/'operation_errors.json',errors);save(out/'connectivity_disagreements.json',bond_diffs)
 rows=[];excluded=[];reference_disagreements=[]
 historic=pd.read_csv(ROOT/'raw_sources/ni_hda_raw.body').set_index('Unnamed: 0')
 for ligand in ligands:
  mid=ligand['molecule_id'];g=frame[frame.molecule_id==mid];expected=len(ligand.get('conformer_ids',[]))
  if expected==0:
   excluded.append({'molecule_id':mid,'reason':'no_DFT_geometry_available','all_properties':True});continue
  assert len(g)==expected
  published_path=ROOT/'primary/published_dft'/f'{mid}.body'
  raw=json.loads(published_path.read_text())
  published={r['property']:r for r in raw} if isinstance(raw,list) else {}
  for column,(_,_,prop,units) in FIELDS.items():
   values=g[column].to_numpy(dtype=float)
   if not np.isfinite(values).all():
    excluded.append({'molecule_id':mid,'property':column,'reason':'incomplete_SUT_property','available':int(np.isfinite(values).sum()),'expected':expected});continue
   for reduction in ['min','max','delta','vburminconf']:
    ref=published.get(prop,{}).get(reduction)
    if ref is None:
     excluded.append({'molecule_id':mid,'property':column,'reduction':reduction,'reason':'published_reference_missing'});continue
    if reduction=='min':idx=int(np.argmin(values));value=float(values[idx])
    elif reduction=='max':idx=int(np.argmax(values));value=float(values[idx])
    elif reduction=='delta':idx=None;value=float(values.max()-values.min())
    else:
     bv=g.buried_volume.to_numpy(dtype=float)
     if not np.isfinite(bv).all():
      excluded.append({'molecule_id':mid,'property':column,'reduction':reduction,'reason':'vburminconf_unavailable'});continue
     idx=int(np.argmin(bv));value=float(values[idx])
    source_ref=float(ref)
    if column=='percent_buried_volume':ref=source_ref*100/(4*math.pi*3.5**3/3)
    old_key=f'{prop}_{reduction}'
    if old_key in historic and pd.notna(historic.loc[mid,old_key]):
     previous=float(historic.loc[mid,old_key])
     if previous!=source_ref:reference_disagreements.append({'molecule_id':mid,'property':prop,'reduction':reduction,'api':source_ref,'primary_csv':previous,'difference':source_ref-previous})
    rows.append({'molecule_id':mid,'descriptor':column,'reduction':reduction,'metric':f'{column}_{reduction}','units':units,'reference':float(ref),'sut':value,'residual':value-float(ref),'absolute_error':abs(value-float(ref)),'conformers':expected,'sut_extremizer_conformer_id':int(g.iloc[idx].conformer_id) if idx is not None else None,'P_H_count':int(g.P_H_count.max()),'topology_category':g.iloc[0].topology_category})
 comparisons=pd.DataFrame(rows);comparisons.to_csv(out/'comparisons.csv',index=False)
 save(out/'excluded_comparisons.json',excluded);save(out/'reference_snapshot_differences.json',reference_disagreements)
 all_metrics={name:{**metrics(g.reference,g.sut),'units':g.iloc[0].units} for name,g in comparisons.groupby('metric',sort=True)}
 save(out/'metrics.json',all_metrics)
 worst={name:g.nlargest(20,'absolute_error').astype(object).where(pd.notna(g.nlargest(20,'absolute_error')),None).to_dict('records') for name,g in comparisons.groupby('metric',sort=True)}
 save(out/'top20_by_descriptor.json',worst)
 headline=worst['max_delta_qvbur_min'];save(out/'headline_top20.json',headline)
 coverage={'universe':len(ligands),'available_DFT_ligands':sum(bool(r.get('dft_available')) for r in ligands),'available_DFT_conformers':len(frame),'fresh_SUT_records':len(seen),'complete_ligands_by_property':{c:int(comparisons[(comparisons.descriptor==c)&(comparisons.reduction=='min')].molecule_id.nunique()) for c in FIELDS},'operation_error_count':len(errors),'operation_errors_by_message':dict(collections.Counter((r['operation']+': '+r['error']) for r in errors)),'independent_topology_SUT_connectivity_disagreements':len(bond_diffs),'metrics_count':len(all_metrics),'no_outliers_removed':True,'Boltzmann_reductions':'Not reproduced: primary geometry exports and conformer-data endpoint do not provide per-conformer energies/weights. Published ligand-level Boltzmann values cannot serve as their own inputs.'}
 save(out/'coverage.json',coverage)
 for group_name,names in [('buried_volume_min',[f'{c}_min' for c in list(FIELDS)[:8]]),('sterimol_extrema',[f'{c}_{red}' for red in ['min','max'] for c in ['sterimol_l','sterimol_b1','sterimol_b5']]),('pyramidalization_extrema',[f'{c}_{red}' for red in ['min','max'] for c in ['pyr_p','pyr_alpha']])]:
  fig,axes=plt.subplots(len(names),2,figsize=(11,2.3*len(names)),squeeze=False)
  for i,name in enumerate(names):
   g=comparisons[comparisons.metric==name];m=all_metrics[name]
   a,b=axes[i];a.scatter(g.reference,g.sut,s=6,alpha=.45)
   lo=min(g.reference.min(),g.sut.min());hi=max(g.reference.max(),g.sut.max());a.plot([lo,hi],[lo,hi],color='black',lw=.7)
   a.set(xlabel='Published reference',ylabel='StericX',title=f"{name}: N={m['N']}, R²={m['R2_one_to_one']:.5f}")
   b.scatter(g.reference,g.residual,s=6,alpha=.45);b.axhline(0,color='black',lw=.7)
   b.set(xlabel='Published reference',ylabel='SUT − reference',title=f"MAE={m['MAE']:.5g}; max |error|={m['maximum_absolute_error']:.5g}")
  fig.tight_layout();fig.savefig(out/f'{group_name}_residuals.png',dpi=150);plt.close(fig)
 save(out/'manifest.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'script_sha256':sha(Path(__file__)),'frozen_SUT_manifest_sha256':sha(ROOT/'sut/manifest_frozen.json'),'frozen_preparation_sha256':sha(ROOT/'prepared_all/manifest.json'),'metrics_definition':'R²=1−Σ(SUT−reference)²/Σ(reference−mean(reference))²; slope/intercept fit SUT against reference solely as diagnostic, never calibrate SUT; every complete ligand retained.','outputs':{p.name:sha(p) for p in out.iterdir() if p.is_file()}})
 print(json.dumps(coverage,indent=2));print('HEADLINE',all_metrics['max_delta_qvbur_min']);print('TOP20',[(r['molecule_id'],r['absolute_error']) for r in headline])
if __name__=='__main__':main()
