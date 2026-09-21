#!/usr/bin/env python3
"""Interpret already-frozen SUT/reference observations; no production kernels."""
import collections, datetime, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from kraken_analyze import ROOT, FIELDS, metrics, sha, save
OUT=ROOT/'interpretation'

def finite_metrics(rows,x='reference',y='sut'):
 g=pd.DataFrame(rows).dropna(subset=[x,y]); return metrics(g[x],g[y]) if len(g) else {'N':0}

def main():
 OUT.mkdir(exist_ok=False)
 reference={r['id']:r for r in map(json.loads,(ROOT/'full_analytic_reference/reference_raw.jsonl').open())}
 diagnostic={r['id']:r for r in map(json.loads,(ROOT/'morfeus_outliers/reference_raw.jsonl').open())}
 selected={r['request']['id']:r for r in map(json.loads,(ROOT/'morfeus_outliers/selected_inputs.jsonl').open())}
 inventory={r['id']:r for r in map(json.loads,(ROOT/'prepared_all/inventory.jsonl').open()) if 'id' in r}
 direct=[];analytic=[];primary=[];error_details=[]
 with (OUT/'same_center_analytic_raw.jsonl').open('w') as dst:
  for s in map(json.loads,(ROOT/'sut/stdout.jsonl').open()):
   key=s['id'];meta=inventory[key];a=reference[key];base={'id':key,'molecule_id':meta['molecule_id'],'conformer_id':meta['conformer_id']}
   if isinstance(s['center'],list):
    center=np.asarray(s['center']);p=np.asarray([a['position'] for a in s['atoms']]);r=np.asarray([a['vdw_radius'] for a in s['atoms']]);u=p[meta['donor']]-center;u=u/np.linalg.norm(u);v=p-center;z=v@u
    # Direct dot products and perpendicular residuals, no quaternion rotation.
    values={'sterimol_l':float(np.max(z+r)+.4),'sterimol_b5':float(np.max(np.linalg.norm(v-z[:,None]*u,axis=1)+r))}
    row={**base,'reference':values,'sut':{f:s['sterimol_coordination'].get(k) for f,k in [('sterimol_l','l'),('sterimol_b5','b5')]},'axis':u.tolist()};dst.write(json.dumps(row,allow_nan=False)+'\n')
    for f,val in values.items():analytic.append({**base,'descriptor':f,'reference':val,'sut':row['sut'][f],'axis_cos_Z':float(u[2])})
   for op,label in [('sterimol_primary','sterimol'),('pyr_primary_nearest','pyr_nearest'),('pyr_topology','pyr_topology')]:
    if 'error' in a[op]:error_details.append({**base,'operation':op,'error':a[op]['error']});continue
    for f,val in a[op].items():
     operation,field,_,_=FIELDS[f];primary.append({**base,'descriptor':f,'variant':label,'reference':val,'sut':s[operation].get(field)})
   if key not in diagnostic:continue
   d=diagnostic[key]
   for variant in ['matched_sut_001','matched_sut_0001','primary_001','primary_0001','ster_primary_center_bondi','ster_primary_center_paton','ster_sut_center_sut_radii','ster_sut_center_paton','pyr_f32_topology','pyr_f64_topology','pyr_f64_primary_neighbors']:
    if variant not in d:continue
    if 'error' in d[variant]:error_details.append({**base,'operation':variant,'error':d[variant]['error']});continue
    for f,val in d[variant].items():
     if f=='orientations':continue
     operation,field,_,_=FIELDS[f]
     # Same-center matched SUT comparisons must use its FIRST orientation for total/near/far.
     chosen=d[variant]['orientations'][0][f] if variant.startswith('matched_sut') and f in ['buried_volume','near_vbur','far_vbur'] else val
     direct.append({**base,'descriptor':f,'variant':variant,'reference':chosen,'sut':s[operation].get(field),'center_shift_A':d.get('center_shift_A'),'center_angle_degrees':d.get('center_angle_degrees')})
 save(OUT/'analytic_raw_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'file_sha256':sha(OUT/'same_center_analytic_raw.jsonl'),'source_sut_sha256':sha(ROOT/'sut/stdout.jsonl'),'script_sha256':sha(Path(__file__)),'method':'Independent closed-form dot product axial and perpendicular projection using exact observed f32 coordinates/center/radii promoted to f64. No SUT kernel; +0.40 L convention. Raw values frozen before error aggregation.'})
 for name,rows in [('same_center_analytic',analytic),('full_primary_conformers',primary),('direct_morfeus',direct)]:
  frame=pd.DataFrame(rows);frame['residual']=frame.sut-frame.reference;frame['absolute_error']=frame.residual.abs();frame.to_csv(OUT/(name+'.csv'),index=False)
  groups=['variant','descriptor'] if 'variant' in frame else ['descriptor'];summary={}
  for keys,g in frame.groupby(groups):
   k='/'.join(keys) if isinstance(keys,tuple) else str(keys);valid=g.dropna(subset=['sut','reference']);m=finite_metrics(valid.to_dict('records'));m['maximum_case']=valid.loc[valid.absolute_error.idxmax()].to_dict() if len(valid) else None;summary[k]=m
  save(OUT/(name+'_metrics.json'),summary)
 save(OUT/'reference_errors.json',error_details)
 # Complete primary reference ensemble reductions vs publication, independently of SUT success.
 prim=pd.DataFrame(primary);comparisons=pd.read_csv(ROOT/'analysis/comparisons.csv');pub={}
 for mid in sorted(prim.molecule_id.unique()):pub[mid]={r['property']:r for r in json.loads((ROOT/'primary/published_dft'/f'{mid}.body').read_text())}
 aggregates=[];counts=collections.Counter(v['molecule_id'] for v in inventory.values())
 for (mid,variant,f),g in prim.groupby(['molecule_id','variant','descriptor']):
  expected=counts[mid]
  if len(g)!=expected:raise ValueError(f'Incomplete reference ensemble {mid}/{variant}/{f}: {len(g)}/{expected}')
  prop=FIELDS[f][2];values=g.reference.to_numpy();sr=pub[mid].get(prop,{})
  for red in ['min','max','delta']:
   if sr.get(red) is None:continue
   val=float(values.min() if red=='min' else values.max() if red=='max' else values.max()-values.min())
   aggregates.append({'molecule_id':int(mid),'variant':variant,'descriptor':f,'reduction':red,'metric':f+'_'+red,'reference':float(sr[red]),'independent':val,'residual':val-float(sr[red]),'absolute_error':abs(val-float(sr[red])),'conformers':len(g)})
 agg=pd.DataFrame(aggregates);agg.to_csv(OUT/'primary_ensemble_vs_published.csv',index=False)
 summary={f'{v}/{m}':finite_metrics(g.to_dict('records'),y='independent') for (v,m),g in agg.groupby(['variant','metric'])};save(OUT/'primary_ensemble_metrics.json',summary)
 # Complete diagnostic ensembles were requested for the twenty headline outliers.
 top=json.loads((ROOT/'analysis/headline_top20.json').read_text());dossiers=[]
 for row in top:
  mid=row['molecule_id'];entries=[r for r in diagnostic.values() if int(r['id'].split(':')[1])==mid];variants={}
  for v in ['matched_sut_001','matched_sut_0001','primary_001','primary_0001']:
   good=[(r['id'],r[v]['max_delta_qvbur']) for r in entries if v in r and 'error' not in r[v]]
   if len(good)!=row['conformers']:variants[v]={'error':'Incomplete independent ensemble','available':len(good),'expected':row['conformers']};continue
   case,val=min(good,key=lambda z:z[1]);variants[v]={'value':val,'extremizer_id':case,'residual_vs_published':val-row['reference']}
  center_angles=[r['center_angle_degrees'] for r in entries if 'center_angle_degrees' in r]
  dossiers.append({**row,'reference_definition':'Original Kraken DFT: raw P−bonded-neighbor vector sum, normalized to2.28Å; Bondi1.17 without H; Morfeus density0.001; maximum adjacent quadrant difference across3substituent frames, then ensemble minimum.','inputs':[{'id':r['id'],'SDF':f"primary/conformers/{r['id'].split(':')[2]}.body",'SDF_sha256':sha(ROOT/'primary/conformers'/f"{r['id'].split(':')[2]}.body"),'primary_neighbors':r['primary_neighbors'],'SDF_neighbors':r['topological_neighbors'],'primary_donor':r['primary_donor'],'input_donor':r['input_donor'],'primary_center':r['primary_center'],'sut_center':r['sut_center']} for r in entries],'reference_variants':variants,'center_angle_degrees_range':[min(center_angles),max(center_angles)] if center_angles else None,'likely_cause':'Documented raw-bond-vector versus normalized-bond-vector center convention plus default density0.01 versus primary0.001; exact causal contributions given by untuned reference variants. Residual after matching primary convention remains subject to SDF rounding/current reference version.','scientific_magnitude_A3':row['absolute_error'],'affected_claims':['K02','K06','K07'],'proposed_correction':'Correct electronic-LMO attribution in documentation; explicitly choose/document whether to reproduce Kraken raw-vector center and density. No production changes made; any future convention change requires separate scientific review.'})
 save(OUT/'headline_outlier_dossiers.json',dossiers)
 # All descriptor top20 entries are retained with their selected-conformer independent diagnostics.
 all_top=json.loads((ROOT/'analysis/top20_by_descriptor.json').read_text());all_dossiers={}
 for metric,rows in all_top.items():
  family=rows[0]['descriptor'];records=[]
  for row in rows:
   cid=row.get('sut_extremizer_conformer_id');key=f"KRAKEN:{row['molecule_id']}:{int(cid)}" if cid is not None else None
   checks={k:v for k,v in diagnostic.get(key,{}).items() if k not in ['id']};independent_aggregate=[v for v in aggregates if v['molecule_id']==row['molecule_id'] and v['metric']==metric]
   records.append({**row,'input_id':key,'input_path':f'primary/conformers/{int(cid)}.body' if cid is not None else 'Full ensemble in prepared_all/inventory.jsonl','reference_definition':FIELDS[family][2]+' / '+row['reduction']+'; original source and REPORT.md convention table','selected_conformer_diagnostics':checks,'independent_complete_ensemble':independent_aggregate,'likely_cause':('Radius1.20 versus1.09H, raw versus unit-vector center, angular grid forB1, quaternion approximation forL/B5; variants separate these effects.' if family.startswith('sterimol') else 'Neighbor definition, coordinate precision and scalar-triple/sign convention; exact frozen independent results retained.' if family.startswith('pyr') else 'Center definition, finite grid density and frame/order; reductions can select different conformers. Delta outliers require both extrema; all underlying measurements are retained.'),'scientific_magnitude':row['absolute_error'],'affected_claims':['K04'] if family.startswith('sterimol') else ['K05'] if family.startswith('pyr') else ['K02','K06','K07'],'proposed_correction':'Specify published versus StericX-specific conventions and quantify deviations; inspect exact retained geometry before any production correction.','limitation':'Selected SUT extremizer is diagnostic, not necessarily the primary extremizer. Complete primary ensemble references are provided for Sterimol/pyramidalization; complete independent BV ensembles are provided for headline20 only.'})
  all_dossiers[metric]=records
 save(OUT/'all_top20_outlier_dossiers.json',all_dossiers)
 save(OUT/'manifest.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'script_sha256':sha(Path(__file__)),'source_manifests':{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'sut/manifest_frozen.json',ROOT/'morfeus_outliers/manifest_frozen.json',ROOT/'full_analytic_reference/manifest_frozen.json',ROOT/'analysis/manifest.json']},'outputs':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()},'no_production_change':True})
 print(json.dumps({'analytic':json.loads((OUT/'same_center_analytic_metrics.json').read_text()),'headline_primary_max_residual':max(abs(d['reference_variants']['primary_0001']['residual_vs_published']) for d in dossiers)},indent=2))
if __name__=='__main__':main()
