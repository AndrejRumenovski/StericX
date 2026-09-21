#!/usr/bin/env python3
"""Independent checks of actual largest axis errors and published-geometry mismatch."""
import datetime,itertools,json
import numpy as np
import pandas as pd
from kraken_morfeus_outliers import ROOT,save,sha,ster,pyr
OUT=ROOT/'focused_diagnosis';OUT.mkdir(exist_ok=False)
a=json.loads((ROOT/'interpretation/same_center_analytic_metrics.json').read_text());axis_ids={v['maximum_case']['id'] for v in a.values()}
errors=json.loads((ROOT/'analysis/operation_errors.json').read_text());error_ids={r['id'] for r in errors}
agg=pd.read_csv(ROOT/'interpretation/primary_ensemble_vs_published.csv');chosen=agg[(agg.variant=='pyr_nearest') & agg.reduction.isin(['min','max'])].nlargest(20,'absolute_error');mids=set(chosen.molecule_id)
requests={r['id']:r for r in map(json.loads,(ROOT/'prepared_all/requests.jsonl').open()) if r['id'] in axis_ids|error_ids or int(r['id'].split(':')[1]) in mids}
sut={r['id']:r for r in map(json.loads,(ROOT/'sut/stdout.jsonl').open()) if r['id'] in requests}
full={r['id']:r for r in map(json.loads,(ROOT/'full_analytic_reference/reference_raw.jsonl').open()) if r['id'] in requests}
save(OUT/'plan_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'script_sha256':sha(__import__('pathlib').Path(__file__)),'axis_ids':sorted(axis_ids),'pyramidalization_ligands':sorted(map(int,mids)),'method':'Axis: direct Morfeus with identical observed f32 coordinates/center/radii. Pyramidalization: exhaust4096corner perturbations of donor and3neighbors at±0.00005Å (4decimal export box); empirical envelope only, not rigorous interior extremum bound. No SUT calls or tuned thresholds.','seed':None,'sut_sha256':sha(ROOT/'sut/stdout.jsonl')})
axis=[];pyrrows=[];bond=[]
for key,r in requests.items():
 s=sut[key];e=[v['element'] for v in r['atoms']];c=np.array([v['position'] for v in r['atoms']]);d=r['donor']
 if key in axis_ids:
  f=np.array([v['position'] for v in s['atoms']]);rs=[v['vdw_radius'] for v in s['atoms']];val=ster(e,f,d,np.array(s['center']),rs)
  axis.append({'id':key,'request':r,'sut':s['sterimol_coordination'],'independent_morfeus':val,'independent_projection':{k:v['maximum_case'] for k,v in a.items() if v['maximum_case']['id']==key},'definition':'Axial support=max_i((x_i−center)·unit_axis+r_i)+0.40; B5=max_i(norm((x_i−center)−axial_i unit_axis)+r_i).','cause':'glam0.30.10 rotation-arc approximates nearparallel or antiparallel axes using ±1 threshold; actual nonzero tilt is discarded. Direct Morfeus and dot projection do not share this branch.','proposal':'Independently justify exact frame construction or direct projections in a future production change; retain these original geometries as regression witnesses.'})
 if key in error_ids:
  neigh=s['bonded_neighbors'];bond.append({'id':key,'SDF_neighbors':r['neighbors'],'primary_neighbors':full[key]['primary_neighbors'],'sut_neighbors':neigh,'bond_geometry':[{'index':v['index'],'element':e[v['index']],'distance':float(np.linalg.norm(c[v['index']]-c[d])),'sut_threshold':1.3*(s['atoms'][d]['covalent_radius']+s['atoms'][v['index']]['covalent_radius']),'in_SDF_graph':v['index'] in r['neighbors']} for v in neigh],'sut_errors':{op:s[op]['error'] for op in ['buried_volume','sterimol_coordination'] if 'error' in s[op]},'independent':full[key],'proposal':'Use explicit topology where known; distance heuristic has no universal chemical guarantee. A symmetric nonzero occupied sphere is scientifically valid and should not alone be labeled a degenerate frame; investigate separately before a correction.'})
 if int(key.split(':')[1]) in mids:
  neighbors=r['neighbors'];small=np.array([c[d]]+[c[i] for i in neighbors]);vals=[]
  for signs in itertools.product([-1,1],repeat=12):
   x=small+np.array(signs).reshape(4,3)*.00005;v=pyr(['P']+[e[i] for i in neighbors],x,0,[1,2,3]);vals.append([v['pyr_p'],v['pyr_alpha']])
  v=np.asarray(vals);pyrrows.append({'id':key,'unperturbed':full[key]['pyr_primary_nearest'],'sample_count':len(v),'pyr_p_corner_range':[float(v[:,0].min()),float(v[:,0].max())],'pyr_alpha_corner_range':[float(v[:,1].min()),float(v[:,1].max())],'max_corner_deviation_alpha':float(np.max(np.abs(v[:,1]-full[key]['pyr_primary_nearest']['pyr_alpha']))),'max_corner_deviation_P':float(np.max(np.abs(v[:,0]-full[key]['pyr_primary_nearest']['pyr_p'])))})
save(OUT/'axis_counterexamples.json',axis);save(OUT/'operation_failure_dossiers.json',bond);save(OUT/'pyramidalization_coordinate_rounding.json',pyrrows)
save(OUT/'manifest.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'files':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()},'no_production_change':True})
print('Done',len(axis),len(bond),len(pyrrows),flush=True)
