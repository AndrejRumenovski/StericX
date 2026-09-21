#!/usr/bin/env python3
"""Frozen diagnostic selection and explicit, untuned Morfeus convention matrix."""
from __future__ import annotations
import ast, concurrent.futures, datetime, hashlib, json, math, os
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from morfeus import BuriedVolume, Sterimol, Pyramidalization
from morfeus.utils import get_radii
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'morfeus_outliers'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n')
def emit(f,v):f.write(json.dumps(v,separators=(',',':'),allow_nan=False)+'\n');f.flush()
# Only the rcov dictionary and get_conmat function are extracted; no original workflow executes.
source=ast.parse((ROOT/'raw_sources/official_PL_dft_library_201027.body').read_text())
selected=[n for n in source.body if (isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='rcov' for t in n.targets)) or (isinstance(n,ast.FunctionDef) and n.name=='get_conmat')]
np.math=math # Original np.math.exp alias was removed in NumPy 2; numerical function is unchanged.
namespace={'np':np};exec(compile(ast.Module(body=selected,type_ignores=[]),'<frozen primary get_conmat>','exec'),namespace)
get_conmat=namespace['get_conmat']
def protect(fun,*a,**kw):
 try:return fun(*a,**kw)
 except Exception as e:return {'error':type(e).__name__+': '+str(e)}
def bv(elements,coordinates,donor,neighbors,center,density,radii=None):
 coords=np.vstack((coordinates,center));elems=elements+['Pd'];metal=len(elems)
 rs=None if radii is None else list(radii)+[1.63]
 orientations=[]
 for neighbor in neighbors:
  obj=BuriedVolume(elems,coords.copy(),metal,excluded_atoms=[metal],z_axis_atoms=[donor+1],xz_plane_atoms=[neighbor+1],density=density,radius=3.5,radii=rs,radii_type='bondi',radii_scale=1.17,include_hs=False)
  obj.octant_analysis();q=list(obj.quadrants['buried_volume'].values());o=list(obj.octants['buried_volume'].values())
  orientations.append({'plane_atom':neighbor,'buried_volume':float(obj.buried_volume),'qvbur_min':float(min(q)),'qvbur_max':float(max(q)),'max_delta_qvbur':float(max(abs(q[j]-q[j-1]) for j in range(4))),'ovbur_min':float(min(o)),'ovbur_max':float(max(o)),'near_vbur':float(sum(o[4:])),'far_vbur':float(sum(o[:4]))})
 result={k:(min(o[k] for o in orientations) if k in ['qvbur_min','ovbur_min'] else max(o[k] for o in orientations)) for k in ['qvbur_min','qvbur_max','max_delta_qvbur','ovbur_min','ovbur_max']}
 result.update({k:orientations[-1][k] for k in ['buried_volume','near_vbur','far_vbur']})
 result['orientations']=orientations
 return result
def ster(elements,coords,donor,center,radii):
 obj=Sterimol(elements+['Pd'],np.vstack((coords,center)),len(elements)+1,donor+1,radii=list(radii)+[1.63],n_rot_vectors=3600)
 return {'sterimol_l':float(obj.L_value),'sterimol_b1':float(obj.B_1_value),'sterimol_b5':float(obj.B_5_value)}
def pyr(elements,coords,donor,neighbors):
 o=Pyramidalization(coords,donor+1,elements=elements,neighbor_indices=[i+1 for i in neighbors]);return {'pyr_p':float(o.P),'pyr_alpha':float(o.alpha)}
def run(record):
 req=record['request'];sut=record['sut'];d=req['donor'];elements=[a['element'] for a in req['atoms']];coords=np.array([a['position'] for a in req['atoms']]);f32=np.array([a['position'] for a in sut['atoms']]);radii=[a['vdw_radius'] for a in sut['atoms']]
 connectivity=get_conmat(elements,coords);neighbors=np.flatnonzero(connectivity[d]).tolist();topology=req['neighbors']
 vec=sum((coords[d]-coords[i] for i in neighbors),start=np.zeros(3));center=coords[d]+2.28*vec/np.linalg.norm(vec)
 result={'id':req['id'],'primary_neighbors':neighbors,'topological_neighbors':topology,'primary_center':center.tolist(),'sut_center':sut['center'],'primary_donor':next((i for i,e in enumerate(elements) if e=='P' and sum(connectivity[i])<=3),None),'input_donor':d,'atom_mapping':'Unchanged SDF atom order; appended Pd dummy excluded from BV; all indices zero-based in this record.'}
 center_sut=np.array(sut['center']) if isinstance(sut['center'],list) else None
 if center_sut is not None:
  u=(center_sut-f32[d]);v=(center-coords[d]);result['center_shift_A']=float(np.linalg.norm(center-center_sut));result['center_angle_degrees']=float(np.degrees(np.arccos(np.clip(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)),-1,1))))
  order=[req['reference']]+[r['index'] for r in sut['bonded_neighbors'] if r['index']!=req['reference']]
  result['matched_sut_001']=protect(bv,elements,f32,d,order,center_sut,.01,radii)
  result['matched_sut_0001']=protect(bv,elements,f32,d,order,center_sut,.001,radii)
 result['primary_001']=protect(bv,elements,coords,d,neighbors,center,.01)
 result['primary_0001']=protect(bv,elements,coords,d,neighbors,center,.001)
 bondi=np.array(get_radii(elements,radii_type='bondi'));paton=bondi.copy();paton[paton==1.2]=1.09
 result['ster_primary_center_bondi']=protect(ster,elements,coords,d,center,bondi)
 result['ster_primary_center_paton']=protect(ster,elements,coords,d,center,paton)
 if center_sut is not None:
  result['ster_sut_center_sut_radii']=protect(ster,elements,f32,d,center_sut,radii)
  result['ster_sut_center_paton']=protect(ster,elements,f32,d,center_sut,paton)
 result['pyr_f32_topology']=protect(pyr,elements,f32,d,topology)
 result['pyr_f64_topology']=protect(pyr,elements,coords,d,topology)
 result['pyr_f64_primary_neighbors']=protect(pyr,elements,coords,d,neighbors)
 return result

def main():
 OUT.mkdir(exist_ok=False)
 top=json.loads((ROOT/'analysis/top20_by_descriptor.json').read_text());headline={r['molecule_id'] for r in top['max_delta_qvbur_min']}
 single={f"KRAKEN:{r['molecule_id']}:{r['sut_extremizer_conformer_id']}" for name,rows in top.items() if not name.endswith('_delta') for r in rows if r['sut_extremizer_conformer_id'] is not None}
 # pandas nullable integer metadata is represented as a float; normalize only IDs.
 single={':'.join(k.split(':')[:2]+[str(int(float(k.split(':')[2])))]) for k in single}
 errors={r['id'] for r in json.loads((ROOT/'analysis/operation_errors.json').read_text())}
 requests={r['id']:r for r in map(json.loads,(ROOT/'prepared_all/requests.jsonl').open()) if int(r['id'].split(':')[1]) in headline or r['id'] in single or r['id'] in errors}
 selected=[]
 with (ROOT/'sut/stdout.jsonl').open() as f,(OUT/'selected_inputs.jsonl').open('w') as dst:
  for line in f:
   row=json.loads(line)
   if row['id'] in requests:
    pair={'request':requests[row['id']],'sut':row};selected.append(pair);emit(dst,pair)
 save(OUT/'selection_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'source_SUT_manifest_sha256':sha(ROOT/'sut/manifest_frozen.json'),'analysis_manifest_sha256':sha(ROOT/'analysis/manifest.json'),'selected_inputs_sha256':sha(OUT/'selected_inputs.jsonl'),'script_sha256':sha(Path(__file__)),'primary_code_sha256':sha(ROOT/'raw_sources/official_PL_dft_library_201027.body'),'N':len(selected),'headline_ligands_all_conformers':sorted(headline),'selection':'All conformers of the 20 largest absolute headline min(maxΔQvbur) outliers; SUT-selected conformers for all available per-descriptor top20 min/max/vburminconf entries; every SUT operation error. No result filtering.','variant_keys':{'matched_sut_001':'density0.01; actual SUT f32 coordinates, radii and center; SUT neighbor/reference order.','matched_sut_0001':'same, density0.001','primary_001':'SDF float64 coords, original primary get_conmat and raw-vector center; Bondi1.17, density0.01','primary_0001':'same, original published density0.001'},'reference_limit':'Morfeus0.8.0 is independently executed contemporary reference implementation, not physical truth and not frozen 2021 Morfeus software. Raw API SDF coordinates rounded to4decimals; original calculation coordinate precision unavailable.','numpy_compatibility':'np.math=math restores removed NumPy alias used by original get_conmat; no formula changed.'})
 print('Frozen selection',len(selected),flush=True)
 with (OUT/'reference_raw.jsonl').open('w') as dst,concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
  for i,row in enumerate(pool.map(run,selected),1):
   emit(dst,row)
   if i%25==0:print('Completed',i,'/',len(selected),flush=True)
 save(OUT/'manifest_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'reference_raw_sha256':sha(OUT/'reference_raw.jsonl'),'selection_manifest_sha256':sha(OUT/'selection_frozen.json'),'record_count':len(selected),'frozen_before_interpretation':True})
if __name__=='__main__':main()
