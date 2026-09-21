#!/usr/bin/env python3
"""Fresh full-corpus Morfeus Sterimol and pyramidalization under primary conventions."""
import concurrent.futures, datetime, json, math
from pathlib import Path
import numpy as np
from kraken_morfeus_outliers import ROOT, namespace, sha, save, emit, protect, ster, pyr, get_radii
OUT=ROOT/'full_analytic_reference'
def run(req):
 d=req['donor'];elems=[a['element'] for a in req['atoms']];coords=np.array([a['position'] for a in req['atoms']]);rcov=namespace['rcov'];neighbors=[]
 # Only donor row of the original connection matrix is needed; identical scalar operations.
 for i,e in enumerate(elems):
  if i!=d and e in rcov:
   distance=np.linalg.norm(coords[i]-coords[d]);ratio=((rcov[elems[d]]+rcov[e])*(4.0/3.0))/distance
   if 1.0/(1.0+math.exp(-16.0*(ratio-1.0)))>.85:neighbors.append(i)
 vec=np.zeros(3)
 for i in neighbors:vec+=coords[d]-coords[i]
 center=coords[d]+2.28*vec/np.linalg.norm(vec);radii=np.array(get_radii(elems,radii_type='bondi'));radii[radii==1.2]=1.09
 # Original Morfeus Pyramidalization uses nearest3 (excluding dummy), not explicit graph.
 nearest=np.argsort(np.linalg.norm(coords-coords[d],axis=1))[1:4].tolist()
 return {'id':req['id'],'primary_neighbors':neighbors,'topological_neighbors':req['neighbors'],'pyr_nearest_neighbors':nearest,'primary_center':center.tolist(),'sterimol_primary':protect(ster,elems,coords,d,center,radii),'pyr_primary_nearest':protect(pyr,elems,coords,d,nearest),'pyr_topology':protect(pyr,elems,coords,d,req['neighbors'])}
def main():
 OUT.mkdir(exist_ok=False);source=ROOT/'prepared_all/requests.jsonl';records=list(map(json.loads,source.open()))
 save(OUT/'plan_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'input_sha256':sha(source),'script_sha256':sha(Path(__file__)),'helper_sha256':sha(Path(__file__).with_name('kraken_morfeus_outliers.py')),'N':len(records),'selection':'All31,721 available DFT conformers, no residual-dependent exclusions.','convention':'Primary raw-vector Pd center; scalar donor-row connectivity from frozen Kraken get_conmat; Morfeus0.8 Sterimol with Bondi H1.09,n3600,+.40L; Pyr nearest3 matching original default and independently explicit SDF topology.','reference_limit':'Current Morfeus version and4-decimal SDF export do not exactly reconstruct original calculation dependency/coordinate precision.'})
 with (OUT/'reference_raw.jsonl').open('w') as dst,concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
  for i,row in enumerate(pool.map(run,records,chunksize=16),1):
   emit(dst,row)
   if i%1000==0:print('Completed',i,'/',len(records),flush=True)
 save(OUT/'manifest_frozen.json',{'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'raw_sha256':sha(OUT/'reference_raw.jsonl'),'plan_sha256':sha(OUT/'plan_frozen.json'),'N':len(records),'frozen_before_interpretation':True})
if __name__=='__main__':main()
