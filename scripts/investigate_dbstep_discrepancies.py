#!/usr/bin/env python3
"""Triangulate every volume case and resolve the frozen sphere-boundary witness."""
import argparse
from contextlib import redirect_stdout, redirect_stderr
from decimal import Decimal, localcontext
import io
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
from dbstep import sterics
from morfeus import BuriedVolume
from types import SimpleNamespace
from study_conformer_ensembles import f32, sha, write_csv, references
from validate_dbstep_volume import load_cases, grid_key

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
    source=ROOT/'docs/validation/dbstep/volume_v2'
    for name,digest in json.loads((source/'manifest.json').read_text()).items():
        if sha(source/name)!=digest: raise ValueError(f'changed evidence: {name}')
    records=[json.loads(s) for s in (source/'native_stdout.jsonl').read_text().splitlines()]
    cases=load_cases(); grids={}; rows=[]; failures=[]; boundary=None
    stdout,stderr=io.StringIO(),io.StringIO()
    for case,record in zip(cases,records):
        key=grid_key(case['options']['volume'])
        if record['dump']['points'] is not None: grids[key]=f32(record['dump']['points'])
        with redirect_stdout(stdout),redirect_stderr(stderr):
            values,errors,frames=references(case['row'],case['options'])
        if errors: failures.append({'case_id':case['id'],'errors':errors})
        rows.append({'case_id':case['id'],'native':float(f32(record['buried_volume']['buried_volume'])),
                     'morfeus':values['buried_volume'],'morfeus_frames':frames,
                     'classification':'numerical integration difference',
                     'note':'independently aligned and generated Morfeus projection grid; f64 radii scaling'})
        if case['id'].endswith('sphere_radius=2.5'):
            boundary=(case,record,grids[key])
    case,record,points=boundary
    radius=float(f32(case['options']['volume']['sphere_radius'])); r2=radius**2
    d2=np.sum(points**2,axis=1); outside=np.flatnonzero(d2>r2)
    exact=[]
    with localcontext() as ctx:
        ctx.prec=80
        for i in outside:
            ds=sum(Decimal.from_float(float(x))**2 for x in points[i])
            exact.append({'point_index':int(i),'point':points[i].tolist(),
                          'exact_squared_distance':str(ds),'radius_squared':str(Decimal.from_float(r2)),
                          'excess_squared_distance':str(ds-Decimal.from_float(r2))})
    retained=[at for at in record['atoms'] if at['element'].upper()!='H']
    scaled=np.array([np.float32(at['vdw_radius'])*np.float32(case['options']['volume']['radii_scale']) for at in retained])
    options=SimpleNamespace(grid=1.0,radius=radius,verbose=False,debug=False)
    masks=[]; requests=[]
    for frame in record['dump']['orientations']:
        coords=f32([at['position'] for at in frame['aligned_atoms']])
        occupied,tree,_=sterics.occupied(points,coords,scaled.astype(float),np.zeros(3),options)
        db_mask=np.zeros(len(points),dtype=bool)
        for xyz,r in zip(coords,scaled): db_mask[tree.query_ball_point(xyz,float(r))]=True
        delta=points.astype(np.float32)[:,None,:]-coords.astype(np.float32)[None,:,:]
        sq=delta*delta
        ds=(sq[:,:,0]+sq[:,:,1])+sq[:,:,2]
        native_model=np.any(ds<=np.array([at['radius_squared'] for at in frame['aligned_atoms']],dtype=np.float32),axis=1)
        vs=np.float32(np.float32(np.float32(4)*np.float32(math.pi))*np.float32(np.float32(radius)*np.float32(radius)*np.float32(radius))/np.float32(3))
        reconstructed=np.float32(np.float32(native_model.sum())/np.float32(len(points)))*vs
        if reconstructed.tobytes()!=np.float32(frame['buried_volume']).tobytes():
            failures.append({'case_id':case['id'],'plane':frame['plane'],'error':'IEEE count model does not reproduce native volume'})
        for i in outside:
            requests.append({'id':f"{frame['plane']}/{int(i)}",'op':'raw_occupancy','sphere_radius':radius,
                             'points':[points[i].tolist()],'atoms':frame['aligned_atoms']})
        b=BuriedVolume(['Ni',*[at['element'] for at in retained]],np.vstack([np.zeros(3),coords]),1,
                       radii=np.r_[0.,scaled.astype(float)],include_hs=True,radius=radius,density=float(f32(case['options']['volume']['density'])))
        masks.append({'plane':frame['plane'],'dbstep_full_grid_occupied':int(db_mask.sum()),
                      'dbstep_sphere_occupied':int((db_mask&(d2<=r2)).sum()),
                      'native_ieee_model_occupied':int(native_model.sum()),
                      'atom_occupancy_disagreement_indices':np.flatnonzero(db_mask!=native_model).tolist(),
                      'outside_occupied_indices':outside[db_mask[outside]].tolist(),
                      'native_public_volume':float(f32(frame['buried_volume'])),
                      'morfeus_sphere_count':len(b._sphere.points),'morfeus_occupied_count':len(b._buried_points),
                      'morfeus_volume':float(b.buried_volume)})
    build=json.loads((ROOT/'.stericx/paton/dbstep_build_v2/manifest.json').read_text())
    executable=build['binaries']['observer']
    if sha(executable['path'])!=executable['sha256']: raise ValueError('changed observer')
    raw=''.join(json.dumps(r)+'\n' for r in requests)
    (a.output/'boundary_requests.jsonl').write_text(raw)
    run=subprocess.run([executable['path']],input=raw,text=True,capture_output=True)
    (a.output/'boundary_native_stdout.jsonl').write_text(run.stdout)
    (a.output/'boundary_native_stderr.txt').write_text(run.stderr)
    outputs=[json.loads(s) for s in run.stdout.splitlines()]
    if run.returncode or len(outputs)!=len(requests): failures.append({'error':'incomplete native boundary observations'})
    for request,value in zip(requests,outputs):
        plane,index=map(int,request['id'].split('/'))
        expected=index in next(m['outside_occupied_indices'] for m in masks if m['plane']==plane)
        if bool(value['buried_volume'])!=expected: failures.append({'point':request['id'],'error':'native boundary occupancy differs'})
    report={'case_id':case['id'],'native_sphere_count':len(points),'dbstep_sphere_count':len(points)-len(outside),
            'boundary_points':exact,'frames':masks,'classification':'numerical integration difference',
            'earliest_divergence':'f32 versus f64 evaluation-sphere membership before atom occupancy',
            'production_change_justified':False,'original_gate_remains_failed':True,'failures':failures}
    (a.output/'boundary_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    (a.output/'morfeus_all_cases.json').write_text(json.dumps(rows,indent=2)+'\n')
    (a.output/'reference_stdout.txt').write_text(stdout.getvalue()); (a.output/'reference_stderr.txt').write_text(stderr.getvalue())
    (a.output/'receipt.json').write_text(json.dumps({'command':sys.argv,'script_sha256':sha(__file__),
          'source_manifest_sha256':sha(source/'manifest.json'),'observer':executable,
          'reference_environment_sha256':sha(ROOT/'docs/validation/dbstep/environment.json'),
          'helpers':{name:sha(ROOT/'scripts'/name) for name in ['validate_dbstep_volume.py','study_conformer_ensembles.py']}},indent=2)+'\n')
    (a.output/'manifest.json').write_text(json.dumps({f.name:sha(f) for f in sorted(a.output.iterdir()) if f.is_file()},indent=2)+'\n')
    print(json.dumps({'cases':len(rows),'boundary_points':len(outside),'frames':masks,'failures':failures},indent=2))
    return bool(failures)

if __name__=='__main__': raise SystemExit(main())
