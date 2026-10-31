#!/usr/bin/env python3
"""Independent DBSTEP cube parsing, density occupation, public volume and convergence."""
import argparse
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime, timezone
import io
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
import numpy as np
from scipy.spatial import cKDTree
from dbstep import Dbstep, parse_data, sterics, constants
from study_conformer_ensembles import sha

ROOT=Path(__file__).resolve().parents[1]
BOHR=0.529177210544
EPS=float(np.finfo(float).eps)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--engine',type=Path,default=ROOT/'target/debug/examples/density_report')
    args=parser.parse_args(); out=args.output.resolve(); out.mkdir(parents=True,exist_ok=False)
    inputs=out/'inputs';inputs.mkdir()
    official=ROOT/'.stericx/paton/reference/DBSTEP'
    shutil.copyfile(official/'LICENSE.txt',inputs/'DBSTEP_LICENSE.txt')
    cases=[]
    for name in ['Ne_medium.cube','benzene_coarse.cube','CH2CMe3_100.cube']:
        path=inputs/name;shutil.copyfile(official/'tests/cube_files'/name,path)
        for radius in [1.,1.5,2.]: cases.append({'path':path,'radius':radius,'threshold':.0016,'kind':'official_density'})
    for spacing in [.25,.125,.0625]:
        n=round(4/spacing)+1
        path=inputs/f'gaussian_{spacing}.cube'
        axis=-2.+np.arange(n)*spacing
        x,y,z=np.meshgrid(axis,axis,axis,indexing='ij')
        density=np.exp(-(x*x+y*y+z*z)).ravel()
        with path.open('w') as stream:
            stream.write(f'analytic isotropic Gaussian\ntotal electron density model\n1 -2 -2 -2\n{n} {spacing} 0 0\n{n} 0 {spacing} 0\n{n} 0 0 {spacing}\n1 1 0 0 0\n')
            for i in range(0,len(density),6): stream.write(' '.join(f'{v:.17e}' for v in density[i:i+6])+'\n')
        cases.append({'path':path,'radius':1.8*BOHR,'threshold':math.exp(-1.),'kind':'analytic_gaussian','spacing_bohr':spacing})
    requests=[];parsed=[]
    for case in cases:
        cube=parse_data.CubeParser(str(case['path']),'cube');parsed.append(cube)
        center=np.asarray(cube.CARTESIANS[0])/constants.BOHR_TO_ANG*BOHR
        requests.append({'path':str(case['path']),'declaration':{'field':'total_electron_density','unit':'electrons_per_bohr_cubed',
            'provenance':f"{case['kind']}; SHA256={sha(case['path'])}; original quantum method not inferred"},
            'center_angstrom':center.tolist(),'radius_angstrom':case['radius'],'threshold_electrons_per_bohr3':case['threshold']})
    (out/'requests.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in requests))
    receipt={'created_utc':datetime.now(timezone.utc).isoformat(),'command':sys.argv,'controller_sha256':sha(__file__),
        'engine':str(args.engine.resolve()),'engine_sha256':sha(args.engine),
        'sources':{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/lib.rs',ROOT/'examples/density_report.rs',*sorted((ROOT/'src/density').glob('*.rs'))]},
        'protocol_sha256':sha(ROOT/'docs/methods/electron_density_buried_volume.md'),
        'reference_environment_sha256':sha(ROOT/'docs/validation/dbstep/environment.json'),
        'fixtures':{p.name:sha(p) for p in sorted(inputs.iterdir())}}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    process=subprocess.run([str(args.engine.resolve())],input=(out/'requests.jsonl').read_text(),text=True,capture_output=True)
    (out/'native_stdout.jsonl').write_text(process.stdout);(out/'native_stderr.txt').write_text(process.stderr)
    observed=[json.loads(line) for line in process.stdout.splitlines()]
    if process.returncode or len(observed)!=len(cases): raise RuntimeError('incomplete native density capture')
    rows=[];failures=[];convergence=[]
    for case,request,cube,response in zip(cases,requests,parsed,observed):
        identity=f"{case['path'].name}/R={case['radius']}"
        print(identity,flush=True)
        stdout,stderr=io.StringIO(),io.StringIO()
        try:
            if response['status']!='success': raise ValueError(response)
            native=response['value']; meta=native['grid']; dims=[cube.xdim,cube.ydim,cube.zdim]
            if meta['dimensions']!=dims: raise ValueError('parser dimensions differ')
            for key,legacy in [('origin_angstrom',cube.ORIGIN),('vectors_angstrom',cube.INCREMENTS)]:
                common=np.asarray(legacy)/constants.BOHR_TO_ANG*BOHR
                if not np.allclose(meta[key],common,rtol=64*EPS,atol=64*EPS): raise ValueError('parser geometry differs')
            old_origin=np.asarray(cube.ORIGIN)-np.asarray(cube.CARTESIANS[0])
            old_axes=[np.linspace(old_origin[i],old_origin[i]+(dims[i]-1)*cube.SPACING,dims[i]) for i in range(3)]
            old_grid=np.array(np.meshgrid(*old_axes,indexing='ij')).reshape(3,-1).T
            radius=case['radius'];threshold=case['threshold']
            margin=min(min(-a[0],a[-1]) for a in old_axes)
            if margin<1.1*radius: raise ValueError('public DBSTEP would resize grid; comparison unavailable')
            old_options=SimpleNamespace(grid=cube.SPACING,radius=radius,isoval=threshold,verbose=False,debug=False)
            occupied,whole_volume,tree=sterics.occupied_dens(old_grid,cube.DENSITY,old_options)
            fraction,_=sterics.buried_vol(occupied,tree,np.zeros(3),radius,0.,old_options)
            n_sphere=len(tree.query_ball_point(np.zeros(3),radius));n_occupied=int(np.count_nonzero(np.sum(occupied**2,axis=1)<=radius**2))
            with redirect_stdout(stdout),redirect_stderr(stderr):
                public=Dbstep.dbstep(str(case['path']),surface='density',volume=True,atom1=1,r=radius,isoval=threshold,quiet=True,norot=True)
            # A second official support-function lane uses the current CODATA unit convention.
            new_origin=np.asarray(meta['origin_angstrom'])-np.asarray(request['center_angstrom'])
            new_spacing=meta['vectors_angstrom'][0][0]
            new_axes=[new_origin[i]+np.arange(dims[i])*new_spacing for i in range(3)]
            new_grid=np.array(np.meshgrid(*new_axes,indexing='ij')).reshape(3,-1).T
            opt=SimpleNamespace(grid=new_spacing,radius=radius,isoval=threshold,verbose=False,debug=False)
            occ,_,new_tree=sterics.occupied_dens(new_grid,cube.DENSITY,opt)
            common_percent,_=sterics.buried_vol(occ,new_tree,np.zeros(3),radius,0.,opt)
            ns=len(new_tree.query_ball_point(np.zeros(3),radius));no=int(np.count_nonzero(np.sum(occ**2,axis=1)<=radius**2))
            corrected_volume=n_occupied*new_spacing**3
            checks={'public_support_percent':abs(float(public.bur_vol)-fraction)<=64*EPS*max(1.,abs(fraction)),
                'common_sphere_count':native['sphere_nodes']==ns,'common_occupied_count':native['occupied_nodes']==no,
                'public_sphere_count':native['sphere_nodes']==n_sphere,'public_occupied_count':native['occupied_nodes']==n_occupied,
                'common_percent':abs(native['percent_buried_volume']-common_percent)<=64*EPS*max(1.,abs(common_percent)),
                'public_percent':abs(native['percent_buried_volume']-float(public.bur_vol))<=64*EPS*max(1.,abs(float(public.bur_vol))),
                'unit_corrected_volume':abs(native['buried_volume_angstrom3']-corrected_volume)<=64*EPS*max(1.,abs(corrected_volume))}
            checks={key:bool(value) for key,value in checks.items()}
            if not all(checks.values()): failures.append({'case':identity,'checks':checks})
            row={'case':identity,'kind':case['kind'],'native':native,'public_percent':float(public.bur_vol),
                'support_percent':fraction,'public_sphere_nodes':n_sphere,'public_occupied_nodes':n_occupied,
                'common_sphere_nodes':ns,'common_occupied_nodes':no,'common_percent':common_percent,
                'reference_bohr_angstrom':constants.BOHR_TO_ANG,'reference_spacing_angstrom':cube.SPACING,
                'reference_grid_buried_volume':n_occupied*cube.SPACING**3,
                'reference_whole_molecule_grid_volume_NOT_buried':float(public.occ_vol),
                'unit_corrected_grid_buried_volume':corrected_volume,'domain_margin_angstrom':margin,
                'checks':checks,'stdout':stdout.getvalue(),'stderr':stderr.getvalue()}
            rows.append(row)
            if case['kind']=='analytic_gaussian':
                true_volume=4*math.pi*BOHR**3/3
                delta=math.sqrt(3)*case['spacing_bohr']*BOHR/2
                bound=4*math.pi*((BOHR+delta)**3-max(0.,BOHR-delta)**3)/3
                error=abs(native['buried_volume_angstrom3']-true_volume)
                convergence.append({'spacing_bohr':case['spacing_bohr'],'native':native['buried_volume_angstrom3'],
                    'analytic':true_volume,'absolute_error':error,'boundary_layer_bound':bound,'passed':error<=bound})
        except (Exception,SystemExit) as error:
            failures.append({'case':identity,'error':repr(error),'stdout':stdout.getvalue(),'stderr':stderr.getvalue()})
    if not all(r['passed'] for r in convergence):failures.append({'error':'analytic boundary-layer gate failed'})
    summary={'cases':len(cases),'completed':len(rows),'failures':failures,'convergence':convergence,
             'admitted':len(rows)==len(cases) and not failures,'scope':'declared scalar cubic grids, covered spheres, explicit threshold; no physical threshold validation'}
    (out/'reference_raw.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (out/'manifest.json').write_text(json.dumps({str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return not summary['admitted']

if __name__=='__main__': raise SystemExit(main())
