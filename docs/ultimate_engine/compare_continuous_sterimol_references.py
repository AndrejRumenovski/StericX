#!/usr/bin/env python3
"""Matched-radii Morfeus and DBSTEP comparisons to continuous support.

Finite scans are checked against their mathematical angular error envelope,
not treated as continuous exact answers. DBSTEP's absolute-axis L is retained.
"""
import argparse
import hashlib
import json
import math
import inspect
import shutil
from pathlib import Path
import numpy as np
from morfeus import Sterimol
from dbstep import calculator, sterics

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--validation',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output;out.mkdir(exist_ok=False)
shutil.copyfile(__file__,out/'controller.py')
(out/'morfeus_calculate_source.py').write_text(inspect.getsource(Sterimol.calculate))
(out/'dbstep_rotate_source.py').write_text(inspect.getsource(calculator.rotate_mol))
cases=json.loads((a.validation/'cases.json').read_text())
observed={r['id']:r for r in map(json.loads,(a.validation/'native.jsonl').read_text().splitlines())}
result={'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.validation/'cases.json',a.validation/'native.jsonl',Path(__file__)]},'rows':[],'failures':[]}
for case in cases:
    coordinates=np.array([case['dummy']]+[r['position'] for r in case['atoms']],dtype=float)
    radii=np.array([0.]+[r['radius'] for r in case['atoms']],dtype=float)
    n=observed[case['id']]['measurements'][0]['result']
    row={'id':case['id'],'continuous':n}
    try:
        m=Sterimol(['Ni']+['C']*len(case['atoms']),coordinates.copy(),dummy_index=1,attached_index=case['base']+2,radii=radii,n_rot_vectors=3600)
        relative=coordinates-coordinates[0]
        aligned=calculator.rotate_mol(relative.copy(),1,relative[case['base']+1].copy())
        length,b5,b1,_=sterics.get_classic_sterimol(aligned,radii,['Bq']+['C']*len(case['atoms']))
        axis=relative[case['base']+1]/np.linalg.norm(relative[case['base']+1])
        row['dbstep_with_upstream_rounded_alignment']={'l_absolute_axis':float(length),'b1':float(b1),'b5':float(b5)}
        # Upstream rotation rounds coordinates to eight decimals. Preserve that
        # lane above; isolate its descriptor kernel using an unrounded frame.
        seed=np.eye(3)[int(np.argmin(np.abs(axis)))]
        x=seed-np.dot(seed,axis)*axis;x/=np.linalg.norm(x)
        y=np.cross(axis,x)
        common=np.column_stack([relative@x,relative@y,relative@axis])
        length,b5,b1,_=sterics.get_classic_sterimol(common,radii,['Bq']+['C']*len(case['atoms']))
        transverse=relative-np.outer(relative@axis,axis)
        rho=float(np.max(np.linalg.norm(transverse,axis=1)))
        scale=max(float(np.max(np.linalg.norm(relative,axis=1)+radii)),1e-300)
        margin=128*np.finfo(float).eps*scale
        row['morfeus']={'l_raw':float(m.L_value_uncorrected),'b1':float(m.B_1_value),'b5':float(m.B_5_value)}
        row['dbstep']={'l_absolute_axis':float(length),'b1':float(b1),'b5':float(b5)}
        row['arithmetic_margin']=margin
        row['angular_error_bounds']={}
        for library,count in [('morfeus',3600),('dbstep',361)]:
            values=row[library];bound=rho*math.pi/(count-1)
            row['angular_error_bounds'][library]=bound
            assert -margin <= values['b1']-n['b1'] <= bound+margin,(library,'B1',values['b1']-n['b1'],bound,margin)
            if library=='morfeus':
                # Morfeus maximizes the same finite angular support grid for B5.
                # At the nearest angle to the winning disc's radial direction,
                # support is short by at most rho*(1-cos(pi/(N-1))). This is a
                # convention-specific discretization bound, not a looser float
                # tolerance. Preserve the original invalid exact-B5 comparison.
                b5_bound=rho*2*math.sin(math.pi/(2*(count-1)))**2
                row['angular_error_bounds']['morfeus_b5']=b5_bound
                assert -margin <= n['b5']-values['b5'] <= b5_bound+margin,('Morfeus sampled B5',values['b5'],n['b5'],b5_bound,margin)
            else:
                assert abs(values['b5']-n['b5'])<=margin,(library,'B5',values['b5'],n['b5'],margin)
        assert abs(row['morfeus']['l_raw']-n['l'])<=margin,('Morfeus L',row['morfeus']['l_raw'],n['l'],margin)
        row['passed']=True
    except Exception as e:
        row['passed']=False;row['error']=repr(e);result['failures'].append(case['id'])
    result['rows'].append(row)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(case['id'],row['passed'],row.get('error',''),flush=True)
result['complete']=True
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
