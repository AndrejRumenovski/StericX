#!/usr/bin/env python3
"""Unchanged external APIs with explicitly disclosed unrounded DBSTEP alignment."""
import argparse, json, time
from pathlib import Path
import numpy as np
from morfeus import Sterimol
from dbstep import sterics

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('input',type=Path);p.add_argument('tool',choices=['morfeus','dbstep'])
p.add_argument('iterations',type=int);a=p.parse_args();assert a.iterations>0
cases=[]
for case in json.loads(a.input.read_text()):
    coords=np.array([case['dummy']]+[r['position'] for r in case['atoms']],dtype=float)
    radii=np.array([0.]+[r['radius'] for r in case['atoms']],dtype=float)
    cases.append((case,coords,radii,['Ni']+['C']*len(case['atoms']),['Bq']+['C']*len(case['atoms'])))

def batch():
    rows=[]
    for case,coords,radii,elements,db_elements in cases:
        if a.tool=='morfeus':
            m=Sterimol(elements,coords.copy(),dummy_index=1,attached_index=case['base']+2,radii=radii,n_rot_vectors=361)
            l,b1,b5=m.L_value_uncorrected,m.B_1_value,m.B_5_value
        else:
            relative=coords-coords[0]
            axis=relative[case['base']+1]/np.linalg.norm(relative[case['base']+1])
            seed=np.eye(3)[int(np.argmin(np.abs(axis)))]
            x=seed-np.dot(seed,axis)*axis;x/=np.linalg.norm(x);y=np.cross(axis,x)
            common=np.column_stack([relative@x,relative@y,relative@axis])
            l,b5,b1,_=sterics.get_classic_sterimol(common,radii,db_elements)
        rows.append({'id':case['id'],'l':float(l),'b1':float(b1),'b5':float(b5)})
    return rows

expected=batch();started=time.perf_counter_ns()
for _ in range(a.iterations):rows=batch()
elapsed=time.perf_counter_ns()-started
assert rows==expected
print(json.dumps({'rows':rows,'elapsed_ns':elapsed,'iterations':a.iterations}))
