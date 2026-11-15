#!/usr/bin/env python3
"""Construct exact f32 inputs, validate f64 Sterimol against independent search."""
import argparse
import hashlib
import json
import random
import shutil
import statistics
import struct
import subprocess
from pathlib import Path
import mpmath as mp
from continuous_support_reference import minimum_support

ROOT = Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--engine',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve();out.mkdir(exist_ok=False)
for source in [Path(__file__),Path(__file__).with_name('continuous_support_reference.py')]:
    shutil.copyfile(source,out/source.name)
f32=lambda x:struct.unpack('f',struct.pack('f',x))[0]
cases=[]
def add(name,discs):
    cases.append({'id':name,'base':0,'dummy':[0.,0.,-1.],
                  'atoms':[{'position':[0.,0.,0.],'radius':0.}]+[
                      {'position':[f32(x),f32(y),0.],'radius':f32(r)} for x,y,r in discs]})
add('single',[(2,0,3)])
add('single_clipped',[(2,0,.5)])
add('opposed',[(2,0,.5),(-2,0,.5)])
add('square',[(2,2,.5),(-2,2,.5),(-2,-2,.5),(2,-2,.5)])
add('nested',[(0,0,3),(1,0,1)])
add('duplicates',[(1,.25,2),(1,.25,2),(1,.25,1)])
add('near_tangent_support',[(0,1,1),(1,1,f32(2-2**-23)),(-1,-1,.1)])
add('tiny',[(1e-30,2e-30,3e-30),(-2e-30,1e-30,2e-30)])
add('large',[(1e30,2e30,3e30),(-2e30,1e30,2e30)])
rng=random.Random(20261115)
for i in range(18):
    n=[3,8,32][i%3]
    add(f'random_{i}_{n}',[(rng.uniform(-4,4),rng.uniform(-4,4),rng.uniform(0,2.5)) for _ in range(n)])
source_observations=ROOT/'.stericx/ultimate_engine/baseline_observations/geometry/stdout.jsonl'
source_requests=ROOT/'.stericx/profiling/scientifically_validated_optimization/scientific_oracle_v3/requests/geometry.jsonl'
requests={d['id']:d for d in map(json.loads,source_requests.read_text().splitlines())}
for row in list(map(json.loads,source_observations.read_text().splitlines()))[:12]:
    cases.append({'id':'real_'+row['id'],'base':requests[row['id']]['donor'],'dummy':row['center'],
                  'atoms':[{'position':atom['position'],'radius':atom['vdw_radius']} for atom in row['atoms']]})
(out/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
command=['taskset','-c','2',str(a.engine.resolve()),str(out/'cases.json')]
with (out/'native.jsonl').open('w') as stdout, (out/'native.stderr').open('w') as stderr:
    subprocess.run(command,stdout=stdout,stderr=stderr,check=True)
native=[json.loads(line) for line in (out/'native.jsonl').read_text().splitlines()]
receipt={'command':command,'engine_sha256':hashlib.sha256(a.engine.read_bytes()).hexdigest(),
         'input_sha256':hashlib.sha256((out/'cases.json').read_bytes()).hexdigest(),
         'sources':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('continuous_support_reference.py'),source_observations,source_requests]},
         'scope':'80-decimal angular interval search, no production event enumeration',
         'rows':[]}
for case,observed in zip(cases,native,strict=True):
    assert case['id']==observed['id']
    with mp.workdps(80):
        origin=mp.matrix(case['dummy']);axis=mp.matrix(case['atoms'][case['base']]['position'])-origin
        z=axis/mp.norm(axis)
        seed=mp.matrix([1 if i==min(range(3),key=lambda j:abs(z[j])) else 0 for i in range(3)])
        x=seed-mp.fdot(seed,z)*z;x/=mp.norm(x)
        y=mp.matrix([z[1]*x[2]-z[2]*x[1],z[2]*x[0]-z[0]*x[2],z[0]*x[1]-z[1]*x[0]])
        discs=[];lengths=[];radials=[]
        for atom in case['atoms']:
            relative=mp.matrix(atom['position'])-origin;radius=mp.mpf(atom['radius'])
            axial=mp.fdot(relative,z);px=mp.fdot(relative,x);py=mp.fdot(relative,y)
            discs.append((px,py,radius));lengths.append(axial+radius)
            radials.append(mp.sqrt(px*px+py*py)+radius)
        expected_l=max(lengths);expected_b5=max(radials)
    ref=minimum_support(discs,relative_width=1e-12)
    assert ref['converged'],(case['id'],ref)
    scale=float(max(abs(v) for d in discs for v in d))
    arithmetic_margin=128*2**-52*scale
    row={'id':case['id'],'reference':ref,'arithmetic_margin':arithmetic_margin,'methods':[]}
    for measurement in observed['measurements']:
        value=measurement['result']['b1']
        exact=measurement['method']=='continuous_events_f64'
        if exact:
            assert ref['lower']-arithmetic_margin <= value <= ref['upper']+arithmetic_margin, (case['id'],value,ref)
        low_error=max(0.,ref['lower']-value,value-ref['upper'])
        high_error=max(abs(value-ref['lower']),abs(value-ref['upper']))
        row['methods'].append({**measurement,'median_ns':statistics.median(measurement['elapsed_ns'])/measurement['iterations'],
                               'absolute_error_lower':low_error,'absolute_error_upper':high_error})
    # Independent 80-decimal projection and axial/radial formulas.
    with mp.workdps(80):
        continuous=observed['measurements'][0]['result']
        assert abs(continuous['l']-float(expected_l)) <= 32*2**-52*max(abs(float(expected_l)),1e-300)
        assert abs(continuous['b5']-float(expected_b5)) <= 32*2**-52*max(abs(float(expected_b5)),1e-300)
    receipt['rows'].append(row)
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(case['id'],'PASS',ref['intervals'],flush=True)
receipt['complete']=True
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
