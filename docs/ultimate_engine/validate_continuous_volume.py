#!/usr/bin/env python3
"""Independent rational/analytic validation and accuracy-cost sampling of volume."""
import argparse
import hashlib
import json
import math
import random
import shutil
import statistics
import struct
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from dbstep_resolution_reference import pi_interval,sphere_intersection_interval

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--engine',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve();out.mkdir(exist_ok=False)
shutil.copyfile(__file__,out/'controller.py')
shutil.copyfile(ROOT/'scripts/dbstep_resolution_reference.py',out/'independent_reference.py')
f32=lambda x:struct.unpack('f',struct.pack('f',x))[0]
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
pil,pih=pi_interval()
assert Fraction.from_float(math.nextafter(math.pi,-math.inf))<pil<pih<Fraction.from_float(math.nextafter(math.pi,math.inf))
(out/'pi_certificate.json').write_text(json.dumps({'machin_lower':str(pil),'machin_upper':str(pih),
    'float_lower':math.nextafter(math.pi,-math.inf).hex(),'float_upper':math.nextafter(math.pi,math.inf).hex()},indent=2)+'\n')
def enclose_fraction(value):
    v=float(value);q=Fraction.from_float(v)
    return (math.nextafter(v,-math.inf) if q>value else v,
            math.nextafter(v,math.inf) if q<value else v)
def coefficient_interval(c):
    return enclose_fraction(c*pil)[0],enclose_fraction(c*pih)[1]
def lens_coefficient(radius,atom_radius,separation):
    r,s,d=map(Fraction.from_float,map(float,(radius,atom_radius,separation)))
    if d>=r+s:return Fraction(0)
    if d<=abs(r-s):return Fraction(4,3)*min(r,s)**3
    return (r+s-d)**2*(d*d+2*d*(r+s)-3*(r-s)**2)/(12*d)
cases=[];references={};requests=[]
def add(name,r,centers,radii,coefficient=None):
    atoms=[{'center':c,'radius':v} for c,v in zip(centers,radii,strict=True)]
    cases.append({'id':name,'region':{'center':[0.,0.,0.],'radius':r},'atoms':atoms})
    if coefficient is not None: bounds=coefficient_interval(coefficient)
    elif not atoms: bounds=(0.,0.)
    elif len(atoms)==1:
        nonzero=[abs(x) for x in centers[0] if x!=0.]
        assert len(nonzero)<=1,'single-ball reference needs explicit exact scalar distance'
        bounds=sphere_intersection_interval(r,radii[0],nonzero[0] if nonzero else 0.)
    elif len(atoms)==2:
        # The original two-ball fixture is wholly inside the integration ball.
        assert all(math.dist(c,[0,0,0])+v<r for c,v in zip(centers,radii,strict=True))
        coefficient=(lens_coefficient(r,radii[0],0.)+lens_coefficient(r,radii[1],0.)
                     -lens_coefficient(radii[0],radii[1],abs(centers[0][0]-centers[1][0])))
        bounds=coefficient_interval(coefficient)
    else: raise AssertionError(name)
    references[name]=bounds
for c in json.loads((ROOT/'docs/ultimate_engine/analytical_cases.json').read_text()):
    add(c['id'],f32(c['radius']),[[f32(v) for v in x] for x in c['centers']],[f32(r) for r in c['radii']])
add('three_identical',2.,[[0.,0.,0.]]*3,[1.]*3,Fraction(4,3))
add('three_overlapping_chain',3.5,[[-1.,0.,0.],[0.,0.,0.],[1.,0.,0.]],[1.]*3,Fraction(19,6))
add('three_disjoint',4.,[[-2.,0.,0.],[0.,0.,0.],[2.,0.,0.]],[.5]*3,Fraction(1,2))
for case in cases:
    for n in [16,32,64,128,256]:
        requests.append({**case,'id':case['id']+f'/midpoint_{n}',
                         'method':{'method':'midpoint_slices','subdivisions':n},'samples':5})
    budgets=[1023,8191,65535] if len(case['atoms'])>=3 else [65535]
    for budget in budgets:
        requests.append({**case,'id':case['id']+f'/reference_{budget}',
                         'method':{'method':'reference','absolute_tolerance':1e-5,'max_boxes':budget},'samples':3})
rng=random.Random(20261115)
for i in range(100):
    r=rng.randint(1,64)/16;s=rng.randint(1,48)/16;d=rng.randint(0,128)/16
    add(f'random_lens_{i}',r,[[d,0.,0.]],[s])
for scale in [1e-200,1e-150,1e-100,1e-30,1.,1e30,1e99]:
    for label,d in [('overlap',scale),('external_below',math.nextafter(2*scale,0.)),
                    ('external_equal',2*scale),('external_above',math.nextafter(2*scale,math.inf))]:
        add(f'boundary_{scale}_{label}',scale,[[d,0.,0.]],[scale])
for case in cases[11:]:
    requests.append({**case,'id':case['id']+'/reference_stress',
        'method':{'method':'reference','absolute_tolerance':max(math.ulp(0.),case['region']['radius']**3*1e-10),'max_boxes':1023},'samples':1})
for original in cases[:11]:
    # Exact signed coordinate permutation and dyadic translation. This tests
    # mathematical rigid invariance; finite XY quadrature need not be identical.
    transform=lambda p:[p[2]+1000.,p[0]-500.,-p[1]+250.]
    variant={**original,'id':original['id']+'_rigid','region':{**original['region'],'center':transform([0.,0.,0.])},
             'atoms':[{'center':transform(atom['center']),'radius':atom['radius']} for atom in original['atoms']]}
    for old,new in zip(original['atoms'],variant['atoms'],strict=True):
        actual=[Fraction(new['center'][i])-Fraction(variant['region']['center'][i]) for i in range(3)]
        assert actual==[Fraction(old['center'][2]),Fraction(old['center'][0]),-Fraction(old['center'][1])]
    references[variant['id']]=references[original['id']]
    for method in [{'method':'midpoint_slices','subdivisions':128},
                   {'method':'reference','absolute_tolerance':1e-5,'max_boxes':65535}]:
        requests.append({**variant,'id':variant['id']+'/'+method['method'],'method':method,'samples':3})
        permutation={**original,'id':original['id']+'_permuted','atoms':list(reversed(original['atoms']))}
        references[permutation['id']]=references[original['id']]
        requests.append({**permutation,'id':permutation['id']+'/'+method['method'],'method':method,'samples':3})
(out/'requests.json').write_text(json.dumps(requests,indent=2)+'\n')
command=['taskset','-c','2',str(a.engine.resolve()),str(out/'requests.json')]
with (out/'native.jsonl').open('w') as stdout,(out/'native.stderr').open('w') as stderr:
    subprocess.run(command,stdout=stdout,stderr=stderr,check=True)
native=[json.loads(line) for line in (out/'native.jsonl').read_text().splitlines()]
receipt={'command':command,'engine_sha256':sha(a.engine),'input_sha256':sha(out/'requests.json'),
         'reference_sha256':sha(out/'independent_reference.py'),'controller_sha256':sha(out/'controller.py'),
         'pi_neighbor_bracket_verified':True,'rows':[],'failures':[]}
for request,observed in zip(requests,native,strict=True):
    assert observed['id']==request['id']
    name=request['id'].split('/')[0];lo,hi=references[name];result=observed['result']
    value=result['volume_angstrom3']
    row={**observed,'reference_interval':[lo,hi],'median_ns':statistics.median(observed['elapsed_ns']),
         'absolute_error_lower':max(0.,lo-value,value-hi),'absolute_error_upper':max(abs(lo-value),abs(hi-value))}
    if result['lower_bound_angstrom3'] is not None:
        contained=result['lower_bound_angstrom3']<=lo and result['upper_bound_angstrom3']>=hi
        error=result['absolute_error_bound_angstrom3']
        consistent=(result['target_met']==(error<=request['method']['absolute_tolerance']))
        row['independent_interval_contained']=contained;row['target_flag_consistent']=consistent
        row['work_budget_obeyed']=result['work_units']<=request['method']['max_boxes']
        if not (contained and consistent and row['work_budget_obeyed']):receipt['failures'].append(request['id'])
    receipt['rows'].append(row)
receipt['complete']=True
by_id={row['id']:row for row in receipt['rows']}
receipt['permutation_checks']=[]
for case in cases[:11]:
    for suffix,original_suffix in [('midpoint_slices','midpoint_128'),('reference','reference_65535')]:
        left=by_id[case['id']+'/'+original_suffix]['result'];right=by_id[case['id']+'_permuted/'+suffix]['result']
        same=left==right
        receipt['permutation_checks'].append({'id':case['id'],'method':suffix,'exact':same})
        if not same:receipt['failures'].append(case['id']+'_permutation/'+suffix)
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'rows':len(receipt['rows']),'failures':receipt['failures'],
                  'unmet_targets':[r['id'] for r in receipt['rows'] if r['result']['target_met'] is False]},indent=2))
