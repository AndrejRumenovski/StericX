"""Independent analytical unions plus frozen molecular convergence checks."""
from fractions import Fraction as F
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
from campaign import ROOT,EVIDENCE,save,sha
sys.path.insert(0,str(ROOT/'scripts'))
from dbstep_resolution_reference import pi_interval,sphere_intersection_interval,_float_enclosure


def sqrt_fraction(q, bits=160):
    """Exact integer square-root enclosure; no floating-point square root."""
    k=math.isqrt((q.numerator << (2*bits))//q.denominator)
    scale=1<<bits
    low=F(k,scale)
    return low,low if low*low==q else F(k+1,scale)


def distance2(a,b):
    return sum((F(x)-F(y))**2 for x,y in zip(a,b,strict=True))


def triangle_free_union(region,atoms):
    """Inclusion-exclusion ONLY after proving every triple has a disjoint pair."""
    pil,pih=pi_interval()
    coefficient=sum(F(4,3)*F(a['radius'])**3 for a in atoms)
    lower,upper=coefficient*pil,coefficient*pih
    edges=set()
    for i,a in enumerate(atoms):
        assert distance2(a['center'],region['center']) <= (F(region['radius'])-F(a['radius']))**2
        for j,b in enumerate(atoms[:i]):
            q=distance2(a['center'],b['center'])
            if q >= (F(a['radius'])+F(b['radius']))**2:continue
            edges.add((j,i))
            dl,dh=sqrt_fraction(q)
            lens_low=sphere_intersection_interval(a['radius'],b['radius'],dh)[0]
            lens_high=sphere_intersection_interval(a['radius'],b['radius'],dl)[1]
            lower-=F(lens_high);upper-=F(lens_low)
    for i,j,k in itertools.combinations(range(len(atoms)),3):
        assert not all(e in edges for e in [(i,j),(i,k),(j,k)]),'Triple overlaps invalidate pair-only formula'
    return [_float_enclosure(lower)[0],_float_enclosure(upper)[1]]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,required=True);p.add_argument('--output',required=True)
    a=p.parse_args();out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
    rng=random.Random(20261117)
    requests=[]
    for case in range(80):
        count=3+case%10
        v=[rng.uniform(-1,1) for _ in range(3)];norm=math.sqrt(sum(x*x for x in v));v=[x/norm for x in v]
        origin=[rng.choice([0.,16.,-32.,1024.]) for _ in range(3)]
        atoms=[{'center':[origin[j]+v[j]*(i-(count-1)/2)*.9 for j in range(3)],
                'radius':rng.randint(7,9)/16} for i in range(count)]
        region={'center':origin,'radius':count*.5+1.0}
        reference=triangle_free_union(region,atoms)
        requests.append({'id':f'triangle_free_{case}','region':region,'atoms':atoms,
                         'reference':reference,'tolerance':.01,'budget':262143})
    # Exact coordinate permutations/sign changes preserve the same represented geometry.
    for case in requests[:8].copy():
        def transform(p):return [p[2],-p[0],p[1]]
        requests.append({**case,'id':case['id']+'_rotated','region':{**case['region'],'center':transform(case['region']['center'])},
                         'atoms':[{**atom,'center':transform(atom['center'])} for atom in case['atoms']]})
        requests.append({**case,'id':case['id']+'_permuted','atoms':list(reversed(case['atoms']))})
    molecules=json.loads((EVIDENCE/'molecular_chord_requests.json').read_text())
    for case in molecules:
        for budget in [1023,8191,65535,262143]:
            requests.append({**case,'id':case['id']+f'/budget_{budget}','budget':budget,'tolerance':1e-12})
        for tolerance in [.1,.05,.01]:
            requests.append({**case,'id':case['id']+f'/target_{tolerance}','tolerance':tolerance,'budget':1048575})
    save(out/'requests.json',requests)
    with (out/'observed.jsonl').open('x') as log,(out/'stderr').open('x') as err:
        subprocess.run(['taskset','-c','2',str(a.engine.resolve()),str(out/'requests.json')],stdout=log,stderr=err,check=True)
    observations=list(map(json.loads,(out/'observed.jsonl').read_text().splitlines()))
    failures=[]
    for q,r in zip(requests,observations,strict=True):
        assert q['id']==r['id']
        if 'reference' in q:
            lo,hi=q['reference']
            if not r['lower']<=lo<=hi<=r['upper']:failures.append([q['id'],'analytic_containment',q['reference'],r])
        if not r['lower']<=r['volume']<=r['upper']:failures.append([q['id'],'estimate_outside'])
        if r['target_met']!=(r['absolute_error_bound']<=q['tolerance']):failures.append([q['id'],'target_flag'])
        if r['cells']>q['budget']:failures.append([q['id'],'budget'])
    by_id={r['id']:r for r in observations}
    for case in molecules:
        last=None
        for budget in [1023,8191,65535,262143]:
            r=by_id[case['id']+f'/budget_{budget}']
            if last and not (last['lower']<=r['lower']<=r['upper']<=last['upper']):failures.append([r['id'],'not_nested'])
            last=r
    for case in requests[:8]:
        def science(r):return {k:v for k,v in r.items() if k not in ['id','elapsed_ns']}
        if science(by_id[case['id']])!=science(by_id[case['id']+'_permuted']):failures.append([case['id'],'permutation'])
    receipt={'complete':True,'rows':len(observations),'failures':failures,
        'input_sha256':sha(out/'requests.json'),'engine_sha256':sha(a.engine),'controller_sha256':sha(Path(__file__)),
        'analytical_reference':'Exact rational coordinates/radii and integer-sqrt distance enclosures, rational Machin pi, pair inclusion-exclusion only for a verified triangle-free overlap graph.',
        'targets':{str(t):[{'id':c['id'],**by_id[c['id']+f'/target_{t}']} for c in molecules] for t in [.1,.05,.01]},
        'unmet':[r['id'] for r in observations if not r['target_met']]}
    save(out/'summary.json',receipt)
    print(json.dumps({'rows':len(observations),'failures':failures,'molecular_targets':{k:sum(r['target_met'] for r in v) for k,v in receipt['targets'].items()}},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
