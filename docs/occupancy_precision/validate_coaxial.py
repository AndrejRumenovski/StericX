"""Exact piecewise-polynomial references with genuine many-sphere intersections."""
import argparse,itertools,json,random,subprocess,sys
from fractions import Fraction as F
from pathlib import Path
from campaign import ROOT,EVIDENCE,save,sha
sys.path.insert(0,str(ROOT/'scripts'))
from dbstep_resolution_reference import pi_interval,_float_enclosure

def exact_coaxial(region,atoms):
    # Cross-sections normal to the center line are concentric disks. Every
    # squared disk radius is -z²+2cz+(r²-c²). All pair crossings are rational.
    center,radius=map(F,region);spheres=[(center,radius),*[(F(c),F(r)) for c,r in atoms]]
    polynomials=[(2*c,r*r-c*c) for c,r in spheres]
    lower,upper=center-radius,center+radius;cuts={lower,upper}
    for c,r in spheres:cuts.update(x for x in [c-r,c+r] if lower<x<upper)
    for (a,b),(c,d) in itertools.combinations(polynomials,2):
        if a!=c:
            x=(d-b)/(a-c)
            if lower<x<upper:cuts.add(x)
    coefficient=F(0)
    cuts=sorted(cuts)
    for lo,hi in zip(cuts,cuts[1:]):
        mid=(lo+hi)/2
        winner=max(range(1,len(polynomials)),key=lambda i:polynomials[i][0]*mid+polynomials[i][1])
        if polynomials[0][0]*mid+polynomials[0][1]<polynomials[winner][0]*mid+polynomials[winner][1]:winner=0
        a,b=polynomials[winner]
        if -mid*mid+a*mid+b<=0:continue
        coefficient+=-(hi**3-lo**3)/3+a*(hi**2-lo**2)/2+b*(hi-lo)
    pi=pi_interval();return [_float_enclosure(coefficient*pi[0])[0],_float_enclosure(coefficient*pi[1])[1]]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--engine',type=Path,required=True);p.add_argument('--output',required=True);a=p.parse_args()
    out=EVIDENCE/a.output;out.mkdir(exist_ok=False);rng=random.Random(34850);requests=[]
    for i in range(36):
        balls=[(rng.randint(-20,20)/8,rng.randint(4,13)/8) for _ in range(3+i%15)]
        region=(rng.randint(-2,2)/4,2.5+i%3/2)
        reference=exact_coaxial(region,balls)
        for axis in range(3):
            origin=[1024.,-32.,16.] if i%2 else [0.,0.,0.]
            def sphere(c,r):
                center=origin.copy();center[axis]+=c;return {'center':center,'radius':r}
            requests.append({'id':f'coaxial_{i}/axis_{axis}','region':sphere(*region),'atoms':[sphere(*s) for s in balls],'tolerance':.01,'budget':1048575,'reference':reference})
    save(out/'requests.json',requests)
    with (out/'observed.jsonl').open('x') as log,(out/'stderr').open('x') as err:subprocess.run([str(a.engine.resolve()),str(out/'requests.json')],stdout=log,stderr=err,check=True)
    results=list(map(json.loads,(out/'observed.jsonl').read_text().splitlines()));failures=[]
    for q,r in zip(requests,results,strict=True):
        lo,hi=q['reference']
        if not r['lower']<=lo<=hi<=r['upper']:failures.append([q['id'],'containment',q['reference'],r])
        if r['target_met']!=(r['absolute_error_bound']<=q['tolerance']):failures.append([q['id'],'target_flag'])
    receipt={'complete':True,'cases':len(results),'failures':failures,'unmet':[r['id'] for r in results if not r['target_met']],
        'reference':'Exact rational piecewise quadratic disk-envelope integral times rational Machin-pi enclosure; supports arbitrary triple and higher intersections. Three exact coordinate-axis permutations; translated binary-exact centers.',
        'input_sha256':sha(out/'requests.json'),'engine_sha256':sha(a.engine),'controller_sha256':sha(Path(__file__))}
    save(out/'summary.json',receipt);print(json.dumps(receipt))
    if failures:raise SystemExit(1)
if __name__=='__main__':main()
