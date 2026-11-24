"""Independent rational cross-section envelope and cap references for regions."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,itertools,json,random,subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from dbstep_resolution_reference import pi_interval,_float_enclosure

def coefficient(region,atoms,positive):
    c,r=map(F,region);spheres=[(c,r),*[(F(c),F(r)) for c,r in atoms]]
    lower,upper=c-r,c+r
    if positive:lower=max(lower,c)
    else:upper=min(upper,c)
    if not atoms:return F(0)
    polynomials=[(2*c,r*r-c*c) for c,r in spheres]
    cuts={lower,upper}
    for c,r in spheres:cuts.update(x for x in (c-r,c+r) if lower<x<upper)
    for (a,b),(c,d) in itertools.combinations(polynomials,2):
        if a!=c:
            x=(d-b)/(a-c)
            if lower<x<upper:cuts.add(x)
    total=F(0);cuts=sorted(cuts)
    for lo,hi in zip(cuts,cuts[1:]):
        mid=(lo+hi)/2
        i=max(range(1,len(polynomials)),key=lambda i:polynomials[i][0]*mid+polynomials[i][1])
        if polynomials[0][0]*mid+polynomials[0][1]<polynomials[i][0]*mid+polynomials[i][1]:i=0
        a,b=polynomials[i]
        if -mid*mid+a*mid+b>0:total+=-(hi**3-lo**3)/3+a*(hi**2-lo**2)/2+b*(hi-lo)
    return total

def interval(c):
    lo,hi=pi_interval()
    return [_float_enclosure(c*lo)[0],_float_enclosure(c*hi)[1]]

def make_case(name,axial_region,balls,axis,origin,tol=.02,budget=262144):
    def sphere(c,r):
        center=origin.copy();center[axis]+=c
        return dict(center=center,radius=r)
    low=coefficient(axial_region,balls,False);high=coefficient(axial_region,balls,True)
    signs=[(1,1,1),(-1,1,1),(-1,-1,1),(1,-1,1),(1,1,-1),(-1,1,-1),(-1,-1,-1),(1,-1,-1)]
    octants=[(high if side[axis]>0 else low)/4 for side in signs]
    quadrants=[octants[i]+octants[i+4] for i in range(4)]
    denom=F(4,3)*F(axial_region[1])**3
    return dict(id=name,region=sphere(*axial_region),atoms=[sphere(*s) for s in balls],
        tolerance=tol,budget=budget,references=dict(total=interval(high+low),quadrants=list(map(interval,quadrants)),octants=list(map(interval,octants))),
        percent_references=dict(total=float(100*(high+low)/denom),quadrants=[float(400*x/denom) for x in quadrants],octants=[float(800*x/denom) for x in octants]))

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--engine',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    rng=random.Random(20261124);qs=[]
    for i in range(24):
        balls=[(rng.randint(-16,16)/8,rng.randint(2,12)/8) for _ in range(1+i%9)]
        for axis in range(3):
            qs.append(make_case(f'coaxial_{i}_axis{axis}',(0.,2.5),balls,axis,[1024.,-32.,16.] if i%2 else [0.,0.,0.]))
    for i,d in enumerate([0.,0.5,1.999,2.,2.0001,-0.5]):
        for axis in range(3):qs.append(make_case(f'cap_contact_{i}_axis{axis}',(0.,1.),[(d,1.)],axis,[0.,0.,0.]))
    for b in [4,32,1024,16384]:qs.append(make_case(f'budget_{b}',(0.,2.),[(-.5,1.),(.25,1.),(.5,1.)],2,[0.,0.,0.],1e-8,b))
    qs.append(make_case('empty',(0.,2.),[],2,[0.,0.,0.]))
    # Eight disjoint balls strictly inside eight octants: independent exact volumes.
    signs=[(1,1,1),(-1,1,1),(-1,-1,1),(1,-1,1),(1,1,-1),(-1,1,-1),(-1,-1,-1),(1,-1,-1)]
    spheres=[dict(center=list(s),radius=(i+1)/32) for i,s in enumerate(signs)]
    os=[F(4,3)*F(s['radius'])**3 for s in spheres];qv=[os[i]+os[i+4] for i in range(4)];v=sum(os);den=F(36)
    qs.append(dict(id='eight_disjoint_octants',region=dict(center=[0.,0.,0.],radius=3.),atoms=spheres,tolerance=.01,budget=262144,references=dict(total=interval(v),quadrants=list(map(interval,qv)),octants=list(map(interval,os))),percent_references=dict(total=float(100*v/den),quadrants=[float(400*x/den) for x in qv],octants=[float(800*x/den) for x in os])))
    path=a.output/'requests.json';path.write_text(json.dumps(qs,indent=2)+'\n')
    with (a.output/'observed.jsonl').open('x') as out,(a.output/'process_metrics.txt').open('x') as err:
        subprocess.run(['/usr/bin/time','-v',str(a.engine.resolve()),str(path)],stdout=out,stderr=err,check=True)
    rs=list(map(json.loads,(a.output/'observed.jsonl').read_text().splitlines()));failures=[];unmet=[]
    for q,row in zip(qs,rs,strict=True):
        assert q['id']==row['id'];r=row.get('result')
        if r is None:failures.append([q['id'],'error',row]);continue
        if not r['target_met']:unmet.append(q['id'])
        for kind in ['total','quadrants','octants']:
            ms=[r[kind]] if kind=='total' else r[kind]
            refs=[q['references'][kind]] if kind=='total' else q['references'][kind]
            prs=[q['percent_references'][kind]] if kind=='total' else q['percent_references'][kind]
            for i,(m,ref,pr) in enumerate(zip(ms,refs,prs,strict=True)):
                if not m['lower_bound_angstrom3']<=ref[0]<=ref[1]<=m['upper_bound_angstrom3']:failures.append([q['id'],kind,i,'containment',m,ref])
                if m['target_met']!=(m['absolute_error_bound_angstrom3']<=q['tolerance']):failures.append([q['id'],kind,i,'target_flag'])
                pct=m['percent_buried_volume']
                if not pct['lower_bound']<=pr<=pct['upper_bound']:failures.append([q['id'],kind,i,'percent',pct,pr])
        if r['cells_evaluated']>q['budget'] or r['cells_evaluated']!=4+2*r['subdivisions']:failures.append([q['id'],'work_count'])
        for i,m in enumerate(r['quadrants']):
            bins=[r['octants'][i],r['octants'][i+4]]
            l=sum(F(x['lower_bound_angstrom3']) for x in bins);h=sum(F(x['upper_bound_angstrom3']) for x in bins)
            if not F(m['lower_bound_angstrom3'])<=l<=h<=F(m['upper_bound_angstrom3']):failures.append([q['id'],'regional_sum',i])
    receipt=dict(cases=len(qs),regional_integrals=13*len(qs),failures=failures,unmet=unmet,input_sha256=sha(path),observed_sha256=sha(a.output/'observed.jsonl'),engine_sha256=sha(a.engine),controller_sha256=sha(Path(__file__)),reference='Exact rational quadratic envelope integration split at the regional plane, times rational Machin-pi enclosure; sphere cap/union/disjoint analytical references. Percent references cancel pi algebraically.')
    (a.output/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
    if failures:raise SystemExit(1)
if __name__=='__main__':main()
