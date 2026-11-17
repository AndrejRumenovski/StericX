"""Check copied production chord predicates against exact rational extrema."""
import argparse,json,os,random,shutil,subprocess
from fractions import Fraction as F
from pathlib import Path
from campaign import EVIDENCE,DOCS,sha,save
from validate_chords import sqrt_fraction
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',required=True);p.add_argument('--output',required=True);a=p.parse_args()
source=EVIDENCE/a.build;out=EVIDENCE/a.output;out.mkdir(exist_ok=False);(out/'src').mkdir()
for name in ['Cargo.toml','Cargo.lock']:
    (out/name).write_text((source/'observer'/name).read_text().replace('stericx-ultimate-observer','stericx-predicate-probe'))
prefix=(source/'source/src/geometry/continuous_volume.rs').read_text().replace('//!','//')
(out/'src/observed.rs').write_text(prefix+'\npub fn inspect(q:&serde_json::Value)->serde_json::Value {certified::inspect_predicate(q)}\n')
(out/'src/certified.rs').write_text((source/'source/src/geometry/continuous_volume/certified.rs').read_text()+'\n'+(DOCS/'predicate_probe.rs').read_text())
(out/'src/main.rs').write_text('''pub use steric_x::geometry::Molecule;
pub struct BuriedVolumeCalculator;
#[derive(Debug)] pub struct DescriptorError(pub String);
#[allow(dead_code,unused_imports)] mod observed {include!("observed.rs");}
fn main(){let q:Vec<serde_json::Value>=serde_json::from_str(&std::fs::read_to_string(std::env::args().nth(1).unwrap()).unwrap()).unwrap();for r in q {println!("{}",observed::inspect(&r));}}
''')
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/predicates',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--release','--offline','--locked','--manifest-path',str(out/'Cargo.toml')]
with (out/'build.log').open('x') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
engine=out/'probe';shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/stericx-predicate-probe',engine)
rng=random.Random(725291);requests=[]
for i in range(6000):
    scale=2.**rng.choice([-500,-200,-30,-1,0,5,50,200,300])
    origin=[rng.uniform(-4,4)*scale for _ in range(3)]
    center=[x+rng.uniform(-2,2)*scale for x in origin]
    radius=rng.choice([.5,1.,2.])*scale
    lo=[rng.uniform(-2,2)*scale for _ in range(2)]
    hi=[x+rng.choice([2.**-30,2.**-15,.01,.5])*scale for x in lo]
    if i%4==0:
        lo[0]=center[0]-origin[0]+radius
        hi[0]=lo[0]
        lo[1]=center[1]-origin[1];hi[1]=lo[1]
    requests.append({'id':i,'origin':origin,'sphere':{'center':center,'radius':radius},'lo':lo,'hi':hi})
save(out/'requests.json',requests)
with (out/'observed.jsonl').open('x') as f:subprocess.run([str(engine),str(out/'requests.json')],stdout=f,check=True)
failures=[];regular=empty=0
def contained(bounds,low,high):
    return bounds[0] is not None and bounds[1] is not None and F(bounds[0])<=low<=high<=F(bounds[1])
for q,line in zip(requests,(out/'observed.jsonl').open(),strict=True):
    r=json.loads(line);s=q['sphere'];center=[F(x)-F(y) for x,y in zip(s['center'],q['origin'])]
    lower=upper=F(s['radius'])**2
    for j in range(2):
        l=F(q['lo'][j])-center[j];h=F(q['hi'][j])-center[j]
        minimum=0 if l<=0<=h else min(l*l,h*h)
        lower-=max(l*l,h*h);upper-=minimum
    if not contained(r['q'],lower,upper):failures.append([q['id'],'q'])
    if r['span'] is None:
        empty+=1
        if upper>0:failures.append([q['id'],'false_empty'])
        continue
    # Exact integer sqrt with enough absolute precision for every tested scale.
    lowroot=sqrt_fraction(max(F(0),lower),bits=1200)
    highroot=sqrt_fraction(max(F(0),upper),bits=1200)
    span=r['span'];z=center[2]
    if not contained(span['low'],z-highroot[1],z-lowroot[0]):failures.append([q['id'],'lower_endpoint'])
    if not contained(span['high'],z+lowroot[0],z+highroot[1]):failures.append([q['id'],'upper_endpoint'])
    if span['regular']:
        regular+=1
        if lower<=0:failures.append([q['id'],'false_regular'])
        midpoint=[(F(q['lo'][j])+F(q['hi'][j]))/2 for j in range(2)]
        for j,m in enumerate(midpoint):
            if not contained(r['mid'][j],m,m):failures.append([q['id'],'midpoint'])
        mq=F(s['radius'])**2-sum((midpoint[j]-center[j])**2 for j in range(2))
        sl,sh=sqrt_fraction(mq,bits=1200)
        if not contained(r['mid_endpoints'][0],z-sh,z-sl):failures.append([q['id'],'midpoint_lower'])
        if not contained(r['mid_endpoints'][1],z+sl,z+sh):failures.append([q['id'],'midpoint_upper'])
receipt={'complete':True,'cases':len(requests),'regular':regular,'empty':empty,'failures':failures,
    'reference':'Exact Fraction extrema and integer-sqrt enclosures (1200 fractional bits), no production arithmetic reused.',
    'hashes':{str(p):sha(p) for p in [Path(__file__),DOCS/'predicate_probe.rs',out/'src/observed.rs',out/'src/certified.rs',out/'requests.json',engine]},'build_command':cmd}
save(out/'summary.json',receipt);print(json.dumps({k:v for k,v in receipt.items() if k not in ['hashes','build_command']}))
if failures:raise SystemExit(1)
