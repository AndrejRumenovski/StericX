"""Alternating precision candidates with exact certificate comparisons."""
import argparse,json,os,statistics,subprocess
from pathlib import Path
from campaign import EVIDENCE,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--output',required=True);p.add_argument('--pairs',type=int,default=8);p.add_argument('--scope',default='Exploratory API timing on CPU 4; legacy independent-reference processes are confined to CPUs 0/1. No whole-host quiet gate.');a=p.parse_args()
out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
inputs=json.loads((EVIDENCE/'molecular_chord_requests.json').read_text());requests=[]
for t in [.05,.01]:
    for q in inputs:requests.append({**q,'id':q['id']+f'/target_{t}','tolerance':t,'budget':1048575})
save(out/'requests.json',requests)
receipt={'scope':a.scope+' One warmup process per candidate followed by eight alternating pairs. Exact certificate/cell/evaluation equality required.','bindings':{str(p):sha(p) for p in [a.baseline,a.candidate,Path(__file__),out/'requests.json']},'cpu_before':Path('/proc/stat').read_text(),'runs':[],'complete':False}
expected={};env=dict(os.environ,RAYON_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for rep in range(-1,a.pairs):
    for label,engine in ([('baseline',a.baseline),('candidate',a.candidate)] if rep%2==0 else [('candidate',a.candidate),('baseline',a.baseline)]):
        folder=out/f'{rep}_{label}';folder.mkdir()
        with (folder/'stdout.jsonl').open('x') as log,(folder/'stderr').open('x') as err:subprocess.run(['taskset','-c','4',str(engine.resolve()),str(out/'requests.json')],env=env,stdout=log,stderr=err,check=True)
        rows=list(map(json.loads,(folder/'stdout.jsonl').read_text().splitlines()))
        for q,r in zip(requests,rows,strict=True):
            assert q['id']==r['id'];science={k:v for k,v in r.items() if k not in ['id','elapsed_ns','allocations']}
            assert science==expected.setdefault(q['id'],science),(q['id'],label,'certificate changed')
        receipt['runs'].append({'rep':rep,'warmup':rep<0,'variant':label,'rows':rows});save(out/'summary.json',receipt)
    print(rep,flush=True)
summary=[]
for q in requests:
    values={label:[next(r['elapsed_ns'] for r in run['rows'] if r['id']==q['id']) for run in receipt['runs'] if run['variant']==label and not run['warmup']] for label in ['baseline','candidate']}
    summary.append({'id':q['id'],'samples_ns':values,'speedup':statistics.median(values['baseline'])/statistics.median(values['candidate'])})
receipt['summary']=summary;receipt['complete']=True;receipt['cpu_after']=Path('/proc/stat').read_text();save(out/'summary.json',receipt)
print(json.dumps({'median_speedup':statistics.median(r['speedup'] for r in summary),'range':[min(r['speedup'] for r in summary),max(r['speedup'] for r in summary)]}))
