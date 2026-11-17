"""Gated common-grid radius-profile timing with fresh/warm scopes separated."""
import argparse,gzip,json,os,subprocess,sys
from pathlib import Path
from campaign import ROOT,EVIDENCE,sha,save
sys.path.insert(0,str(ROOT/'scripts'));import profile_stericx as h
sys.path.insert(0,str(ROOT/'docs/ultimate_engine'))
from host_activity import measure,postflight
from benchmark_artifacts import compact_metrics
p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['preflight','timing']);p.add_argument('--output',required=True);a=p.parse_args()
out=EVIDENCE/a.output;out.mkdir(exist_ok=False);gate=EVIDENCE/'competitive_radius_gate/summary.json';g=json.loads(gate.read_text());assert g['complete'] and g['passed']
native=EVIDENCE/'competitive_radius_build/radius_workflow';worker=Path(__file__).with_name('radius_worker.py');python=ROOT/'.stericx/paton/reference/venv/bin/python'
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1');env.pop('STERICX_PROFILE_PATH',None)
receipt={'phase':a.phase,'scope':g['scope'],'bindings':{str(p):sha(p) for p in [gate,native,worker,Path(__file__)]},'runs':[],'complete':False}
affinities={1:[2],2:[2,3],4:[0,1,2,3],6:list(range(6)),12:list(range(12))}
if a.phase=='timing':
    receipt['host_preflight']=measure();save(out/'summary.json',receipt)
    if not receipt['host_preflight']['quiet']:raise SystemExit('Quiet window required before radius-profile timing.')
launcher=h.build_launcher(out/'launcher') if a.phase=='timing' else None
configs=[(12,t) for t in affinities] if a.phase=='preflight' else [(1,1),(12,1),(12,6)]
for count,threads in configs:
    source=g['rows'][:count];paths=[r['path'] for r in source]
    for scope in (['preflight'] if a.phase=='preflight' else ['fresh_process','warm_workflow']):
        for rep in ([0] if a.phase=='preflight' else range(-1,8)):
            for tool in (['native','morfeus','dbstep'] if rep%2==0 else ['dbstep','morfeus','native']):
                iterations=1 if a.phase=='preflight' else (20 if count==1 else 3) if scope=='warm_workflow' else 0
                cmd=[str(native),str(iterations),*paths] if tool=='native' else [str(python),str(worker),'--tool',tool,'--workers',str(threads),'--warm-iterations',str(iterations),*paths]
                folder=out/f'n{count}_t{threads}_{scope}_{rep}_{tool}';child_env=dict(env,RAYON_NUM_THREADS=str(threads))
                if launcher:
                    metrics=compact_metrics(h.measure(cmd,child_env,set(affinities[threads]),folder,launcher),folder);assert metrics.get('exit_code',metrics.get('returncode',0))==0;raw=(folder/'stdout.txt').read_bytes()
                else:
                    folder.mkdir();proc=subprocess.run(['taskset','-c',','.join(map(str,affinities[threads])),*cmd],env=child_env,capture_output=True,check=True);raw=proc.stdout;(folder/'stdout.txt').write_bytes(raw);(folder/'stderr.txt').write_bytes(proc.stderr);metrics=None
                value=json.loads(raw);assert value['rows']==[r[tool] for r in source],(tool,count,threads)
                receipt['runs'].append({'count':count,'threads':threads,'scope':scope,'rep':rep,'warmup':rep<0,'tool':tool,'metrics':metrics,'warm_elapsed_ns':value['warm_elapsed_ns'],'exact':True,'command':cmd,'stdout_sha256':sha(folder/'stdout.txt')});save(out/'summary.json',receipt)
        print(count,threads,scope,flush=True)
for path,digest in receipt['bindings'].items():assert sha(Path(path))==digest
if a.phase=='timing':receipt['host_postflight']=postflight(out)
receipt['complete']=True;save(out/'summary.json',receipt)
