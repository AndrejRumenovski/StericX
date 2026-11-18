"""Fixed 0.05 Å³ batch-scaling supplement, preserving every individual certificate."""
import argparse,json,os,subprocess,sys
from pathlib import Path
from campaign import ROOT,EVIDENCE,DOCS,sha,save
sys.path.insert(0,str(ROOT/'scripts'));import profile_stericx as harness
sys.path.insert(0,str(ROOT/'docs/ultimate_engine'))
from host_activity import measure,postflight
from benchmark_artifacts import compact_metrics
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);p.add_argument('--preflight',action='store_true');p.add_argument('--start-at',default='certified_scaling');a=p.parse_args()
out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
engine=EVIDENCE/'certified_batch_build/certified_batch'
requests=[{**q,'tolerance':.05,'budget':1048575} for q in json.loads((EVIDENCE/'molecular_chord_requests.json').read_text())];assert len(requests)==12
save(out/'requests.json',requests)
reference=EVIDENCE/'accuracy_preflight_v2/summary.json'
expected={r['source_id']:{k:v for k,v in r['result'].items() if k!='elapsed_seconds'} for r in json.loads(reference.read_text())['runs'] if r['configuration']=='certified_0.05'}
expected_rows=[{'id':q['id'],'result':expected[q['id']]} for q in requests]
receipt={'complete':False,'preflight_only':a.preflight,'scope':'Twelve independent original molecules at certified 0.05 Å³; external Rayon batch pool, no internal parallelism or shared geometry cache. Pool construction/input parsing/one warmup are outside internal batch timing; process metrics include them. Eight alternating worker-count orders plus retained warmup round. Ordered results must equal each frozen certificate exactly.',
 'bindings':{str(p):sha(p) for p in [Path(__file__),engine,reference,out/'requests.json',EVIDENCE/'certified_batch_build/receipt.json']},'runs':[]}
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');env.pop('STERICX_PROFILE_PATH',None)
affinity={1:[2],2:[2,3],4:list(range(4)),6:list(range(6)),12:list(range(12))}
if not a.preflight:
    prior=json.loads((EVIDENCE/'certified_scaling_preflight/summary.json').read_text());assert prior['complete'] and len(prior['runs'])==5
    receipt['host_preflight']=measure();save(out/'summary.json',receipt)
    if not receipt['host_preflight']['quiet']:raise SystemExit('Quiet gate failed before certified batch scaling')
launcher=harness.build_launcher(out/'launcher') if not a.preflight else None
for rep in ([0] if a.preflight else range(-1,8)):
    for threads in (list(affinity) if rep%2==0 else list(reversed(affinity))):
        folder=out/f't{threads}_rep{rep}';cmd=[str(engine),str(out/'requests.json'),str(threads)]
        if launcher:
            metrics=compact_metrics(harness.measure(cmd,env,set(affinity[threads]),folder,launcher),folder);assert metrics.get('returncode',metrics.get('exit_code',0))==0
            raw=(folder/'stdout.txt').read_bytes()
        else:
            folder.mkdir();r=subprocess.run(['taskset','-c',','.join(map(str,affinity[threads])),*cmd],env=env,capture_output=True,check=True);raw=r.stdout;(folder/'stdout.txt').write_bytes(raw);(folder/'stderr.txt').write_bytes(r.stderr);metrics=None
        value=json.loads(raw);assert value['rows']==expected_rows and value['threads']==threads
        receipt['runs'].append({'threads':threads,'rep':rep,'warmup':rep<0,'exact':True,'elapsed_ns':value['elapsed_ns'],'metrics':metrics,'stdout_sha256':sha(folder/'stdout.txt')});save(out/'summary.json',receipt)
    print(rep,len(receipt['runs']),flush=True)
for path,h in receipt['bindings'].items():assert sha(Path(path))==h
if not a.preflight:receipt['host_postflight']=postflight(out)
receipt['complete']=True;save(out/'summary.json',receipt)
if not a.preflight and not receipt['host_postflight']['quiet']:raise SystemExit('Quiet gate failed after certified batch scaling; attempt retained')
