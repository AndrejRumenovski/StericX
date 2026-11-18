"""Fixed accuracy/work ladder with separate scientific, native and allocation runs."""
import argparse,json,os,statistics,subprocess,sys
from pathlib import Path
from campaign import ROOT,EVIDENCE,sha,save
sys.path.insert(0,str(ROOT/'scripts'));import profile_stericx as harness
sys.path.insert(0,str(ROOT/'docs/ultimate_engine'))
from host_activity import measure,postflight
from benchmark_artifacts import compact_metrics
p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['preflight','timing','memory']);p.add_argument('--output',required=True);p.add_argument('--observer',default='precision_observer_v3');a=p.parse_args()
out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
engine=EVIDENCE/a.observer/('precision_alloc' if a.phase=='memory' else 'precision')
molecules=json.loads((EVIDENCE/'molecular_chord_requests.json').read_text())
fine={r['id'].replace('/target_0.001',''):r for r in map(json.loads,(EVIDENCE/'c1_extended_budget/observed.jsonl').read_text().splitlines()) if r['id'].endswith('/target_0.001')}
analytical=json.loads((EVIDENCE/'c1_coaxial/requests.json').read_text())
inputs=[{**q,'reference':[fine[q['id']]['lower'],fine[q['id']]['upper']],'reference_scope':'Separate 0.001 certificate from the same C1 algorithm; not an independent exact molecular value.'} for q in molecules]
inputs += [{**analytical[i],'reference_scope':'Independent exact rational polynomial integral times Machin-pi enclosure.'} for i in [0,1,2,24,54,105]]
configs=[('certified_'+str(t),{'kind':'certified','tolerance':t,'budget':1048575}) for t in [.1,.05,.01]]
configs += [('reference_'+str(b),{'kind':'reference','tolerance':.01,'budget':b}) for b in [262143,1048575]]
configs += [('midpoint_'+str(n),{'kind':'midpoint','subdivisions':n}) for n in [128,256,512]]
cases=[]
for q in inputs:
    for label,c in configs:
        case={**q,**c,'id':q['id']+'/'+label,'source_id':q['id'],'configuration':label,'warmups':1 if a.phase=='timing' else 0,'samples':1 if a.phase in ['timing','memory'] else 0,'iterations':1};cases.append(case)
save(out/'cases.json',cases)
receipt={'phase':a.phase,'scope':'Public API latency: one warmup, eight alternating configuration orders. All samples retained. Intrinsic certificates and observed error intervals are separate classes; midpoint timing excludes reference construction and cannot establish time to a certificate. 3D boxes and XY cells are different work units.','cases_sha256':sha(out/'cases.json'),'engine_sha256':sha(engine),'controller_sha256':sha(Path(__file__)),'runs':[],'complete':False}
env=dict(os.environ,RAYON_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1');env.pop('STERICX_PROFILE_PATH',None)
launcher=harness.build_launcher(out/'launcher') if a.phase in ['timing','memory'] else None
if a.phase=='timing':
    receipt['host_preflight']=measure();save(out/'summary.json',receipt)
    if not receipt['host_preflight']['quiet']:raise SystemExit('Quiet window required; no timing started.')
expected={}
for rep in range(8 if a.phase=='timing' else 1):
    for index in (range(len(cases)) if rep%2==0 else reversed(range(len(cases)))):
        q=cases[index];folder=out/f'case_{index:03d}_rep_{rep}';request=out/f'request_{index:03d}.json'
        if not request.exists():save(request,[q])
        cmd=[str(engine),str(request)]
        if launcher:
            metrics=compact_metrics(harness.measure(cmd,env,{2},folder,launcher),folder);data=json.loads((folder/'stdout.txt').read_text());assert metrics.get('exit_code',metrics.get('returncode',0))==0
        else:
            folder.mkdir();r=subprocess.run(['taskset','-c','2',*cmd],capture_output=True,env=env,check=True);(folder/'stdout.txt').write_bytes(r.stdout);(folder/'stderr.txt').write_bytes(r.stderr);data=json.loads(r.stdout);metrics=None
        result=data['output']['result'];science={k:v for k,v in result.items() if k!='elapsed_seconds'}
        assert science==expected.setdefault(q['id'],science),(q['id'],'nondeterministic')
        lo,hi=q['reference'];v=result['volume_angstrom3'];observed_error=[max(0.,lo-v,v-hi),max(abs(v-lo),abs(v-hi))]
        lower=result.get('lower_bound_angstrom3');upper=result.get('upper_bound_angstrom3');bound=result.get('absolute_error_bound_angstrom3')
        if lower is not None:
            assert lower<=v<=upper and max(lower,lo)<=min(upper,hi),(q['id'],'enclosure')
            assert result['target_met']==(bound<=q['tolerance'])
            if q['reference_scope'].startswith('Independent'):assert lower<=lo<=hi<=upper,(q['id'],'analytical containment')
        row={'id':q['id'],'source_id':q['source_id'],'configuration':q['configuration'],'rep':rep,'result':result,'metrics':metrics,'elapsed_ns':data['elapsed_ns'],'allocations':data['allocations'],'observed_error_interval':observed_error,'reference_scope':q['reference_scope']}
        receipt['runs'].append(row);save(out/'summary.json',receipt)
    print(a.phase,rep+1,len(receipt['runs']),flush=True)
receipt['complete']=True
if a.phase=='timing':receipt['host_postflight']=postflight(out)
save(out/'summary.json',receipt)
