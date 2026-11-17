"""Final sequential controlled measurements after complete scientific admission."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
from campaign import ROOT,EVIDENCE,DOCS,sha,save
sys.path.insert(0,str(ROOT/'docs/ultimate_engine'))
from host_activity import measure,postflight
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);p.add_argument('--start-at',default='cli_ab');p.add_argument('--gate-only',action='store_true');a=p.parse_args()
out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
def read(path):return json.loads((EVIDENCE/path).read_text())
assert read('final_exact_comparison/summary.json')['exact']
assert read('final_engineering/summary.json')['all_commands_succeeded']
assert read('final_arm/summary.json')['complete'] and not read('final_arm/summary.json')['failed_runs']
assert read('final_science_review/review.json')['legacy_reference_preservation_pass']
for path in ['final_cli_preflight/summary.json','accuracy_preflight/summary.json','accuracy_memory/summary.json','competitive_preflight/summary.json','competitive_warm_preflight/summary.json','competitive_radius_preflight/summary.json']:
    assert read(path)['complete'],path
for path in ['competitive_gate/summary.json','competitive_radius_gate/summary.json']:
    assert read(path)['passed'],path
python=ROOT/'.venv/bin/python';old=ROOT/'docs/ultimate_engine';final=EVIDENCE/'final_profile';baseline=EVIDENCE/'baseline'
env=dict(os.environ,RAYON_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');env.pop('STERICX_PROFILE_PATH',None)
commands=[
    ('cli_ab',[python,old/'final_paired_benchmark.py','--baseline',baseline,'--candidate',final,'--output',out/'cli_ab']),
    ('kernel_ab',[python,old/'compare_candidates.py','--baseline',baseline,'--candidate',final,'--output',out/'kernel_ab','--pairs','8','--scope','Final sequential quiet-window A/B, CPU 2; all native samples retained.']),
    ('kernel_profile',[python,old/'profile_campaign.py','--build',final,'--output',out/'kernel_profile','--phase','kernels']),
    ('cli_profile',[python,old/'profile_campaign.py','--build',final,'--output',out/'cli_profile','--phase','cli']),
    ('accuracy',[python,DOCS/'accuracy_campaign.py','timing','--output',out/'accuracy']),
    ('competitive',[python,old/'competitive_measure.py','timing','--gate',EVIDENCE/'competitive_gate/summary.json','--preflight',EVIDENCE/'competitive_preflight/summary.json','--native',final/'bin/volume_workflow','--output',out/'competitive']),
    ('radius_competitive',[python,DOCS/'radius_benchmark.py','timing','--output',out/'radius_competitive']),
]
assert a.start_at in [name for name,_ in commands]
commands=commands[[name for name,_ in commands].index(a.start_at):]
receipt={'controller_sha256':sha(Path(__file__)),'source_build_sha256':sha(final/'build.json'),'baseline_build_sha256':sha(baseline/'build.json'),'scope':'No builds or independent reference work overlaps this measurement window. Quiet endpoints use the prospectively frozen <=0.15 busy cores and <=0.05 I/O-wait cores gate. They do not guarantee zero interruptions between endpoints.','runs':[],'complete':False}
if a.gate_only:
    receipt['host_preflight']=measure();receipt['gate_only']=True;save(out/'receipt.json',receipt);print(json.dumps(receipt['host_preflight']));raise SystemExit(0 if receipt['host_preflight']['quiet'] else 2)
for name,command in commands:
    row={'name':name,'command':list(map(str,command)),'status':'pending'};receipt['runs'].append(row);save(out/'receipt.json',receipt)
    row['host_preflight']=measure();save(out/'receipt.json',receipt)
    if not row['host_preflight']['quiet']:
        row['status']='host_not_quiet_before_measurement';save(out/'receipt.json',receipt);raise SystemExit('Quiet-window confirmation needed before '+name)
    start=time.monotonic()
    with (out/(name+'.log')).open('x') as log:r=subprocess.run(row['command'],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
    row.update(returncode=r.returncode,seconds=time.monotonic()-start,log_sha256=sha(out/(name+'.log')),status='completed' if r.returncode==0 else 'failed');save(out/'receipt.json',receipt);print(name,r.returncode,flush=True);r.check_returncode()
    row['host_postflight']=postflight(out/name);save(out/'receipt.json',receipt)
    if not row['host_postflight']['quiet']:
        row['status']='host_not_quiet_after_measurement';save(out/'receipt.json',receipt);raise SystemExit('Attempt retained; quiet window failed after '+name)
receipt['complete']=True;save(out/'receipt.json',receipt)
