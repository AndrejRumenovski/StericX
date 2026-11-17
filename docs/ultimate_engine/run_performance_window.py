#!/usr/bin/env python3
"""Sequential final measurements after scientific admission; no builds or reference jobs overlap."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
from host_activity import measure, postflight
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/ultimate_engine'
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
def load(relative):return json.loads((BASE/relative).read_text())
assert load('final_exact_comparison/summary.json')['exact']
assert load('final_science/execution/references/journal.json')['complete']
assert load('final_science/volume_execution_compare/receipt.json')['complete']
assert load('final_science_review/review.json')['legacy_reference_preservation_pass']
assert load('final_engineering_resume/receipt.json')['all_commands_succeeded']
assert load('final_numerical_replay/summary.json')['complete']
assert load('paired_python_preflight_v2/summary.json')['complete']
assert load('competitive_gate/summary.json')['passed']
assert load('competitive_warm_preflight/summary.json')['complete']
assert load('final_cli_ab_preflight/summary.json')['complete']
python=ROOT/'.venv/bin/python';binding_python=ROOT/'.stericx/paton/python-build/bin/python'
docs=ROOT/'docs/ultimate_engine';candidate=BASE/'final_profile';baseline=BASE/'baseline'
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS='1')
env.pop('STERICX_PROFILE_PATH',None)
commands=[
 ('cli_ab',[python,docs/'final_paired_benchmark.py','--baseline',baseline,'--candidate',candidate,'--output',out/'cli_ab'],env),
 ('kernel_ab',[python,docs/'compare_candidates.py','--baseline','baseline','--candidate','final_profile','--output',str((out/'kernel_ab').relative_to(BASE)),'--pairs','8','--scope','Final sequential quiet-window kernel A/B; CPU 2; parent receipt records host admission'],env),
 ('kernel_profile',[python,docs/'profile_campaign.py','--build',candidate,'--output',out/'kernel_profile','--phase','kernels'],env),
 ('cli_profile',[python,docs/'profile_campaign.py','--build',candidate,'--output',out/'cli_profile','--phase','cli'],env),
 ('python_ab',[python,docs/'paired_python_benchmark.py','--output',out/'python_ab'],env),
 ('python_native_bridge',[binding_python,ROOT/'scripts/benchmark_paton_python.py','--native',candidate/'native-target/release/examples/paton_ensemble_benchmark',
    '--input',ROOT/'python/examples/data/ligand_401_ensemble.xyz','--options',ROOT/'python/examples/data/ligand_401_options.json','--iterations','40','--output',out/'python_native_bridge'],dict(env,PYTHONPATH=str(candidate/'python-installed'))),
 ('accuracy_cost',[python,docs/'measure_accuracy_cost.py','--output',out/'accuracy_cost'],env),
 ('competitive',[python,docs/'competitive_measure.py','timing','--output',out/'competitive'],env),
]
receipt={'scope':'Final sequential measurement window; engineering/science already passed in declared scopes',
    'controller_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'helper_hashes':{str(docs/name):hashlib.sha256((docs/name).read_bytes()).hexdigest() for name in ['host_activity.py','benchmark_artifacts.py']},
    'runs':[],'complete':False}
def save():(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
for name,command,child_env in commands:
    command=list(map(str,command));start=time.monotonic()
    row={'name':name,'command':command,'status':'started'};receipt['runs'].append(row);save()
    row['host_preflight']=measure();save()
    if not row['host_preflight']['quiet']:
        row['status']='host_not_quiet_before_measurement';save()
        raise SystemExit('Measurement stopped before '+name+': quiet window required.')
    with (out/(name+'.log')).open('w') as log:r=subprocess.run(command,cwd=ROOT,env=child_env,stdout=log,stderr=subprocess.STDOUT)
    row.update(returncode=r.returncode,seconds=time.monotonic()-start,status='completed' if r.returncode==0 else 'stopped',
               log_sha256=hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest());save()
    print(json.dumps(row),flush=True)
    r.check_returncode()
    row['host_postflight']=postflight(out/name);save()
    if not row['host_postflight']['quiet']:
        row['status']='host_not_quiet_after_measurement';save()
        raise SystemExit('Retain '+name+' as attempted measurement; host activity requires review.')
receipt['complete']=True;save()
