#!/usr/bin/env python3
"""Alternate original/final wheels in identical Python environments and CPU affinities."""
import argparse,json,os,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import profile_stericx as h
from host_activity import measure
BASE=ROOT/'.stericx/ultimate_engine'
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
p.add_argument('--preflight-only',action='store_true');a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
receipt={'scope':'Same interpreter; source-bound original/final wheels; imports/input setup excluded from internal timers',
    'preflight_only':a.preflight_only,'runs':[],'complete':False}
def save():(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
if not a.preflight_only:
    receipt['host_preflight']=measure();save()
    if not receipt['host_preflight']['quiet']:raise SystemExit('No timing: quiet window required.')
launcher=h.build_launcher(out/'launcher');env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS='1')
env.pop('STERICX_PROFILE_PATH',None)
python=ROOT/'.stericx/paton/python-build/bin/python';worker=Path(__file__).with_name('python_measure_worker.py')
variants={'baseline':BASE/'baseline/python-installed','candidate':BASE/'final_profile/python-installed'}
receipt['controller_sha256']=h.digest(Path(__file__));receipt['worker_sha256']=h.digest(worker)
shutil.copy2(Path(__file__),out/Path(__file__).name);shutil.copy2(worker,out/worker.name)
cases=[(m,1,n) for m,n in [('bridge',40),('public',40),('json_decode',500),('typed_construction',500),
    ('numpy_export',500),('options_encode',1000),('record_export',200),('per_file_56',3)]]
cases += [('batch_56',t,3) for t in [1,2,4,6,12]]
expected={}
for mode,threads,count in cases:
    affinity=[2] if threads==1 else [2,3] if threads==2 else list(range(threads))
    for pair in ([0] if a.preflight_only else range(8)):
        sides=['baseline','candidate'] if pair%2==0 else ['candidate','baseline']
        for side in sides:
            folder=out/f'{mode}_{threads}_{pair}_{side}'
            command=[str(python),str(worker),mode,str(1 if a.preflight_only else count),'--threads',str(threads)]
            metrics=h.measure(command,dict(env,PYTHONPATH=str(variants[side]),RAYON_NUM_THREADS=str(threads)),set(affinity),folder,launcher)
            row={'mode':mode,'threads':threads,'pair':pair,'side':side,'metrics':metrics};receipt['runs'].append(row);save()
            assert metrics.get('exit_code',metrics.get('returncode',0))==0,(mode,threads,side)
            value=json.loads((folder/'stdout.txt').read_text());row['kernel']=value
            assert value['result_sha256']==expected.setdefault((mode,threads),value['result_sha256']),(mode,side)
            row['exact']=True;save()
if not a.preflight_only:receipt['host_postflight']=measure()
receipt['complete']=True;save()
