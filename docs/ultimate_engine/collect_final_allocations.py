#!/usr/bin/env python3
"""Replay baseline instrumented kernel calls for memory evidence, not native speed."""
import hashlib,json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import profile_stericx as h
base=ROOT/'.stericx/ultimate_engine'
out=base/'final_allocations';out.mkdir()
prior=json.loads((base/'baseline_kernels/journal.json').read_text())
binary=base/'final_profile/bin/observer_alloc'
launcher=h.build_launcher(out/'launcher')
env=dict(os.environ,RAYON_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LC_ALL='C')
env.pop('STERICX_PROFILE_PATH',None)
receipt={'scope':'Instrumented allocation and process high-water measurements; elapsed time is not native performance',
    'binary_sha256':h.digest(binary),'baseline_journal_sha256':h.digest(base/'baseline_kernels/journal.json'),'runs':[]}
for old in prior['runs']:
    if not old['name'].endswith('/3') or old['name'].startswith('batch/'):continue
    command=list(old['argv']);command[0]=str(binary)
    name=old['name'].split('/')[0];folder=out/name
    metrics=h.measure(command,env,{2},folder,launcher)
    assert metrics.get('exit_code',metrics.get('returncode',0))==0
    value=json.loads((folder/'stdout.txt').read_text())
    fingerprint=hashlib.sha256(json.dumps(value['result'],sort_keys=True).encode()).hexdigest()
    assert fingerprint==old['result_sha256'],name
    receipt['runs'].append({'mode':name,'metrics':metrics,'kernel':{k:v for k,v in value.items() if k!='result'},
        'baseline_kernel':old['kernel'],'result_sha256':fingerprint,'exact':True})
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n');print(name,flush=True)
receipt['complete']=True
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
