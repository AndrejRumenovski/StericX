#!/usr/bin/env python3
"""Resume checks after the preserved Cargo archive traversal interruption."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from replay_dbstep_resolution_science import commands
base=ROOT/'.stericx/ultimate_engine'
out=base/'final_engineering_resume';out.mkdir()
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
         RAYON_NUM_THREADS='2',CARGO_BUILD_JOBS='2',RUSTDOCFLAGS='-D warnings',LC_ALL='C',TZ='UTC')
rows=commands('engineering',base/'final_scientific_build',base/'final_science')[2:]
rows[0]=('rustdoc',['cargo','doc','--locked','--offline','--no-deps','--all-features',
    '--manifest-path',str(base/'final_scientific_build/source/Cargo.toml')],ROOT)
receipt={'reason':'Original attempt retained; unchanged final scientific source isolated for rustdoc.',
         'controller_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runs':[]}
for name,argv,cwd in rows:
    start=time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        p=subprocess.run(argv,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
    row={'name':name,'argv':argv,'returncode':p.returncode,'seconds':time.monotonic()-start,
         'log_sha256':hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest()}
    receipt['runs'].append(row)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(row),flush=True)
receipt['complete']=True
receipt['all_commands_succeeded']=all(r['returncode']==0 for r in receipt['runs'])
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
