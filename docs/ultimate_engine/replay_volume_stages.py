#!/usr/bin/env python3
"""Run independent volume computation before its dependent comparisons."""
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from replay_dbstep_resolution_science import commands
base=ROOT/'.stericx/ultimate_engine'
p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['compute','compare']);a=p.parse_args()
science=base/'final_science';out=science/('volume_execution_'+a.phase);out.mkdir()
jobs=commands('volume',base/'final_scientific_build',science)
jobs=jobs[:1] if a.phase=='compute' else jobs[1:]
if a.phase=='compare':
    assert json.loads((science/'execution/references/journal.json').read_text())['complete']
    assert json.loads((science/'volume_execution_compute/receipt.json').read_text())['complete']
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS='2')
receipt={'controller_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'observations_manifest_sha256':hashlib.sha256((science/'observations/manifest.json').read_bytes()).hexdigest(),'runs':[],'complete':False}
for name,argv,cwd in jobs:
    if argv[:2]==['taskset','-c']:argv[2]='2,3'
    start=time.monotonic()
    with (out/(name+'.log')).open('w') as log:r=subprocess.run(argv,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
    row={'name':name,'argv':argv,'returncode':r.returncode,'seconds':time.monotonic()-start,
        'log_sha256':hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest()}
    receipt['runs'].append(row);(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(row),flush=True);r.check_returncode()
receipt['complete']=True;(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
