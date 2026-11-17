"""Fresh independent volume stages with explicit dependency receipts."""
import argparse,os,subprocess,sys,time
from pathlib import Path
from campaign import ROOT,EVIDENCE,sha,save
sys.path.insert(0,str(ROOT/'scripts'))
from replay_dbstep_resolution_science import commands
p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['compute','compare']);a=p.parse_args()
science=EVIDENCE/'final_science';out=science/('volume_execution_'+a.phase);out.mkdir(exist_ok=False)
jobs=commands('volume',EVIDENCE/'final_scientific_build',science)
if a.phase=='compare':
    import json
    assert json.loads((science/'execution/references/journal.json').read_text())['complete']
    assert json.loads((science/'volume_execution_compute/receipt.json').read_text())['complete']
jobs=jobs[:1] if a.phase=='compute' else jobs[1:]
receipt={'controller_sha256':sha(Path(__file__)),'observation_manifest_sha256':sha(science/'observations/manifest.json'),'runs':[],'complete':False}
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS='2')
for name,cmd,cwd in jobs:
    if cmd[:2]==['taskset','-c']:cmd[2]='2,3'
    start=time.monotonic()
    with (out/(name+'.log')).open('x') as f:r=subprocess.run(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT)
    receipt['runs'].append({'name':name,'command':cmd,'returncode':r.returncode,'seconds':time.monotonic()-start,'log_sha256':sha(out/(name+'.log'))});save(out/'receipt.json',receipt);r.check_returncode();print(name,r.returncode,flush=True)
receipt['complete']=True;save(out/'receipt.json',receipt)
