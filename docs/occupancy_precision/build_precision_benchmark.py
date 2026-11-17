"""Supplement frozen final source with matched precision observers."""
import argparse,os,shutil,subprocess
from pathlib import Path
from campaign import EVIDENCE,DOCS,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args()
base=EVIDENCE/'final_profile';out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
shutil.copytree(base/'observer',out/'observer');observer=out/'observer';shutil.copy2(DOCS/'precision_benchmark.rs',observer/'src/bin/precision.rs')
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/final_profile',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
receipt={'final_source_receipt_sha256':sha(base/'build.json'),'observer_sha256':sha(observer/'src/bin/precision.rs'),'commands':[],'binaries':{}}
for name,features in [('precision',[]),('precision_alloc',['--features','profiling'])]:
    cmd=['cargo','build','--locked','--offline','--release','--manifest-path',str(observer/'Cargo.toml'),'--bin','precision',*features]
    with (out/(name+'.log')).open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
    receipt['commands'].append({'command':cmd,'returncode':r.returncode});save(out/'receipt.json',receipt);r.check_returncode();shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/precision',out/name);receipt['binaries'][name]=sha(out/name)
save(out/'receipt.json',receipt)
