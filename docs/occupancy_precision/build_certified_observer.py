"""Build and preserve a source-bound public certified API observer."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
from campaign import build,EVIDENCE,DOCS,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('name');a=p.parse_args()
build(a.name,quick=True)
out=EVIDENCE/a.name
observer=out/'observer'
source=observer/'src/bin/certified.rs'
shutil.copy2(DOCS/'certified_observer.rs',source)
target=Path('/tmp/stericx-occupancy-precision-target/candidates')
env=dict(os.environ,CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
receipt={'observer_sha256':sha(source),'commands':[],'binaries':{}}
for label,features in [('certified',[]),('certified_alloc',['--features','profiling'])]:
    cmd=['cargo','build','--release','--offline','--locked','--manifest-path',str(observer/'Cargo.toml'),'--bin','certified',*features]
    with (out/(label+'_build.log')).open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
    receipt['commands'].append({'argv':cmd,'returncode':r.returncode});save(out/'certified_build.json',receipt)
    r.check_returncode();shutil.copy2(target/'release/certified',out/'bin'/label)
    receipt['binaries'][label]=sha(out/'bin'/label)
save(out/'certified_build.json',receipt)
print(out,flush=True)
