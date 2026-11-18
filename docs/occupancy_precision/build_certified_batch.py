"""Build a batch-scaling observer against the unchanged final engine."""
import os,shutil,subprocess
from pathlib import Path
from campaign import EVIDENCE,DOCS,sha,save
base=EVIDENCE/'final_profile_v2';out=EVIDENCE/'certified_batch_build';out.mkdir(exist_ok=False)
shutil.copytree(base/'observer',out/'observer')
source=out/'observer/src/bin/certified_batch.rs';shutil.copy2(DOCS/'certified_batch.rs',source)
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/final_profile_v2',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--locked','--offline','--release','--manifest-path',str(out/'observer/Cargo.toml'),'--bin','certified_batch']
with (out/'build.log').open('x') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
receipt={'source_build_sha256':sha(base/'build.json'),'observer_sha256':sha(source),'command':cmd,'returncode':r.returncode};save(out/'receipt.json',receipt);r.check_returncode()
shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/certified_batch',out/'certified_batch');receipt['binary_sha256']=sha(out/'certified_batch');save(out/'receipt.json',receipt)
