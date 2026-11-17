"""Freeze a separate three-radius common-lattice admission experiment."""
import os,shutil,subprocess
from pathlib import Path
from campaign import ROOT,EVIDENCE,DOCS,sha,save
out=EVIDENCE/'competitive_radius_build';out.mkdir(exist_ok=False);adapter=out/'adapter'
shutil.copytree(EVIDENCE/'competitive_build/adapter',adapter,ignore=shutil.ignore_patterns('target'))
source=adapter/'policy_copied_volume.rs';s=source.read_text();needle='((f64::from(x)+3.5)*31.0/7.0)';assert s.count(needle)==1;s=s.replace(needle,'((f64::from(x)+f64::from(cfg.sphere_radius))*31.0/(2.0*f64::from(cfg.sphere_radius)))');source.write_text(s)
shutil.copy2(DOCS/'radius_workflow.rs',adapter/'src/bin/radius_workflow.rs')
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/competitive',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--locked','--offline','--release','--manifest-path',str(adapter/'Cargo.toml'),'--bin','radius_workflow']
with (out/'build.log').open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
receipt={'scope':'Admission experiment: R=2.5,3.5,4.5 and density=0.01*(R/3.5)^3; 32-point side common-lattice profile, distinct from constant-density profile. Exact old policy equations with radius parameterized; reject any failed mask/reduction gate.','source_hashes':{str(p):sha(p) for p in adapter.rglob('*') if p.is_file()},'command':cmd,'returncode':r.returncode};save(out/'receipt.json',receipt);r.check_returncode();shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/radius_workflow',out/'radius_workflow');receipt['binary_sha256']=sha(out/'radius_workflow');save(out/'receipt.json',receipt)
