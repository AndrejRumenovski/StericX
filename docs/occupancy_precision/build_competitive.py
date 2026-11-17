"""Build unchanged common-grid external-comparison adapters against final source."""
import argparse,os,re,shutil,subprocess
from pathlib import Path
from campaign import ROOT,EVIDENCE,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',required=True);p.add_argument('--output',required=True);a=p.parse_args()
build=EVIDENCE/a.build;out=EVIDENCE/a.output;out.mkdir(exist_ok=False);adapter=out/'adapter'
old=ROOT/'.stericx/ultimate_engine/competitive_build_pinned/adapter'
shutil.copytree(old,adapter,ignore=shutil.ignore_patterns('target'))
manifest=adapter/'Cargo.toml';s=manifest.read_text();s=re.sub(r'steric_x\s*=\s*\{\s*path\s*=\s*"[^"]+"\s*\}',f'steric_x = {{ path = "{build / "source"}" }}',s);manifest.write_text(s)
append=ROOT/'docs/benchmarks/morfeus_current/revision_policy/grid_append.rs'
(adapter/'policy_copied_volume.rs').write_bytes((build/'source/src/geometry/buried_volume.rs').read_bytes()+append.read_bytes())
target=Path('/tmp/stericx-occupancy-precision-target/competitive')
env=dict(os.environ,CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--locked','--offline','--release','--manifest-path',str(manifest)]
with (out/'build.log').open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
receipt={'source_build_sha256':sha(build/'build.json'),'sources':{str(p):sha(p) for p in adapter.rglob('*') if p.is_file()},'command':cmd,'returncode':r.returncode,'binaries':{}};save(out/'receipt.json',receipt);r.check_returncode()
(adapter/'target/release').mkdir(parents=True)
for name in ['volume-batch','policy-grid-observer']:
    shutil.copy2(target/'release'/name,adapter/'target/release'/name);receipt['binaries'][name]=sha(adapter/'target/release'/name)
save(out/'receipt.json',receipt);print(out)
