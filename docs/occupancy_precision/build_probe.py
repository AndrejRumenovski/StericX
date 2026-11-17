"""Build an isolated, source-bound occupancy instrumentation adapter."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
from campaign import ROOT, EVIDENCE, DOCS, sha, save

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--build', required=True)
a=p.parse_args()
build=EVIDENCE/a.build
out=EVIDENCE/(a.build+'_probe')
out.mkdir(exist_ok=False)
(out/'src').mkdir()
manifest=(build/'observer/Cargo.toml').read_text().replace('stericx-ultimate-observer','stericx-occupancy-probe')
(out/'Cargo.toml').write_text(manifest)
lock=(build/'observer/Cargo.lock').read_text().replace('stericx-ultimate-observer','stericx-occupancy-probe')
(out/'Cargo.lock').write_text(lock)
prefix=(build/'source/src/geometry/buried_volume.rs').read_bytes()
append=(DOCS/'occupancy_probe.rs').read_bytes()
(out/'src/observed.rs').write_bytes(prefix+b'\n'+append)
(out/'src/main.rs').write_text('''
pub use steric_x::{Molecule,covalent_radius};
#[macro_export] macro_rules! profile_scope { ($($arg:tt)*) => {}; }
#[allow(dead_code)] mod observed { include!("observed.rs"); }
fn main() {
    let path=std::env::args().nth(1).unwrap();
    let requests: Vec<serde_json::Value>=serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    for request in requests { println!("{}",observed::occupancy_probe(&request)); }
}
''')
target=Path('/tmp/stericx-occupancy-precision-target')/a.build
env=dict(os.environ,CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='2')
env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--release','--offline','--locked','--manifest-path',str(out/'Cargo.toml')]
with (out/'build.log').open('x') as log:
    result=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
save(out/'receipt.json',{'source_prefix_sha256':sha(build/'source/src/geometry/buried_volume.rs'),
    'append_sha256':sha(DOCS/'occupancy_probe.rs'),'command':cmd,'returncode':result.returncode})
result.check_returncode()
shutil.copy2(target/'release/stericx-occupancy-probe',out/'probe')
requests=[]
fixtures=list(map(json.loads,(ROOT/'.stericx/ultimate_engine/baseline_observations/geometry/stdout.jsonl').read_text().splitlines()))
for fixture in fixtures[:12]:
    for orientation in fixture['dump']['orientations']:
        requests.append({'id':str(fixture['id'])+'/'+str(orientation['plane']),
                         'atoms':orientation['aligned_atoms'],'iterations':30})
import random
rng=random.Random(20261117)
for count in [1,4,16,64,256,1024]:
    requests.append({'id':f'synthetic_{count}', 'atoms':[
        {'position':[rng.uniform(-12,12) for _ in range(3)],'radius_squared':rng.uniform(.25,6.25)} for _ in range(count)],'iterations':30})
save(out/'requests.json',requests)
print(out,flush=True)
