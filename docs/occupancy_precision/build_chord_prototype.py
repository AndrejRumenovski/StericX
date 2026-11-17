"""Freeze and compile the research chord certificate, retaining every version."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
from campaign import EVIDENCE,DOCS,sha,save

p=argparse.ArgumentParser(description=__doc__);p.add_argument('name');a=p.parse_args()
base=EVIDENCE/'baseline';out=EVIDENCE/a.name
out.mkdir(exist_ok=False);(out/'src').mkdir()
(out/'Cargo.toml').write_text((base/'observer/Cargo.toml').read_text().replace('stericx-ultimate-observer','stericx-chord-prototype'))
(out/'Cargo.lock').write_text((base/'observer/Cargo.lock').read_text().replace('stericx-ultimate-observer','stericx-chord-prototype'))
prefix=(base/'source/src/geometry/continuous_volume.rs').read_text().replace('//!','//')
(out/'src/observed.rs').write_text(prefix+'\n'+(DOCS/'chord_prototype.rs').read_text())
(out/'src/main.rs').write_text('''
pub use steric_x::geometry::Molecule;
pub struct BuriedVolumeCalculator;
#[derive(Debug)] pub struct DescriptorError(pub String);
impl std::fmt::Display for DescriptorError {
 fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result {f.write_str(&self.0)}
}
impl std::error::Error for DescriptorError {}
#[allow(dead_code)] mod observed {include!("observed.rs");}
fn main(){
 let requests:Vec<serde_json::Value>=serde_json::from_str(&std::fs::read_to_string(std::env::args().nth(1).unwrap()).unwrap()).unwrap();
 for request in requests {println!("{}",observed::chord_probe(&request));}
}
''')
target=Path('/tmp/stericx-occupancy-precision-target/chord')
env=dict(os.environ,CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--release','--offline','--locked','--manifest-path',str(out/'Cargo.toml')]
with (out/'build.log').open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
save(out/'receipt.json',{'source_sha256':sha(out/'src/observed.rs'),'prototype_sha256':sha(DOCS/'chord_prototype.rs'),'command':cmd,'returncode':r.returncode})
r.check_returncode();shutil.copy2(target/'release/stericx-chord-prototype',out/'probe')
print(out,flush=True)
