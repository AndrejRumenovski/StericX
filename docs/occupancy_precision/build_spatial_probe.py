"""Freeze spatial/parallel research variants without changing production."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from campaign import EVIDENCE,DOCS,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',required=True);p.add_argument('--output',required=True);a=p.parse_args()
base=EVIDENCE/a.build;out=EVIDENCE/a.output;out.mkdir(exist_ok=False);(out/'src').mkdir()
for name in ['Cargo.toml','Cargo.lock']:(out/name).write_text((base/'observer'/name).read_text().replace('stericx-ultimate-observer','stericx-spatial-probe'))
(out/'src/observed.rs').write_bytes((base/'source/src/geometry/buried_volume.rs').read_bytes()+b'\n'+(DOCS/'spatial_probe.rs').read_bytes())
(out/'src/main.rs').write_text('''pub use steric_x::{Molecule,covalent_radius};
#[macro_export] macro_rules! profile_scope { ($($arg:tt)*) => {}; }
#[allow(dead_code)] mod observed { include!("observed.rs"); }
fn main(){let q:Vec<serde_json::Value>=serde_json::from_str(&std::fs::read_to_string(std::env::args().nth(1).unwrap()).unwrap()).unwrap();for r in q {println!("{}",observed::spatial_probe(&r));}}
''')
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/candidates',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--release','--offline','--locked','--manifest-path',str(out/'Cargo.toml')]
with (out/'build.log').open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
save(out/'receipt.json',{'build_command':cmd,'returncode':r.returncode,'source_sha256':sha(out/'src/observed.rs'),'adapter_sha256':sha(DOCS/'spatial_probe.rs')})
r.check_returncode();shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/stericx-spatial-probe',out/'probe')
requests=json.loads((EVIDENCE/'baseline_probe/requests.json').read_text())
for q in requests:q['iterations']=20
requests+=[{**q,'id':q['id']+'/dense_grid','density':.0005,'iterations':5} for q in requests[-3:]]
save(out/'requests.json',requests);print(out)
