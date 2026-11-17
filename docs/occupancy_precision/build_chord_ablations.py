"""Source-bound ablations: smooth remainder and inherited geometric candidates."""
import argparse,os,shutil,subprocess
from pathlib import Path
from campaign import EVIDENCE,DOCS,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',required=True);p.add_argument('--output',required=True);a=p.parse_args()
base=EVIDENCE/a.build;out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/ablations',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
receipt={'variants':{}}
for variant in ['specialized','all_atoms','range_only']:
    folder=out/variant;(folder/'src').mkdir(parents=True)
    for name in ['Cargo.toml','Cargo.lock']:
        (folder/name).write_text((base/'observer'/name).read_text().replace('stericx-ultimate-observer','stericx-chord-ablation'))
    parent=(base/'source/src/geometry/continuous_volume.rs').read_text().replace('//!','//')
    child=(base/'source/src/geometry/continuous_volume/certified.rs').read_text()
    if variant=='all_atoms':
        needle='&parent.candidates,';assert child.count(needle)==1;child=child.replace(needle,'&(1..spheres.len()).collect::<Vec<_>>(),')
    if variant=='range_only':
        needle='let Some(endpoints) = chord_smooth(&spans, spheres) else {';assert child.count(needle)==1
        child=child.replace(needle,'let Some(endpoints) = None::<Vec<(ChordEndpoint, ChordEndpoint)>> else {')
    (folder/'src/observed.rs').write_text(parent);(folder/'src/certified.rs').write_text(child)
    observer=(DOCS/'certified_observer.rs').read_text().replace('//!','//').replace('use steric_x::geometry::{','use observed::{')
    prefix='''pub use steric_x::geometry::Molecule;
pub struct BuriedVolumeCalculator;
#[derive(Debug)] pub struct DescriptorError(pub String);
impl std::fmt::Display for DescriptorError {fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result {f.write_str(&self.0)}}
impl std::error::Error for DescriptorError {}
#[allow(dead_code,unused_imports)] mod observed {include!("observed.rs");}
'''
    (folder/'src/main.rs').write_text(prefix+observer)
    cmd=['cargo','build','--release','--offline','--locked','--manifest-path',str(folder/'Cargo.toml')]
    with (folder/'build.log').open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
    receipt['variants'][variant]={'command':cmd,'returncode':r.returncode,'sources':{str(p):sha(p) for p in (folder/'src').glob('*.rs')}};save(out/'receipt.json',receipt)
    r.check_returncode();shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/stericx-chord-ablation',folder/'probe');receipt['variants'][variant]['binary_sha256']=sha(folder/'probe');save(out/'receipt.json',receipt)
print(out)
