"""Isolated positive-interval specialization after measured arithmetic profiling."""
import os,shutil,subprocess
from pathlib import Path
from campaign import EVIDENCE,DOCS,sha,save
base=EVIDENCE/'c1_ablations_v1/specialized';out=EVIDENCE/'c3_positive_intervals';out.mkdir(exist_ok=False)
shutil.copytree(base/'src',out/'src')
for name in ['Cargo.toml','Cargo.lock']:shutil.copy2(base/name,out/name)
p=out/'src/certified.rs';s=p.read_text();s=s.replace('.mul(','.chord_mul(').replace('.div_positive(','.chord_div_positive(');s+='\n'+(DOCS/'positive_intervals.rs').read_text();p.write_text(s)
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/ablations',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
cmd=['cargo','build','--locked','--offline','--release','--manifest-path',str(out/'Cargo.toml')]
with (out/'build.log').open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
receipt={'scope':'Only C1 multiplication/division call sites change. Original interval arithmetic remains the fallback and the bitwise reference. Production is unchanged.','command':cmd,'returncode':r.returncode,'sources':{str(p):sha(p) for p in (out/'src').glob('*.rs')}};save(out/'receipt.json',receipt);r.check_returncode();shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/stericx-chord-ablation',out/'probe');receipt['binary_sha256']=sha(out/'probe');save(out/'receipt.json',receipt)
