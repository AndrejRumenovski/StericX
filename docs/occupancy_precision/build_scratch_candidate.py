"""Isolated scratch reuse on top of the positive-interval candidate."""
import os,re,shutil,subprocess
from pathlib import Path
from campaign import EVIDENCE,DOCS,sha,save
base=EVIDENCE/'c3_positive_intervals';out=EVIDENCE/'c4_scratch';out.mkdir(exist_ok=False);shutil.copytree(base/'src',out/'src')
for name in ['Cargo.toml','Cargo.lock']:shutil.copy2(base/name,out/name)
p=out/'src/certified.rs';s=p.read_text()
def replace(old,new,count=1):
    global s
    assert s.count(old)==count,(old,s.count(old));s=s.replace(old,new)
replace('let mut counters = ChordCounters {','let mut scratch = ChordScratch::default();\n    let mut counters = ChordCounters {')
replace('&mut counters,\n    )?;','&mut counters,\n        &mut scratch,\n    )?;')
replace('&mut counters,\n            )?;','&mut counters,\n                &mut scratch,\n            )?;')
replace('counters: &mut ChordCounters,\n) -> Result<ChordCell, DescriptorError>', 'counters: &mut ChordCounters,\n    scratch: &mut ChordScratch,\n) -> Result<ChordCell, DescriptorError>')
replace('fn chord_length(mut ranges: Vec<(f64, f64)>, upper: bool)', 'fn chord_length(ranges: &mut Vec<(f64, f64)>, upper: bool)')
start=s.index('fn chord_smooth(');end=s.index('fn chord_endpoint_value(',start);s=s[:start]+(DOCS/'scratch_smooth.rs').read_text()+'\n'+s[end:]
replace('let mut spans = Vec::with_capacity(indices.len());','let spans = &mut scratch.spans;\n    spans.clear();')
start=s.index('    let lower = chord_length(');end=s.index('    cell.bounds = Interval',start)
s=s[:start]+'''    scratch.ranges.clear();
    scratch.ranges.extend(spans.iter().filter(|s|s.regular).map(|s|(s.low.hi,s.high.lo)));
    let lower = chord_length(&mut scratch.ranges, false);
    scratch.ranges.clear();
    scratch.ranges.extend(spans.iter().map(|s|(s.low.lo,s.high.hi)));
    let upper = chord_length(&mut scratch.ranges, true);
'''+s[end:]
replace('let Some(endpoints) = chord_smooth(&spans, spheres) else {','if !chord_smooth(spans, spheres, &mut scratch.order, &mut scratch.components, &mut scratch.endpoints) {')
replace('for (a, b) in endpoints {','for &(a, b) in &scratch.endpoints {')
s+='\n'+(DOCS/'positive_interval_tests.rs').read_text();p.write_text(s)
env=dict(os.environ,CARGO_TARGET_DIR='/tmp/stericx-occupancy-precision-target/ablations',CARGO_BUILD_JOBS='2');env.pop('RUSTFLAGS',None)
receipt={'scope':'Only reusable per-call scratch storage and equivalent component index ranges change; every floating-point operation and iteration ordering is preserved.','sources':{str(p):sha(p) for p in (out/'src').glob('*.rs')},'commands':[]}
for name,args in [('build',['build','--release']),('tests',['test','--release'])]:
    cmd=['cargo',*args,'--locked','--offline','--manifest-path',str(out/'Cargo.toml')]
    with (out/(name+'.log')).open('x') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
    receipt['commands'].append({'command':cmd,'returncode':r.returncode});save(out/'receipt.json',receipt);r.check_returncode()
shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/stericx-chord-ablation',out/'probe');receipt['binary_sha256']=sha(out/'probe');save(out/'receipt.json',receipt)
