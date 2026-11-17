"""Fresh scalar AArch64/musl functional replay under retained local QEMU."""
import argparse,json,os,re,shlex,shutil,subprocess
from pathlib import Path
from campaign import ROOT,EVIDENCE,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',required=True);p.add_argument('--output',required=True);a=p.parse_args()
build=EVIDENCE/a.build;source=build/'source';out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
target=Path('/tmp/stericx-occupancy-precision-target/arm')
sysroot=Path(subprocess.check_output(['rustc','+1.85.0','--print','sysroot'],text=True).strip())
linker=sysroot/'lib/rustlib/x86_64-unknown-linux-gnu/bin/rust-lld'
qemu=ROOT/'.stericx/ultimate_engine/tools/qemu-user/usr/bin/qemu-aarch64'
env=dict(os.environ,CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='2',CARGO_TARGET_AARCH64_UNKNOWN_LINUX_MUSL_LINKER=str(linker),RUSTFLAGS='-C link-self-contained=yes',RAYON_NUM_THREADS='2')
cmd=['cargo','+1.85.0','test','--locked','--offline','--all-targets','--all-features','--no-run','--message-format=json','--target','aarch64-unknown-linux-musl','--manifest-path',str(source/'Cargo.toml')]
with (out/'build.jsonl').open('x') as f,(out/'build.log').open('x') as err:r=subprocess.run(cmd,env=env,cwd=source,stdout=f,stderr=err)
receipt={'scope':'QEMU AArch64 functional tests; no physical-hardware timing claim. CLI children use a private executable shim; host binfmt is unchanged.','build_command':cmd,'build_returncode':r.returncode,'qemu_sha256':sha(qemu),'linker_sha256':sha(linker),'source_receipt_sha256':sha(build/'build.json'),'runs':[],'complete':False};save(out/'summary.json',receipt);r.check_returncode()
cli=target/'aarch64-unknown-linux-musl/debug/stericx';preserved=out/'stericx-arm';shutil.copy2(cli,preserved)
cli.write_text('#!/bin/sh\nexec '+shlex.quote(str(qemu))+' '+shlex.quote(str(preserved))+' "$@"\n');cli.chmod(0o755)
paths=[]
for line in (out/'build.jsonl').open():
    row=json.loads(line)
    if row.get('reason')=='compiler-artifact' and row.get('executable') and row.get('profile',{}).get('test'):paths.append(Path(row['executable']))
assert paths
for i,binary in enumerate(paths):
    cmd=[str(qemu),str(binary),'--test-threads=2']
    with (out/f'test_{i}.log').open('x') as f:r=subprocess.run(cmd,cwd=source,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=300)
    log=out/f'test_{i}.log';count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',log.read_text())))
    receipt['runs'].append({'command':cmd,'returncode':r.returncode,'passed_tests':count,'binary_sha256':sha(binary),'log_sha256':sha(log)});save(out/'summary.json',receipt);print(binary.name,r.returncode,count,flush=True)
receipt['complete']=True;receipt['passed_tests']=sum(r['passed_tests'] for r in receipt['runs']);receipt['failed_runs']=[r for r in receipt['runs'] if r['returncode']];save(out/'summary.json',receipt)
if receipt['failed_runs']:raise SystemExit(1)
