#!/usr/bin/env python3
"""Run compiled AArch64/musl tests under QEMU, including child CLI processes.

A private shell shim launches the preserved ARM CLI binary for guest execve;
the host's binfmt configuration is not changed. This is not a timing benchmark.
"""
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
base=ROOT/'.stericx/ultimate_engine'
out=base/'arm_validation';out.mkdir(exist_ok=False)
qemu=base/'tools/qemu-user/usr/bin/qemu-aarch64'
link=Path('/tmp/stericx-ultimate-qemu-aarch64')
if link.is_symlink():assert link.resolve()==qemu.resolve()
else:link.symlink_to(qemu.resolve())
target=base/'arm-target/aarch64-unknown-linux-musl/debug'
cli=target/'stericx';preserved=out/'stericx-arm'
shutil.copyfile(cli,preserved);preserved.chmod(0o755)
cli.write_text('#!/bin/sh\nexec '+shlex.quote(str(link))+' '+shlex.quote(str(preserved))+' "$@"\n');cli.chmod(0o755)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
receipt={'scope':'QEMU AArch64 functional tests; no hardware speed claim',
         'qemu':{'path':str(qemu),'sha256':sha(qemu)},'arm_cli_sha256':sha(preserved),
         'cli_shim_sha256':sha(cli),'runs':[]}
log=(base/'arm_build.log').read_text()
paths=re.findall(r'Executable .* \(([^\n]+)\)',log)
assert paths,'no compiled tests discovered'
env=dict(os.environ,RAYON_NUM_THREADS='2')
env.pop('STERICX_PROFILE_PATH',None)
for i,path in enumerate(paths):
    binary=ROOT/path
    command=[str(link),str(binary),'--test-threads=2']
    with (out/f'{i}.log').open('w') as stream:
        proc=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=300)
    text=(out/f'{i}.log').read_text()
    count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',text)))
    receipt['runs'].append({'command':command,'binary_sha256':sha(binary),'returncode':proc.returncode,
                            'passed_tests':count,'log':str(out/f'{i}.log')})
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(binary.name,proc.returncode,count,flush=True)
receipt['complete']=True
receipt['passed_tests']=sum(r['passed_tests'] for r in receipt['runs'])
receipt['failed_runs']=[r['command'][1] for r in receipt['runs'] if r['returncode']]
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
