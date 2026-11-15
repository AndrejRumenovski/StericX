#!/usr/bin/env python3
"""Complete a frozen profile build, including local Python binding validation artifacts."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path)
    args=parser.parse_args()
    out=args.build.resolve()
    assert json.loads((out/'manifest.json').read_text())['complete']
    receipt_path=out/'complete_build.json'
    assert not receipt_path.exists()
    source=out/'source'
    for directory in ['examples','python','tests']:
        shutil.copytree(ROOT/directory,source/directory,ignore=shutil.ignore_patterns(
            'target','__pycache__','.pytest_cache','.venv','*.egg-info'))
    sources={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()}
    receipt={'sources':sources,'commands':[],'binaries':{},'complete':False}
    def save():
        receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    def run(label,argv,env):
        argv=list(map(str,argv)); start=time.monotonic()
        with (out/(label+'.log')).open('x') as log:
            proc=subprocess.run(argv,env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        receipt['commands'].append({'label':label,'argv':argv,'returncode':proc.returncode,
            'seconds':time.monotonic()-start,'environment':{k:env.get(k) for k in
             ['RUSTFLAGS','CARGO_TARGET_DIR','CARGO_BUILD_JOBS']},'log_sha256':sha(out/(label+'.log'))})
        save(); proc.check_returncode()
    env=dict(os.environ,CARGO_BUILD_JOBS='2')
    env.pop('RUSTFLAGS',None)
    for label,extra in [('native',[]),('diagnostic',['--features','profiling'])]:
        target=out/(label+'-target')
        run(label,['cargo','build','--locked','--offline','--release','--bins','--examples',
            '--manifest-path',source/'Cargo.toml',*extra],
            dict(env,RUSTFLAGS='-C debuginfo=1',CARGO_TARGET_DIR=str(target)))
        receipt['binaries'][label]={'path':str(target/'release/stericx'),'sha256':sha(target/'release/stericx')}
    run('extra_observers',['cargo','build','--locked','--offline','--release',
        '--manifest-path',out/'observer/Cargo.toml','--bins'],dict(env,CARGO_TARGET_DIR=str(out/'target')))
    for name in ['continuous','continuous_volume','accuracy','dump_density']:
        shutil.copy2(out/'target/release'/name,out/'bin'/name)
        receipt['binaries'][name]={'path':str(out/'bin'/name),'sha256':sha(out/'bin'/name)}
    python=ROOT/'.stericx/paton/python-build/bin/python'
    run('python_wheel',[python.parent/'maturin','build','--locked','--offline','--release',
        '--manifest-path',source/'python/Cargo.toml','--interpreter',python,'--out',out/'wheels'],
        dict(env,CARGO_TARGET_DIR=str(out/'python-target')))
    wheels=list((out/'wheels').glob('*.whl')); assert len(wheels)==1
    with zipfile.ZipFile(wheels[0]) as wheel:
        wheel.extractall(out/'python-installed')
    receipt['wheel']={'path':str(wheels[0]),'sha256':sha(wheels[0])}
    assert sources=={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()}
    receipt['complete']=True; save()
    print(json.dumps({'complete':True,'path':str(receipt_path)}))

if __name__=='__main__':main()
