"""Source-bound documentation, Python binding/research, MSRV and cross checks."""
import argparse,importlib,json,os,subprocess,sys,time,zipfile
from pathlib import Path
from campaign import ROOT,EVIDENCE,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',required=True);p.add_argument('--output',required=True);a=p.parse_args()
build=EVIDENCE/a.build;source=build/'source';out=EVIDENCE/a.output;out.mkdir(exist_ok=False)
env=dict(os.environ,CARGO_BUILD_JOBS='2',RAYON_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RUSTDOCFLAGS='-D warnings',LC_ALL='C');env.pop('RUSTFLAGS',None)
receipt={'controller_sha256':sha(Path(__file__)),'source_build_sha256':sha(build/'build.json'),'runs':[],'complete':False}
def run(name,cmd,extra=None,cwd=ROOT):
    cmd=list(map(str,cmd));start=time.monotonic()
    with (out/(name+'.log')).open('x') as log:r=subprocess.run(cmd,cwd=cwd,env=dict(env,**(extra or {})),stdout=log,stderr=subprocess.STDOUT)
    row={'name':name,'command':cmd,'cwd':str(cwd),'environment':extra,'returncode':r.returncode,'seconds':time.monotonic()-start,'log_sha256':sha(out/(name+'.log'))}
    receipt['runs'].append(row);save(out/'summary.json',receipt);print(name,r.returncode,flush=True);return r.returncode==0
run('format',['cargo','fmt','--all','--check'])
run('rustdoc',['cargo','doc','--locked','--offline','--no-deps','--all-features','--manifest-path',source/'Cargo.toml'],{'CARGO_TARGET_DIR':'/tmp/stericx-occupancy-precision-target/doc'})
run('msrv',['cargo','+1.85.0','test','--locked','--offline','--all-targets','--all-features','--manifest-path',source/'Cargo.toml'],{'CARGO_TARGET_DIR':'/tmp/stericx-occupancy-precision-target/msrv'},source)
for target in ['aarch64-apple-darwin','x86_64-pc-windows-msvc']:
    run('check_'+target,['cargo','+1.98.1','check','--locked','--offline','--lib','--all-features','--target',target,'--manifest-path',source/'Cargo.toml'],{'CARGO_TARGET_DIR':'/tmp/stericx-occupancy-precision-target/cross'})
python=ROOT/'.stericx/paton/python-build/bin/python'
if run('wheel',[python.parent/'maturin','build','--locked','--offline','--release','--manifest-path',source/'python/Cargo.toml','--interpreter',python,'--out',out/'wheels'],{'CARGO_TARGET_DIR':'/tmp/stericx-occupancy-precision-target/python'}):
    wheel,=list((out/'wheels').glob('*.whl'))
    with zipfile.ZipFile(wheel) as f:f.extractall(out/'installed')
    receipt['wheel_sha256']=sha(wheel);save(out/'summary.json',receipt)
    run('python_binding',[python,'-m','pytest','-q',source/'python/tests'],{'PYTHONPATH':str(out/'installed')})
run('python_research',[ROOT/'.venv/bin/python','-c',"import sys,unittest;sys.path.append('.stericx/paton/reference/venv/lib/python3.12/site-packages');r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover('tests',pattern='test_*.py'));sys.exit(0 if r.wasSuccessful() and not r.skipped else 1)"])
run('dbstep_upstream',[ROOT/'.stericx/paton/reference/venv/bin/python','-m','pytest','-q','tests'],cwd=ROOT/'.stericx/paton/reference/DBSTEP')
run('archived_morfeus',[ROOT/'.venv/bin/python','docs/benchmarks/morfeus_current/verify.py'])
receipt['complete']=True;receipt['all_commands_succeeded']=all(r['returncode']==0 for r in receipt['runs']);save(out/'summary.json',receipt)
if not receipt['all_commands_succeeded']:raise SystemExit(1)
