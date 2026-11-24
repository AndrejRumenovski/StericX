"""Fresh, logged engineering checks against a frozen scientific source tree."""
import argparse,os,subprocess,time,json,hashlib,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path);a=p.parse_args();source=a.source.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
env=dict(os.environ,CARGO_BUILD_JOBS='2',RAYON_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RUSTDOCFLAGS='-D warnings',LC_ALL='C');env.pop('RUSTFLAGS',None)
receipt=dict(source=str(source),source_hashes={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()},runs=[],complete=False)
def save():(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
def run(name,cmd,extra=None,cwd=ROOT):
 cmd=list(map(str,cmd));start=time.monotonic()
 with (out/(name+'.log')).open('x') as log:r=subprocess.run(cmd,cwd=cwd,env=dict(env,**(extra or {})),stdout=log,stderr=subprocess.STDOUT)
 receipt['runs'].append(dict(name=name,command=cmd,cwd=str(cwd),environment=extra,returncode=r.returncode,seconds=time.monotonic()-start,log_sha256=sha(out/(name+'.log'))));save();print(name,r.returncode,flush=True);return r.returncode==0
run('format',['cargo','fmt','--all','--check'])
for name,args in [('clippy',['clippy','--all-targets','--all-features','--','-D','warnings']),('rustdoc',['doc','--no-deps','--all-features']),('rust_tests',['test','--all-targets','--all-features'])]:
 run(name,['cargo',*args],{'CARGO_TARGET_DIR':'/tmp/stericx-experimental-engineering/modern'},source)
run('msrv',['cargo','+1.85.0','test','--locked','--offline','--all-targets','--all-features'],{'CARGO_TARGET_DIR':'/tmp/stericx-experimental-engineering/msrv'},source)
for target in ['aarch64-apple-darwin','x86_64-pc-windows-msvc']:
 run('check_'+target,['cargo','+1.98.1','check','--locked','--offline','--lib','--all-features','--target',target],{'CARGO_TARGET_DIR':'/tmp/stericx-experimental-engineering/cross'},source)
python=ROOT/'.stericx/paton/python-build/bin/python'
# Scientific replay snapshot does not include the binding; copy it unchanged and hash it.
import shutil
if not (source/'python').exists():shutil.copytree(ROOT/'python',source/'python',ignore=shutil.ignore_patterns('target','__pycache__','.pytest_cache'))
receipt['binding_source_hashes']={str(p.relative_to(source)):sha(p) for p in (source/'python').rglob('*') if p.is_file()};save()
if run('wheel',[python.parent/'maturin','build','--locked','--offline','--release','--manifest-path',source/'python/Cargo.toml','--interpreter',python,'--out',out/'wheels'],{'CARGO_TARGET_DIR':'/tmp/stericx-experimental-engineering/python'}):
 wheel,=list((out/'wheels').glob('*.whl'))
 with zipfile.ZipFile(wheel) as f:f.extractall(out/'installed')
 receipt['wheel_sha256']=sha(wheel);save();run('python_binding',[python,'-m','pytest','-q',source/'python/tests'],{'PYTHONPATH':str(out/'installed')})
run('python_research',[ROOT/'.venv/bin/python','-c',"import sys,unittest;sys.path.append('.stericx/paton/reference/venv/lib/python3.12/site-packages');r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover('tests',pattern='test_*.py'));sys.exit(0 if r.wasSuccessful() and not r.skipped else 1)"])
run('dbstep_upstream',[ROOT/'.stericx/paton/reference/venv/bin/python','-m','pytest','-q','tests'],cwd=ROOT/'.stericx/paton/reference/DBSTEP')
run('archived_morfeus',[ROOT/'.venv/bin/python','docs/benchmarks/morfeus_current/verify.py'])
receipt['complete']=True;receipt['all_commands_succeeded']=all(r['returncode']==0 for r in receipt['runs']);save()
raise SystemExit(0 if receipt['all_commands_succeeded'] else 1)
