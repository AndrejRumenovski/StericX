"""Repeat invalid space-containing preload-path probes with a whitespace-free tool alias."""
from pathlib import Path
import subprocess,os,json,hashlib,time
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/experimental_validation';OUT=BASE/'valgrind_v2';INPUT=BASE/'performance_v1';BIN=Path('/tmp/stericx-experimental-profile/release')
OUT.mkdir(exist_ok=False);lib=Path('/tmp/stericx-experimental-valgrind-lib');prefix=ROOT/'.stericx/profiling/tools/valgrind'
if not lib.exists():lib.symlink_to(prefix/'usr/libexec/valgrind',target_is_directory=True)
assert lib.resolve()==(prefix/'usr/libexec/valgrind').resolve()
env=dict(os.environ,VALGRIND_LIB=str(lib));vgb=prefix/'usr/bin/valgrind.bin';rows=[]
for label,binary,request in [('old_chain',Path('/tmp/stericx-experimental-baseline-target/release/certified'),'one_chain.json'),('new_chain',BIN/'geometric','one_chain.json'),('regions',BIN/'regions','one_region.json'),('profile',BIN/'profiles','one_profile.json')]:
 for tool in (['callgrind','memcheck'] if label in ['regions','profile'] else ['callgrind']):
  name=tool+'_'+label;cmd=[str(vgb),'--tool='+tool,*(['--callgrind-out-file='+str(OUT/(name+'.out'))] if tool=='callgrind' else ['--error-exitcode=99','--leak-check=full','--errors-for-leak-kinds=definite']),str(binary),str(INPUT/request)];start=time.monotonic()
  with (OUT/(name+'.stdout')).open('wb') as stdout,(OUT/(name+'.stderr')).open('wb') as stderr:r=subprocess.run(cmd,env=env,stdout=stdout,stderr=stderr)
  text=(OUT/(name+'.stderr')).read_text();row=dict(name=name,command=cmd,returncode=r.returncode,elapsed_seconds=time.monotonic()-start,preload_failed='cannot be preloaded' in text,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),input_sha256=hashlib.sha256((INPUT/request).read_bytes()).hexdigest());rows.append(row);(OUT/'summary.json').write_text(json.dumps(dict(rows=rows,complete=False),indent=2)+'\n');print(name,r.returncode,flush=True)
  assert r.returncode==0 and not row['preload_failed'],row
(OUT/'summary.json').write_text(json.dumps(dict(rows=rows,complete=True,explanation='First probe used a VALGRIND_LIB containing a space; the dynamic loader splits LD_PRELOAD on spaces. Its memcheck/instruction evidence is rejected. This rerun uses a whitespace-free symlink to the identical tools; no suppressions added.'),indent=2)+'\n')
