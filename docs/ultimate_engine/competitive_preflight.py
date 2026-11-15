#!/usr/bin/env python3
"""Untimed ordered-output preflight at every requested worker count."""
import argparse,gzip,hashlib,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--build',type=Path,required=True);p.add_argument('--gate',type=Path,required=True)
p.add_argument('--output',type=Path,required=True);a=p.parse_args()
out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
gate=json.loads(a.gate.read_text());assert gate['complete'] and gate['passed']
paths=[r['path'] for r in gate['rows']]
receipt={'scope':'Untimed output/order admission; process elapsed time is not benchmark evidence',
    'gate_sha256':hashlib.sha256(a.gate.read_bytes()).hexdigest(),'runs':[]}
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
for threads,cpus in [(1,[2]),(2,[2,3]),(4,list(range(4))),(6,list(range(6))),(12,list(range(12)))]:
    for tool in ['native','morfeus','dbstep']:
        if tool=='native':command=[str(a.build.resolve()/'adapter/target/release/volume-batch'),*paths]
        else:command=[str(ROOT/'.stericx/paton/reference/venv/bin/python'),str(Path(__file__).with_name('competitive_worker.py')),
            '--tool',tool,'--workers',str(threads),*paths]
        command=['taskset','-c',','.join(map(str,cpus)),*command]
        log=out/f'{tool}_{threads}'
        with log.with_suffix('.stderr').open('w') as err:
            result=subprocess.run(command,env=dict(env,RAYON_NUM_THREADS=str(threads)),stdout=subprocess.PIPE,stderr=err)
        log.with_suffix('.stdout.gz').write_bytes(gzip.compress(result.stdout,mtime=0))
        assert result.returncode==0,(tool,threads)
        value=json.loads(result.stdout);rows=value if tool=='native' else value['rows']
        expected=[r[tool] for r in gate['rows']]
        assert rows==expected,(tool,threads)
        receipt['runs'].append({'tool':tool,'workers':threads,'command':command,'returncode':result.returncode,
            'exact':True,'stdout_sha256':hashlib.sha256(result.stdout).hexdigest(),'rows':len(rows)})
        (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(tool,threads,len(rows),flush=True)
receipt['complete']=True
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
