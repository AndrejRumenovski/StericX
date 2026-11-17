#!/usr/bin/env python3
"""Prospectively gated, accuracy-matched twelve-molecule Sterimol workflows."""
import argparse
import hashlib
import json
import os
from decimal import Decimal, localcontext
from pathlib import Path
import subprocess
import sys
from host_activity import measure, postflight

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / '.stericx/ultimate_engine'
BUILD = BASE / 'sterimol_competitive_build'
HERE = Path(__file__).resolve().parent
METHODS = ['sampled361', 'continuous', 'morfeus', 'dbstep']

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def command(method, iterations):
    if method in METHODS[:2]:
        return [str(BUILD/'sterimol_workflow'), str(BUILD/'cases.json'), method, str(iterations)]
    return [str(BASE.parent/'paton/reference/venv/bin/python'),
            str(HERE/'sterimol_reference_worker.py'), str(BUILD/'cases.json'), method, str(iterations)]

def independent_references():
    previous = {r['id']: r['reference'] for r in json.loads((HERE/'n1_accuracy_measurements.json').read_text())['rows']}
    refs = {}
    with localcontext() as ctx:
        ctx.prec = 80
        for c in json.loads((BUILD/'cases.json').read_text()):
            convert = lambda x: Decimal.from_float(float(x))
            d = list(map(convert, c['dummy']))
            xyz = [[convert(x)-y for x,y in zip(a['position'],d)] for a in c['atoms']]
            axis = xyz[c['base']]
            norm = sum(x*x for x in axis).sqrt()
            axis = [x/norm for x in axis]
            lengths, widths = [], []
            for a, v in zip(c['atoms'], xyz):
                r = convert(a['radius'])
                z = sum(x*y for x,y in zip(v,axis))
                # Direct transverse vector avoids subtracting nearly equal squared lengths.
                transverse = [x-z*y for x,y in zip(v,axis)]
                lengths.append(z+r)
                widths.append(sum(x*x for x in transverse).sqrt()+r)
            b = previous[c['id']]
            assert b['converged']
            refs[c['id']] = {'l': [str(max(lengths))]*2, 'b5': [str(max(widths))]*2,
                            'b1': [b['lower_decimal'], b['upper_decimal']]}
    return refs

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['preflight','timing']);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--gate',type=Path,default=BASE/'sterimol_competitive_preflight/summary.json')
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=[Path(__file__),HERE/'sterimol_reference_worker.py',HERE/'STERIMOL_BENCHMARK_PROTOCOL.md',
             HERE/'host_activity.py',HERE/'benchmark_artifacts.py',HERE/'n1_accuracy_measurements.json',
             BUILD/'sterimol_workflow',BUILD/'cases.json',BUILD/'receipt.json']
    receipt={'phase':a.phase,'bindings':{str(p):digest(p) for p in sources},'runs':[],
             'target_angstrom':'0.01','complete':False}
    def save(): (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS='1')
    env.pop('STERICX_PROFILE_PATH',None)
    if a.phase=='preflight':
        receipt['reference']=independent_references()
        receipt['reference_scope']='80-digit independent coordinate formulas and independent angular search bracket; empirical validation, no universal floating-point certificate'
        version_cmd=[str(BASE.parent/'paton/reference/venv/bin/python'),'-c',
                     'import importlib.metadata as m; print({k:m.version(k) for k in ["morfeus-ml","dbstep","numpy","scipy"]})']
        receipt['packages']=subprocess.check_output(version_cmd,text=True,env=env).strip()
        iterations=1;pairs=[0]
    else:
        gate=json.loads(a.gate.read_text());assert gate['complete']
        for path,sha in gate['bindings'].items():assert digest(Path(path))==sha,path
        receipt['gate_sha256']=digest(a.gate);receipt['admitted']=gate['admitted']
        receipt['host_preflight']=measure();save()
        assert receipt['host_preflight']['quiet'],'No timing: host not quiet'
        iterations=100;pairs=range(-1,8)
    for pair in pairs:
        order=METHODS if pair%2==0 else METHODS[::-1]
        for method in order:
            cmd=['taskset','-c','2',*command(method,iterations)]
            folder=out/f'{pair}_{method}';folder.mkdir()
            (folder/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
            with (folder/'stdout.json').open('w') as stdout,(folder/'stderr.txt').open('w') as stderr:
                result=subprocess.run(cmd,env=env,stdout=stdout,stderr=stderr)
            assert result.returncode==0,(method,result.returncode)
            data=json.loads((folder/'stdout.json').read_text())
            row={'method':method,'pair':pair,'warmup':pair<0,'result':data,'stdout_sha256':digest(folder/'stdout.json')}
            if a.phase=='preflight':
                errors=[]
                with localcontext() as ctx:
                    ctx.prec=80
                    for value in data['rows']:
                        for quantity in ['l','b1','b5']:
                            observed=Decimal.from_float(value[quantity])
                            lo,hi=map(Decimal,receipt['reference'][value['id']][quantity])
                            upper=max(abs(observed-lo),abs(observed-hi))
                            errors.append({'id':value['id'],'quantity':quantity,'absolute_error_upper':str(upper),
                                           'passed':upper<=Decimal('0.01')})
                row['errors']=errors;row['admitted']=all(e['passed'] for e in errors)
            else:
                expected=next(r for r in gate['runs'] if r['method']==method)['result']['rows']
                assert data['rows']==expected,(pair,method)
            receipt['runs'].append(row);save()
    for path,sha in receipt['bindings'].items():assert digest(Path(path))==sha,path
    if a.phase=='preflight':receipt['admitted']={r['method']:r['admitted'] for r in receipt['runs']}
    else:receipt['host_postflight']=postflight(out)
    receipt['complete']=True;save()

if __name__=='__main__':main()
