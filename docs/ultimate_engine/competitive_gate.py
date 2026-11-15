#!/usr/bin/env python3
"""Fresh scientific admission using unchanged exact-decimal/IEEE references."""
import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
H=ROOT/'docs/benchmarks/morfeus_current'
sys.path.insert(0,str(H/'revision_policy'))
import reference as independent
spec=importlib.util.spec_from_file_location('historical_evaluator',H/'revision_policy/evaluate.py')
evaluator=importlib.util.module_from_spec(spec);spec.loader.exec_module(evaluator)
sys.path.insert(0,str(Path(__file__).resolve().parent))
import competitive_worker as worker

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=a.build.resolve()/'adapter/target/release'
    manifest=json.loads((H/'input_manifest.json').read_text())
    paths=[H/r['filename'] for r in manifest['corpus']]
    for path,row in zip(paths,manifest['corpus'],strict=True):assert sha(path)==row['sha256']
    bound=[Path(__file__),Path(worker.__file__),Path(independent.__file__),
        H/'revision_policy/evaluate.py',H/'EQUIVALENCE_POLICY.md',H/'input_manifest.json',
        worker.HISTORICAL,binary/'volume-batch',binary/'policy-grid-observer',
        a.build/'receipt.json',a.build/'adapter/policy_copied_volume.rs']
    receipt={'scope':'Exact finite-grid estimator; DBSTEP occupancy-kernel adapter',
        'bindings':{str(p):sha(p) for p in bound},'rows':[],'failures':[],'complete':False}
    def save():(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    save()
    env=dict(os.environ,RAYON_NUM_THREADS='1')
    def run(name,command):
        with (out/(name+'.stdout')).open('w') as stream,(out/(name+'.stderr')).open('w') as err:
            r=subprocess.run(list(map(str,command)),env=env,stdout=stream,stderr=err)
        assert r.returncode==0,(name,r.returncode)
    run('native',[binary/'volume-batch',*paths])
    native=json.loads((out/'native.stdout').read_text());assert len(native)==56
    run('probe',[binary/'policy-grid-observer',*paths])
    with (out/'probe.stdout').open() as probes:
        for i,(path,native_public,line) in enumerate(zip(paths,native,probes,strict=True)):
            n=json.loads(line);ref=independent.certify(path);m=evaluator.reference_tool(path)
            dbframes,dbvalues=worker.dbstep(str(path),details=True)
            nf,nv,nvolume=independent.reduce_counts(ref['frames'],24)
            mf,mv,mvolume=independent.reduce_counts(ref['frames'],53)
            issues=[];arrays={}
            def check(condition,kind):
                if not condition:issues.append(kind)
            for key in ['donor','neighbors']:check(n[key]==ref[key]==m[key],key)
            check(n['sphere_volume']==nvolume,'native_sphere_volume')
            for tool,observed,expected in [('native',native_public,nv),('morfeus',m['values'],mv),('dbstep',dbvalues,mv)]:
                for key,value in expected.items():check(observed[key]==value,f'{tool}:{key}')
            idx=np.array([r['index'] for r in n['points']]);ids=idx[:,0]*1024+idx[:,1]*32+idx[:,2]
            order=np.argsort(ids);check(np.array_equal(ids[order],independent.IDS),'native_grid')
            check(np.array_equal(np.array([r['region'] for r in n['points']])[order],independent.REGIONS),'native_regions')
            for j,(rf,npf,mpf,dpf,nmodel,mmodel) in enumerate(zip(ref['frames'],n['frames'],m['frames'],dbframes,nf,mf,strict=True)):
                check(rf['plane']==npf['plane']==mpf['plane']==dpf['plane'],f'plane{j}')
                check(not np.any(rf['unknown']),f'uncertain{j}')
                def ids_of(points):
                    ijk=np.rint((points+3.5)*31/7).astype(int)
                    return ijk[:,0]*1024+ijk[:,1]*32+ijk[:,2]
                check(np.array_equal(np.sort(mpf['grid']),independent.IDS),f'morfeus_grid{j}')
                check(np.array_equal(np.sort(ids_of(dpf['points'])),independent.IDS),f'dbstep_grid{j}')
                masks={'native':np.array(npf['occupied'])[order],
                    'morfeus':np.isin(independent.IDS,mpf['buried']),
                    'dbstep':np.isin(independent.IDS,ids_of(dpf['occupied']))}
                for tool,mask in masks.items():
                    check(np.array_equal(mask,rf['occupied']),f'{tool}:occupancy{j}')
                    arrays[f'{tool}{j}']=mask
                arrays[f'reference{j}']=rf['occupied'];arrays[f'unknown{j}']=rf['unknown']
                for tool,observed,model in [('native',npf,nmodel),('morfeus',mpf,mmodel),('dbstep',dpf,mmodel)]:
                    for key,value in model.items():check(observed[key]==value,f'{tool}:{key}:{j}')
            np.savez_compressed(out/f'masks_{i:02}.npz',**arrays)
            row={'path':str(path),'issues':issues,'native':native_public,'morfeus':m['values'],'dbstep':dbvalues}
            receipt['rows'].append(row)
            receipt['failures'].extend({'path':str(path),'issue':issue} for issue in issues)
            save();print(i+1,len(issues),flush=True)
    for path,digest in receipt['bindings'].items():assert sha(Path(path))==digest,path
    receipt.update(complete=True,passed=not receipt['failures'],structures=56,frames=168,
                   certified_point_decisions=56*3*len(independent.IDS));save()
    raise SystemExit(0 if receipt['passed'] else 1)

if __name__=='__main__':main()
