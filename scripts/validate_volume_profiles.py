#!/usr/bin/env python3
"""Validate each radius profile position against frozen and fresh independent calls."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from dbstep import Dbstep
from study_conformer_ensembles import sha,f32
from validate_dbstep_volume import load_cases,grid_key,matched_case,independent_lattice
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    all_cases=load_cases();cases={c['id']:c for c in all_cases}
    base=cases['SIG-NIHDA-401/0']
    path=(a.output/'401_0.xyz').resolve();path.write_text(base['row']['source_geometry'])
    radii=[3.5,-1.,2.5,4.5,3.5]
    requests=[{'mode':'hard_sphere','path':str(path),'options':base['options'],'radii_angstrom':radii}]
    density_root=ROOT/'docs/validation/dbstep/density_v2'
    density_requests=[json.loads(s) for s in (density_root/'requests.jsonl').read_text().splitlines()]
    for i in [0,3,6]:
        r=copy.deepcopy(density_requests[i]);r.pop('radius_angstrom');r.update(mode='density',radii_angstrom=[1.,1.5,2.,10.,1.])
        requests.append(r)
    raw=''.join(json.dumps(r)+'\n' for r in requests);(a.output/'requests.jsonl').write_text(raw)
    binary=ROOT/'target/debug/examples/volume_profile'
    proc=subprocess.run([str(binary)],input=raw,text=True,capture_output=True)
    (a.output/'native_stdout.jsonl').write_text(proc.stdout);(a.output/'native_stderr.txt').write_text(proc.stderr)
    observed=[json.loads(s) for s in proc.stdout.splitlines()]
    if proc.returncode or len(observed)!=len(requests):raise RuntimeError('incomplete profile output')
    failures=[];reference=[]
    grid_map={};captured={}
    for row,c in zip((ROOT/'docs/validation/dbstep/volume_v2/native_stdout.jsonl').read_text().splitlines(),all_cases):
        r=json.loads(row);captured[c['id']]=r
        if r['dump']['points'] is not None:grid_map[grid_key(c['options']['volume'])]=f32(r['dump']['points'])
    profile=observed[0]['value'][0]
    for sample in profile['samples']:
        radius=sample['requested_radius_angstrom']
        if radius<=0:
            if sample['outcome']['status']!='error':failures.append({'radius':radius,'error':'invalid radius was accepted'})
            continue
        identity='SIG-NIHDA-401/0'+(f'/sphere_radius={radius}' if radius!=3.5 else '')
        case=cases[identity];capture=copy.deepcopy(captured[identity]);value=sample['outcome']['value']
        for field,v in value['summary'].items():
            if np.float32(v).tobytes()!=np.float32(capture['buried_volume'][field]).tobytes():failures.append({'radius':radius,'field':field,'error':'profile differs from independent single calculation'})
        for frame,old in zip(value['orientations'],capture['dump']['orientations']):
            if frame['plane_atom_index']!=old['plane'] or not np.array_equal(f32(frame['axes']),f32(old['basis'])):
                failures.append({'radius':radius,'error':'frame mapping differs'})
            for field in ['buried_volume','near_vbur','far_vbur','quadrants','octants']:
                if not np.array_equal(f32(frame[field]),f32(old[field])):failures.append({'radius':radius,'field':field,'error':'frame values differ'})
        frames,summaries,errors,_=matched_case(case,capture,grid_map[grid_key(case['options']['volume'])])
        api,details=independent_lattice(case,capture,.05)
        known=radius==2.5
        if (errors or any(not s['passed'] for s in summaries)) and not known:failures.append({'radius':radius,'error':'unexpected independent comparison failure'})
        reference.append({'mode':'hard_sphere','radius':radius,'matched_frames':frames,'matched_summaries':summaries,
                          'strict_errors':errors,'known_boundary_difference':known,'independent_lattice':api,'reference_settings':details})
    for request,response in zip(requests[1:],observed[1:]):
        profile=response['value']
        for sample in profile['samples']:
            radius=sample['requested_radius_angstrom']
            if radius==10.:
                if sample['outcome']['status']!='error':failures.append({'error':'out-of-domain density radius was accepted'})
                continue
            value=sample['outcome']['value']
            single=next(json.loads(line)['value'] for line,req in zip((density_root/'native_stdout.jsonl').read_text().splitlines(),density_requests)
                        if req['path']==request['path'] and req['radius_angstrom']==radius)
            if value!=single:failures.append({'path':request['path'],'radius':radius,'error':'density profile changed single-radius result'})
            public=Dbstep.dbstep(request['path'],surface='density',volume=True,atom1=1,r=radius,isoval=request['threshold_electrons_per_bohr3'],quiet=True,norot=True)
            residual=value['percent_buried_volume']-float(public.bur_vol)
            limit=64*float(np.finfo(float).eps)*max(1.,abs(float(public.bur_vol)))
            if abs(residual)>limit:failures.append({'path':request['path'],'radius':radius,'residual':residual})
            reference.append({'mode':'density','path':request['path'],'radius':radius,'native_percent':value['percent_buried_volume'],
                              'dbstep_percent':float(public.bur_vol),'residual':residual,'arithmetic_limit':limit})
    summary={'profile_requests':len(requests),'positions':sum(len(r['radii_angstrom']) for r in requests),
             'reference_comparisons':len(reference),'expected_error_positions':4,'failures':failures,
             'profile_preservation_pass':not failures,'hard_sphere_strict_gate_all_positions':False,
             'limitation':'R=2.5 hard-sphere boundary arithmetic difference remains; see DISCREPANCIES.md'}
    (a.output/'reference_raw.json').write_text(json.dumps(reference,indent=2)+'\n');(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (a.output/'receipt.json').write_text(json.dumps({'command':sys.argv,'script_sha256':sha(__file__),'binary_sha256':sha(binary),
        'sources':{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/profiles.rs',ROOT/'examples/volume_profile.rs']},
        'inputs':{str(root/'manifest.json'):sha(root/'manifest.json') for root in [density_root,ROOT/'docs/validation/dbstep/volume_v2']}},indent=2)+'\n')
    (a.output/'manifest.json').write_text(json.dumps({p.name:sha(p) for p in sorted(a.output.iterdir()) if p.is_file()},indent=2)+'\n')
    print(json.dumps(summary,indent=2));return bool(failures)
if __name__=='__main__':raise SystemExit(main())
