#!/usr/bin/env python3
"""Cross-check molecular continuous-volume enclosures against the prior oracle."""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from dbstep_resolution_reference import certified_union_volume
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--engine',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output;out.mkdir(exist_ok=False)
source=ROOT/'.stericx/ultimate_engine/baseline_observations/geometry/stdout.jsonl'
fixtures=list(map(json.loads,source.read_text().splitlines()))[:12]
requests=[]
for fixture in fixtures:
    atoms=[{'center':atom['position'],'radius':atom['vdw_radius']*1.1699999570846558}
           for atom in fixture['atoms'] if atom['element'].upper()!='H']
    for method in [{'method':'reference','absolute_tolerance':.05,'max_boxes':262143},
                   {'method':'midpoint_slices','subdivisions':128},
                   {'method':'midpoint_slices','subdivisions':256}]:
        requests.append({'id':fixture['id']+'/'+method['method']+str(method.get('subdivisions','')),
                         'region':{'center':fixture['center'],'radius':3.5},'atoms':atoms,
                         'method':method,'samples':1})
(out/'requests.json').write_text(json.dumps(requests,indent=2)+'\n')
with (out/'native.jsonl').open('w') as stdout,(out/'native.stderr').open('w') as stderr:
    subprocess.run([str(a.engine.resolve()),str((out/'requests.json').resolve())],stdout=stdout,stderr=stderr,check=True)
native=list(map(json.loads,(out/'native.jsonl').read_text().splitlines()))
receipt={'engine_sha256':hashlib.sha256(a.engine.read_bytes()).hexdigest(),
         'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'controller_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'reference_sha256':hashlib.sha256((ROOT/'scripts/dbstep_resolution_reference.py').read_bytes()).hexdigest(),
         'scope':'Independent octree/rational-volume oracle versus binary heap/rounded-volume enclosure; overlap is corroboration, not an exact molecular-volume certificate.',
         'rows':[],'failures':[]}
for i in range(0,len(requests),3):
    q=requests[i];start=time.monotonic()
    reference=certified_union_volume([x['center'] for x in q['atoms']],[x['radius'] for x in q['atoms']],3.5,
                                     center=q['region']['center'],max_depth=6)
    result=native[i]['result']
    intersection=[max(reference['lower'],result['lower_bound_angstrom3']),min(reference['upper'],result['upper_bound_angstrom3'])]
    if intersection[0]>intersection[1]:receipt['failures'].append(q['id'])
    row={'id':q['id'],'python_reference':reference,'python_seconds':time.monotonic()-start,
         'common_enclosure':intersection,'native':native[i:i+3]}
    for observation in row['native']:
        v=observation['result']['volume_angstrom3']
        observation['absolute_error_enclosure']=[max(0.,intersection[0]-v,v-intersection[1]),max(abs(intersection[0]-v),abs(intersection[1]-v))]
    receipt['rows'].append(row)
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(q['id'],intersection,result['target_met'],flush=True)
receipt['complete']=True
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
