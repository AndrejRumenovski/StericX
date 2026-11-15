#!/usr/bin/env python3
"""Quiet-window costs with batched analytical calls and unchanged scientific references."""
import argparse,hashlib,json,math,statistics,subprocess
from pathlib import Path
from host_activity import measure
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'.stericx/ultimate_engine'
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
receipt={'scope':'Actual Rust algorithms; setup/parsing outside internal timers; analytical and general-union original geometries selected prospectively, all prior stress evidence retained',
    'host_preflight':measure(),'complete':False,'runs':{}}
def save():(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
save()
if not receipt['host_preflight']['quiet']:raise SystemExit('No timing: quiet window required.')
old={r['id']:r for r in json.loads((ROOT/'docs/ultimate_engine/n2_accuracy_measurements.json').read_text())['rows']}
selected={'empty','centered','full','lens','internal_tangency','external_tangency','near_external_tangency',
          'overlapping_pair','three_disjoint','three_identical','three_overlapping_chain'}
cases=json.loads((BASE/'n2_validation_v4/requests.json').read_text());cases=[c for c in cases if c['id'].split('/')[0] in selected]
for case in cases:
    case['samples']=8;case['iterations']=max(1,min(100000,math.ceil(20_000_000/max(old[case['id']]['median_ns'],1))))
(out/'volume_requests.json').write_text(json.dumps(cases,indent=2)+'\n')
for name,binary,request in [('sterimol',BASE/'final_profile/bin/continuous',BASE/'n1_validation_expanded/cases.json'),
                           ('volume',BASE/'final_profile/bin/continuous_volume',out/'volume_requests.json')]:
    cmd=['taskset','-c','2',str(binary),str(request)]
    with (out/(name+'.jsonl')).open('w') as stream,(out/(name+'.stderr')).open('w') as err:r=subprocess.run(cmd,stdout=stream,stderr=err)
    r.check_returncode();rows=list(map(json.loads,(out/(name+'.jsonl')).read_text().splitlines()))
    receipt['runs'][name]={'command':cmd,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
                          'input_sha256':hashlib.sha256(request.read_bytes()).hexdigest(),'rows':rows}
    if name=='volume':
        for row in rows:
            reference=old[row['id']];assert row['result']==reference['result'],row['id']
            row['ns_per_call']=[v/row['iterations'] for v in row['elapsed_ns']]
            row['median_ns']=statistics.median(row['ns_per_call'])
            for key in ['reference_interval','absolute_error_lower','absolute_error_upper']:row[key]=reference[key]
    else:
        previous={r['id']:r for r in map(json.loads,(BASE/'n1_validation_expanded/native.jsonl').read_text().splitlines())}
        for row in rows:
            for current,prior in zip(row['measurements'],previous[row['id']]['measurements'],strict=True):
                assert current['result']==prior['result'],row['id']
    save()
receipt.update(complete=True,host_postflight=measure());save()
