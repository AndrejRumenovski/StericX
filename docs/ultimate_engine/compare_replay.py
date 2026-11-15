#!/usr/bin/env python3
"""Compare every frozen public/private response without numerical tolerances."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifests = [json.loads((p/'manifest.json').read_text()) for p in
                 (args.baseline, args.candidate)]
    assert all(m['complete_capture'] for m in manifests)
    assert manifests[0]['lanes'].keys() == manifests[1]['lanes'].keys()
    result = {'scope':'Exact parsed JSON including serialized float bits and errors',
              'manifests': {str(p/'manifest.json'):digest(p/'manifest.json') for p in
                            (args.baseline,args.candidate)}, 'lanes':{}, 'differences':[]}
    for lane in manifests[0]['lanes']:
        left, right = [p/lane/'stdout.jsonl' for p in (args.baseline,args.candidate)]
        count = 0
        with left.open() as a, right.open() as b:
            for count, (x,y) in enumerate(itertools.zip_longest(a,b),1):
                if x is None or y is None or json.loads(x) != json.loads(y):
                    result['differences'].append({'lane':lane,'row':count,
                        'baseline':None if x is None else json.loads(x),
                        'candidate':None if y is None else json.loads(y)})
        result['lanes'][lane] = {'rows':count,'baseline_sha256':digest(left),
                                'candidate_sha256':digest(right)}
    result['total_rows'] = sum(row['rows'] for row in result['lanes'].values())
    result['exact'] = not result['differences']
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'rows':result['total_rows'],'differences':len(result['differences'])}))
    raise SystemExit(0 if result['exact'] else 1)

if __name__ == '__main__':
    main()
