#!/usr/bin/env python3
"""Exact parser replay and interleaved, exploratory U1 performance experiment."""
import gzip
import hashlib
import json
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '.stericx/ultimate_engine/u1_validation'
OUT.mkdir(exist_ok=False)
BASE = ROOT / '.stericx/ultimate_engine/baseline/bin'
CAND = ROOT / '.stericx/ultimate_engine/u1/bin'
REFERENCE = ROOT / '.stericx/paton/reference/DBSTEP/tests/cube_files'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = {'binaries': {str(p): sha(p) for b in (BASE, CAND) for p in b.iterdir()},
           'scope': 'exploratory pinned CPU 2; final quiet-window A/B remains separate',
           'inputs': {}, 'equivalence': [], 'timings': [], 'allocations': []}
fixtures = [REFERENCE / n for n in ('Ne_medium.cube', 'benzene_coarse.cube', 'CH2CMe3_100.cube')]
header = 'title\ntotal density\n1 -1 -2 -3 1\n2 1 0 0\n2 0 2 0\n2 0 0 3\n6 6 0.5 1 2\n'
for i, token in enumerate(('0', '-0', '1e-3', '1E-3', '1d-3', '1D-3', '+1.5',
                            '1e-310', '1e-999', 'NaN', 'nan', 'inf', '-inf',
                            '1e999', '-0.01', '1d', 'D', '1d-3d2', 'abc', '1_000')):
    p = OUT / f'token_{i}.cube'
    p.write_text(header + f'0 {token} 2 3 4 5 6 7\n')
    fixtures.append(p)
for i, p in enumerate(fixtures):
    receipt['inputs'][str(p)] = sha(p)
    responses = []
    for label, b in [('baseline', BASE), ('candidate', CAND)]:
        run = subprocess.run([str(b/'dump_density'), str(p)], capture_output=True, check=True)
        with gzip.open(OUT/f'{i}_{label}.json.gz', 'wb') as f:
            f.write(run.stdout)
        responses.append(run.stdout)
    equal = responses[0] == responses[1]
    receipt['equivalence'].append({'input': str(p), 'exact': equal,
                                   'response_sha256': hashlib.sha256(responses[0]).hexdigest()})
    assert equal, p
for pair in range(8):
    order = [('baseline', BASE), ('candidate', CAND)]
    if pair % 2:
        order.reverse()
    for label, b in order:
        command = ['taskset', '-c', '2', str(b/'observer'), 'density_parse',
                   str(fixtures[0]), 'unused', '40']
        run = subprocess.run(command, capture_output=True, check=True)
        result = json.loads(run.stdout)
        receipt['timings'].append({'pair': pair, 'variant': label, 'command': command, **result})
for label, b in [('baseline', BASE), ('candidate', CAND)]:
    command = [str(b/'observer_alloc'), 'density_parse', str(fixtures[0]), 'unused', '10']
    result = json.loads(subprocess.check_output(command))
    receipt['allocations'].append({'variant': label, 'command': command, **result})
medians = {label: statistics.median(r['elapsed_ns']/r['iterations'] for r in receipt['timings']
                                  if r['variant'] == label) for label in ('baseline','candidate')}
receipt['median_ns'] = medians
receipt['candidate_speedup'] = medians['baseline']/medians['candidate']
receipt['complete'] = True
(OUT/'summary.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({'medians_ns': medians, 'speedup': receipt['candidate_speedup'],
                  'exact_cases': len(fixtures)}, indent=2))
