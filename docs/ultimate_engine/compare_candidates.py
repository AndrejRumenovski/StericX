#!/usr/bin/env python3
"""Alternating source-bound kernel A/B experiment, with exact result checks."""
import argparse
import hashlib
import json
import os
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import profile_stericx as harness

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', required=True)
p.add_argument('--candidate', required=True)
p.add_argument('--output', required=True)
p.add_argument('--modes', nargs='+', default=['volume', 'regions', 'ensemble', 'radius_profile', 'sterimol'])
p.add_argument('--pairs', type=int, default=8)
a = p.parse_args()
root = ROOT/'.stericx/ultimate_engine'
out = root/a.output
out.mkdir(exist_ok=False)
launcher = harness.build_launcher(out/'launcher')
inputs = [ROOT/'python/examples/data/ligand_401_ensemble.xyz',
          ROOT/'python/examples/data/ligand_401_options.json']
env = dict(os.environ, LC_ALL='C', RAYON_NUM_THREADS='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
env.pop('STERICX_PROFILE_PATH', None)
binaries = {label: root/name/'bin/observer' for label, name in [('baseline',a.baseline),('candidate',a.candidate)]}
receipt = {'argv': sys.argv, 'scope': 'exploratory desktop, CPU 2', 'cpu_before': Path('/proc/stat').read_text(),
           'hashes': {str(p): harness.digest(p) for p in [*inputs,*binaries.values()]}, 'runs': [], 'summary': {}}
counts = {'volume':150,'regions':150,'ensemble':30,'radius_profile':30,'sterimol':10000}
expected = {}
for mode in a.modes:
    for pair in range(a.pairs):
        labels = ['baseline','candidate'] if pair%2 == 0 else ['candidate','baseline']
        for label in labels:
            directory = out/f'{mode}_{pair}_{label}'
            command = [str(binaries[label]), mode, *map(str,inputs), str(counts[mode])]
            metrics = harness.measure(command, env, {2}, directory, launcher)
            assert metrics.get('exit_code',metrics.get('returncode',0)) == 0, metrics
            result = json.loads((directory/'stdout.txt').read_text())
            key = json.dumps(result['result'],sort_keys=True)
            assert key == expected.setdefault(mode,key), (mode,label,'scientific mismatch')
            receipt['runs'].append({'mode':mode,'pair':pair,'variant':label,'metrics':metrics,
                                     'elapsed_ns':result['elapsed_ns'],'iterations':result['iterations'],
                                     'scientific_sha256':hashlib.sha256(key.encode()).hexdigest()})
            (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    values = {label:[r['elapsed_ns']/r['iterations'] for r in receipt['runs']
                     if r['mode']==mode and r['variant']==label] for label in binaries}
    receipt['summary'][mode] = {label:{'values_ns':v,'median_ns':statistics.median(v),
                                      'mad_ns':statistics.median(abs(x-statistics.median(v)) for x in v)}
                               for label,v in values.items()}
    receipt['summary'][mode]['speedup'] = statistics.median(values['baseline'])/statistics.median(values['candidate'])
    print(mode, receipt['summary'][mode]['speedup'], flush=True)
receipt['cpu_after'] = Path('/proc/stat').read_text()
receipt['complete'] = True
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
