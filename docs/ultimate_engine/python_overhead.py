"""Separate Python conversion, array export and ordered batch overhead."""
import argparse
import hashlib
import json
import os
import statistics
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import stericx as sx
from stericx._types import ensemble_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    os.sched_setaffinity(0, {2})
    ensemble = sx.Ensemble.from_file(root/'python/examples/data/ligand_401_ensemble.xyz')
    options = sx.DescriptorOptions(axis='coordination', donor_index=14)
    encoded = json.dumps(asdict(options))
    raw = ensemble._engine.analyze(encoded, '{"missing":"require_complete"}', 'null', None)
    parsed = json.loads(raw)
    result = ensemble_result(parsed)
    assert result.to_dict() == parsed
    array = result.to_numpy()
    assert len(array) == 6
    for i, row in enumerate(result.to_records()):
        for field in array.dtype.names:
            if array[field].dtype.kind == 'f' and row[field] is not None:
                assert array[field][i] == row[field]
    manifest = root/'.stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2/workloads/manifest.json'
    paths = [p for p in json.loads(manifest.read_text())['workloads']['conformers_56']['argv'] if p.endswith('.xyz')]
    assert len(paths) == 56
    batch_options = sx.DescriptorOptions(axis='coordination')
    def individual():
        return tuple(sx.Molecule.from_file(p).descriptors(batch_options) for p in paths)
    expected = individual()
    for threads in [1, 2, 4, 6, 12]:
        rows = sx.calculate_files(paths, batch_options, threads=threads)
        assert tuple(r.conformers[0] for r in rows) == expected
    lanes = {
        'json_decode': (lambda: json.loads(raw), 500),
        'typed_construction': (lambda: ensemble_result(parsed), 500),
        'numpy_export': (result.to_numpy, 500),
        'options_encode': (lambda: sx._options(options), 1000),
        'parsed_input_record_export': (lambda: ensemble.records, 200),
        'per_file_56': (individual, 3),
    }
    for threads in [1, 2, 4, 6, 12]:
        lanes[f'batch_56_t{threads}'] = (lambda t=threads: sx.calculate_files(paths, batch_options, threads=t), 3)
    journal = {'module':sx.__file__, 'native':sx._native.__file__,
        'native_sha256':hashlib.sha256(Path(sx._native.__file__).read_bytes()).hexdigest(),
        'numpy_version':np.__version__, 'scientific_equivalence':True,
        'cache_policy':'Warmed imports and parsed conversion inputs. File calls include parsing.', 'runs':{}}
    for name, (call, iterations) in lanes.items():
        threads = int(name.rsplit('t',1)[1]) if name.startswith('batch_') else 1
        affinity = {2} if threads==1 else {2,3} if threads==2 else set(range(threads))
        os.sched_setaffinity(0, affinity)
        call()
        values = []
        for _ in range(5):
            start = time.perf_counter_ns()
            for _ in range(iterations):
                value = call()
                del value
            values.append((time.perf_counter_ns()-start)/iterations)
        median = statistics.median(values)
        journal['runs'][name] = {'iterations':iterations,'affinity':sorted(affinity),
            'ns_per_call':values,'median_ns':median,
            'mad_ns':statistics.median(abs(v-median) for v in values)}
        (args.output/'summary.json').write_text(json.dumps(journal,indent=2)+'\n')
    print(json.dumps({k:v['median_ns'] for k,v in journal['runs'].items()}))


if __name__ == '__main__':
    main()
