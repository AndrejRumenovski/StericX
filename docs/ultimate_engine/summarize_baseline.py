"""Recompute baseline tables and analytical errors from retained raw measurements."""
import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from dbstep_resolution_reference import sphere_intersection_interval  # noqa: E402


def distribution(values):
    median = statistics.median(values)
    return {'values':values,'median':median,
            'mad':statistics.median(abs(v-median) for v in values),
            'min':min(values),'max':max(values)}


def main():
    local = ROOT/'.stericx/ultimate_engine'
    output = ROOT/'docs/ultimate_engine'
    cli = json.loads((local/'baseline_cli/journal.json').read_text())
    kernels = json.loads((local/'baseline_kernels/journal.json').read_text())
    summary = {'cli':{}, 'kernels':{}, 'accuracy':[]}
    for threads in [1,2,4,6,12]:
        for name in ['single_small','single_large','conformers_56','ensemble_sdf',
                     'descriptors_10000','parse_ensembles_1000','db_build_ensembles']:
            rows = [r for r in cli['runs'] if r['threads']==threads and
                    r['name'].split('/')[1]==name and not r['instrumented'] and not r['warmup']]
            summary['cli'][f'{name}/t{threads}'] = {
                'wall_ms':distribution([r['metrics']['wall_ns']/1e6 for r in rows]),
                'cpu_ms':distribution([1000*(r['metrics']['cpu_user_s']+r['metrics']['cpu_system_s']) for r in rows]),
                'rss_mib':distribution([r['metrics']['peak_rss_bytes']/2**20 for r in rows]),
                'scientific':rows[0]['scientific']}
    for name in sorted({r['name'].split('/')[0] for r in kernels['runs']} - {'batch'}):
        rows = [r for r in kernels['runs'] if r['name'].split('/')[0]==name]
        summary['kernels'][name] = {
            'us_per_call':distribution([r['kernel']['elapsed_ns']/r['kernel']['iterations']/1000
                        for r in rows if r['kernel']['allocations'] is None])}
        diagnostic = next(r['kernel'] for r in rows if r['kernel']['allocations'] is not None)
        summary['kernels'][name]['allocations_per_call'] = {
            k: v/diagnostic['iterations'] for k,v in diagnostic['allocations'].items()
            if k!='process_peak_live_bytes'}
    for row in json.loads((local/'baseline/accuracy.json').read_text()):
        centers, radii = row['centers'], row['radii']
        if not radii:
            interval = (0.,0.)
        elif len(radii)==1:
            assert centers[0][0]==centers[0][1]==0
            interval = sphere_intersection_interval(row['radius'], radii[0], abs(centers[0][2]))
        else:
            assert row['id']=='overlapping_pair' and radii==[1.25,1.25]
            one = sphere_intersection_interval(3.5,1.25,0.5)
            common = sphere_intersection_interval(1.25,1.25,1.)
            interval = (math.nextafter(2*one[0]-common[1],-math.inf),
                        math.nextafter(2*one[1]-common[0],math.inf))
        reference = sum(interval)/2
        observed = row['result']['summary']['buried_volume']
        summary['accuracy'].append({**{k:row[k] for k in ['id','radius','density']},
            'volume':observed,'analytic_interval':interval,'reference':reference,
            'absolute_error':abs(observed-reference),
            'relative_error':abs(observed-reference)/abs(reference) if reference else None,
            'ms_per_call':distribution([x/1e6 for x in row['times_ns']])})
    (output/'baseline_measurements.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines = ['# Baseline profiling tables', '',
        'Exploratory measurements; three warmed CLI observations per configuration, all retained.', '',
        '| Workload | Threads | Median ms | MAD ms | Range ms | Median RSS MiB |',
        '|---|---:|---:|---:|---:|---:|']
    for key,row in summary['cli'].items():
        name,threads = key.split('/t');w=row['wall_ms']
        lines.append(f"| {name} | {threads} | {w['median']:.3f} | {w['mad']:.3f} | {w['min']:.3f}–{w['max']:.3f} | {row['rss_mib']['median']:.3f} |")
    lines += ['', '| Kernel | Median µs | Allocations/call | Reallocations/call | Requested bytes/call |',
              '|---|---:|---:|---:|---:|']
    for key,row in summary['kernels'].items():
        a=row['allocations_per_call']
        lines.append(f"| {key} | {row['us_per_call']['median']:.3f} | {a['allocation_calls']:g} | {a['reallocation_calls']:g} | {a['requested_bytes']:g} |")
    lines += ['', '| Geometry | Density Å³ | Native volume Å³ | Analytic volume Å³ | Absolute error Å³ | Median ms |',
              '|---|---:|---:|---:|---:|---:|']
    for row in summary['accuracy']:
        lines.append(f"| {row['id']} | {row['density']:.9g} | {row['volume']:.10g} | {row['reference']:.10g} | {row['absolute_error']:.8g} | {row['ms_per_call']['median']:.4f} |")
    (output/'BASELINE_TABLES.md').write_text('\n'.join(lines)+'\n')
    print('wrote',output/'baseline_measurements.json')


if __name__ == '__main__':
    main()
