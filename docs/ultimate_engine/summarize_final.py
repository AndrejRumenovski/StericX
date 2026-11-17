#!/usr/bin/env python3
"""Derive final tables and complete distributions from immutable measured records."""
import copy
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import median, mean

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
BASE=ROOT/'.stericx/ultimate_engine'
WINDOW=BASE/'final_window_6'
inputs={}
def read(p):
    p=Path(p);inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text())
def write(name,d): (HERE/name).write_text(json.dumps(d,indent=2)+'\n')
def stats(values):
    values=list(values);m=median(values)
    return {'values':values,'median':m,'mad':median(abs(v-m) for v in values),
            'min':min(values),'max':max(values),'n':len(values)}
def metric_stats(rows):
    return {key:stats(r['metrics'][key] for r in rows) for key in
            ['wall_ns','cpu_user_s','cpu_system_s','peak_rss_bytes','minor_faults','major_faults',
             'voluntary_context_switches','involuntary_context_switches','cpu_percent_one_core']}
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|',
                     *['| '+' | '.join(map(str,r))+' |' for r in rows]])+'\n'
def quiet(d):
    if isinstance(d,dict):
        if 'quiet' in d:assert d['quiet'] is True
        for v in d.values():quiet(v)
    elif isinstance(d,list):
        for v in d:quiet(v)

def main():
    result={'scope':'Final quiet-window native measurements; warmups retained in raw sources, no trial rejection. All samples in each statistic retained.',
            'statistic':'Median, unscaled median absolute deviation, minimum/maximum; independent processes are trial units. Ratios of medians, not confidence intervals.',
            'hardware_counters':'Unavailable; simulated Callgrind/Cachegrind events are separately labelled.',
            'inputs':inputs}
    receipt=read(WINDOW/'receipt.json');assert receipt['complete'];quiet(receipt)
    assert len(receipt['runs'])==8 and all(r['status']=='completed' and r['returncode']==0 for r in receipt['runs'])
    result['host_gates']={r['name']:{k:r[k] for k in ['host_preflight','host_postflight']} for r in receipt['runs']}
    cli=read(WINDOW/'cli_ab/summary.json');assert cli['complete'];quiet(cli)
    assert len(cli['attempts'])==630 and all(r['status']=='passed' for r in cli['attempts'])
    result['cli']={}
    for key,s in cli['results'].items():
        workload,t=key.split('/t');threads=int(t);row=copy.deepcopy(s)
        for side in ['baseline','candidate']:
            rr=[r for r in cli['attempts'] if r['workload']==workload and r['threads']==threads and r['side']==side and not r['warmup']]
            assert len(rr)==8
            row[side]['process_metrics']=metric_stats(rr)
            assert median(r['metrics']['wall_ns'] for r in rr)==s[side]['median_ns']
        row['throughput_change_percent']=(s['speedup']-1)*100
        row['elapsed_time_reduction_percent']=(1-1/s['speedup'])*100
        result['cli'][key]=row
    result['scaling']={}
    for t in [1,2,4,6,12]:
        row={}
        for side in ['baseline','candidate']:
            duration=cli['results'][f'descriptors_10000/t{t}'][side]['median_ns']
            speed=cli['results']['descriptors_10000/t1'][side]['median_ns']/duration
            row[side]={'files_per_second':1e13/duration,'speedup_vs_one_thread':speed,'parallel_efficiency':speed/t}
        result['scaling'][str(t)]=row
    kab=read(WINDOW/'kernel_ab/summary.json');assert kab['complete'];result['kernel_ab']=kab['summary']
    py=read(WINDOW/'python_ab/summary.json');assert py['complete'];quiet(py)
    assert len(py['runs'])==208 and all(r['exact'] for r in py['runs'])
    grouped=defaultdict(list)
    for r in py['runs']:grouped[(r['mode'],r['threads'])].append(r)
    result['python_ab']={}
    for (mode,t),rows in grouped.items():
        d={}
        for side in ['baseline','candidate']:
            rr=[r for r in rows if r['side']==side];assert len(rr)==8
            d[side]=stats(r['kernel']['elapsed_ns']/r['kernel']['iterations'] for r in rr)
        d['speedup']=d['baseline']['median']/d['candidate']['median']
        result['python_ab'][f'{mode}/t{t}']=d
    layers=read(WINDOW/'python_native_bridge/summary.json');assert layers['complete'] and layers['scientific_results_exact']
    result['python_layers']=layers['results'];result['python_layer_interpretation']=layers['interpretation']
    kp=read(WINDOW/'kernel_profile/journal.json');assert kp['complete'];assert len(kp['runs'])==71
    grouped=defaultdict(list)
    for r in kp['runs']:grouped[r['name'].rsplit('/',1)[0]].append(r)
    result['kernels']={}
    for mode,rows in grouped.items():
        rr=[r for r in rows if r['kernel']['allocations'] is None];assert len(rr)==3
        a=[r for r in rows if r['kernel']['allocations'] is not None]
        result['kernels'][mode]={'ns_per_call':stats(r['kernel']['elapsed_ns']/r['kernel']['iterations'] for r in rr),
                                'process_metrics':metric_stats(rr),'allocation_run':a}
    cp=read(WINDOW/'cli_profile/journal.json');assert cp['complete'];assert len(cp['runs'])==175
    profile=read(WINDOW/'cli_profile/t1/descriptors_10000/diagnostic/profile.json')
    result['profile']=profile
    native=[r for r in cp['runs'] if r['name'].startswith('t1/descriptors_10000/native_') and not r['warmup']]
    assert len(native)==3
    result['profile_native_distribution']=stats(r['metrics']['wall_ns'] for r in native)
    result['profile_overhead_ratio']=profile['session_wall_ns']/result['profile_native_distribution']['median']
    comp=read(WINDOW/'competitive/summary.json');assert comp['complete'];quiet(comp)
    assert len(comp['runs'])==432 and all(r['status']=='passed' for r in comp['runs'])
    grouped=defaultdict(list)
    for r in comp['runs']:
        if not r['warmup']:grouped[(r['size'],r['threads'],r['scope'],r['tool'])].append(r)
    result['competitive']={}
    for (n,t,scope,tool),rows in grouped.items():
        assert len(rows)==8
        values=[r['metrics']['wall_ns'] if scope=='fresh_process' else mean(r['warm_elapsed_ns']) for r in rows]
        result['competitive'][f'n{n}/t{t}/{scope}/{tool}']={'time_ns':stats(values),'process_metrics':metric_stats(rows),
             'warm_internal_samples':[r['warm_elapsed_ns'] for r in rows]}
    result['competitive_warmups']=[r for r in comp['runs'] if r['warmup']]
    result['competitive_scaling']={}
    for scope in ['fresh_process','warm_workflow']:
        for tool in ['native','morfeus','dbstep']:
            one=result['competitive'][f'n56/t1/{scope}/{tool}']['time_ns']['median']
            result['competitive_scaling'][f'{scope}/{tool}']={str(t):{'speedup':one/result['competitive'][f'n56/t{t}/{scope}/{tool}']['time_ns']['median'],
                 'efficiency':one/result['competitive'][f'n56/t{t}/{scope}/{tool}']['time_ns']['median']/t} for t in [1,2,4,6,12]}
    gate=read(BASE/'sterimol_competitive_preflight/summary.json')
    st=read(BASE/'sterimol_competitive_timing/summary.json');assert st['complete'];quiet(st)
    assert len(st['runs'])==36
    result['sterimol_competitive']={}
    for method in ['sampled361','continuous','morfeus','dbstep']:
        rr=[r for r in st['runs'] if r['method']==method and not r['warmup']];assert len(rr)==8
        g=next(r for r in gate['runs'] if r['method']==method)
        assert all(r['result']['rows']==g['result']['rows'] for r in rr)
        result['sterimol_competitive'][method]={'twelve_molecule_batch_ns':stats(r['result']['elapsed_ns']/r['result']['iterations'] for r in rr),
            'target_pass':g['admitted'],'max_absolute_error_angstrom':max(float(e['absolute_error_upper']) for e in g['errors']),
            'errors':g['errors'],'speed_ratio_admitted_against_all_tools':False}
    accuracy=read(WINDOW/'accuracy_cost/summary.json');assert accuracy['complete'];quiet(accuracy)
    n1=read(HERE/'n1_accuracy_measurements.json');n1=copy.deepcopy(n1)
    original_n1_hash=inputs[str((HERE/'n1_accuracy_measurements.json').relative_to(ROOT))]
    measured={r['id']:r for r in accuracy['runs']['sterimol']['rows']};assert len(measured)==len(n1['rows'])==39
    for row in n1['rows']:
        for m in row['methods']:
            new=next(r for r in measured[row['id']]['measurements'] if r['method']==m['method'])
            assert new['result']==m['result']
            m.update(new);m['median_ns']=median(new['elapsed_ns'])/new['iterations']
            m['ns_per_call']=stats(v/new['iterations'] for v in new['elapsed_ns'])
    n1['scope']='Final controlled cost replay; original independent B1 references and residuals unchanged'
    n1['original_validation_sha256']=original_n1_hash;n1['timing_source']=accuracy['runs']['sterimol']
    n1['command']=accuracy['runs']['sterimol']['command']
    n1['engine_sha256']=accuracy['runs']['sterimol']['binary_sha256']
    n1['input_sha256']=accuracy['runs']['sterimol']['input_sha256']
    write('final_n1_accuracy_measurements.json',n1)
    n2={'scope':'Final controlled subset, preserving all earlier 244 validation observations separately',
        'rows':accuracy['runs']['volume']['rows'],'complete':True}
    assert len(n2['rows'])==72
    write('final_n2_accuracy_measurements.json',n2)
    pareto=[]
    for row in n1['rows']:
        reference=(row['reference']['lower']+row['reference']['upper'])/2
        pareto += [{'case':row['id'],'quantity':'B1_angstrom','method':r['method'],'median_ns':r['median_ns'],
                    'absolute_error_upper':r['absolute_error_upper'],'relative_error_upper':r['absolute_error_upper']/abs(reference) if reference else None,
                    'certified_bound':None} for r in row['methods']]
    for r in n2['rows']:
        reference=mean(r['reference_interval'])
        pareto.append({'case':r['id'].split('/')[0],'quantity':'volume_angstrom3','method':r['id'].split('/')[1],
                       'median_ns':r['median_ns'],'absolute_error_upper':r['absolute_error_upper'],
                       'relative_error_upper':r['absolute_error_upper']/abs(reference) if reference else None,
                       'certified_bound':r['result']['absolute_error_bound_angstrom3']})
    for r in pareto:
        peers=[p for p in pareto if p['case']==r['case'] and p['quantity']==r['quantity']]
        r['empirical_pareto']=not any(p['median_ns']<=r['median_ns'] and p['absolute_error_upper']<=r['absolute_error_upper'] and
             (p['median_ns']<r['median_ns'] or p['absolute_error_upper']<r['absolute_error_upper']) for p in peers)
    write('final_accuracy_pareto.json',{'scope':'Per-input empirical frontiers; upper residual estimate is distinct from a certified error bound; no universal optimality', 'rows':pareto})
    result['accuracy_counts']={'sterimol_cases':39,'sterimol_configurations':234,'volume_geometries':len(set(r['id'].split('/')[0] for r in n2['rows'])),'volume_configurations':72}
    write('final_measurements.json',result)
    sections=['# Final controlled measurement tables\n',
      'Derived by `summarize_final.py` from final window 6 and the Sterimol supplement. '
      'Every sample, median, unscaled MAD and range is retained in [final_measurements.json](final_measurements.json). '
      'The figures below are ratios of medians. Eight alternating trials per side follow the declared warmup. '
      'No outliers were removed. Process startup uses warm file caches; there is no claim of controlled cache-cold disk performance.\n',
      '## Whole CLI, milliseconds\n',table(['Workload / threads','Original median','Final median','Final MAD','Final min–max','Speedup'],
      [[k,f"{v['baseline']['median_ns']/1e6:.6f}",f"{v['candidate']['median_ns']/1e6:.6f}",f"{v['candidate']['MAD_ns']/1e6:.6f}",
        f"{v['candidate']['minimum_ns']/1e6:.6f}–{v['candidate']['maximum_ns']/1e6:.6f}",f"{v['speedup']:.4f}×"] for k,v in result['cli'].items()]),
      '## 10,000-file throughput and parallel efficiency\n',table(['Threads','Original files/s','Final files/s','Original efficiency','Final efficiency'],
      [[t,*[f"{r[s]['files_per_second']:,.1f}" for s in ['baseline','candidate']],*[f"{100*r[s]['parallel_efficiency']:.2f}%" for s in ['baseline','candidate']]] for t,r in result['scaling'].items()]),
      '## Isolated paired kernels, microseconds\n',table(['Kernel','Original','Final','Final MAD','Speedup'],
      [[k,f"{v['baseline']['median_ns']/1e3:.3f}",f"{v['candidate']['median_ns']/1e3:.3f}",f"{v['candidate']['mad_ns']/1e3:.3f}",f"{v['speedup']:.4f}×"] for k,v in result['kernel_ab'].items()]),
      '## Python paired costs, microseconds per workflow\n',table(['Mode / threads','Original','Final','Final MAD','Speedup'],
      [[k,f"{v['baseline']['median']/1e3:.3f}",f"{v['candidate']['median']/1e3:.3f}",f"{v['candidate']['mad']/1e3:.3f}",f"{v['speedup']:.4f}×"] for k,v in result['python_ab'].items()]),
      '## Final native / bridge / public layers, microseconds\n',table(['Layer','Median','MAD','Min–max'],
      [[k,f"{v['median_ns']/1e3:.3f}",f"{v['mad_ns']/1e3:.3f}",f"{v['minimum_ns']/1e3:.3f}–{v['maximum_ns']/1e3:.3f}"] for k,v in result['python_layers'].items()]),
      '## Final descriptor/I/O kernel profile, microseconds\n',table(['Mode','Median','MAD','Min–max'],
      [[k,f"{v['ns_per_call']['median']/1e3:.3f}",f"{v['ns_per_call']['mad']/1e3:.3f}",f"{v['ns_per_call']['min']/1e3:.3f}–{v['ns_per_call']['max']/1e3:.3f}"] for k,v in result['kernels'].items()]),
      '## Matched volume workflow, milliseconds\n',
      'Warm workflow: each independent trial contributes its mean internal iteration time (100 iterations for one input; three for 56; one for 1,000). All inner samples are also retained. DBSTEP is the common-grid adapter, not the default DBSTEP workflow.\n',
      table(['Size / workers / scope / tool','Median','MAD','Min–max'],
      [[k,f"{v['time_ns']['median']/1e6:.4f}",f"{v['time_ns']['mad']/1e6:.4f}",f"{v['time_ns']['min']/1e6:.4f}–{v['time_ns']['max']/1e6:.4f}"] for k,v in result['competitive'].items()]),
      '## External 56-conformer scaling efficiency\n',table(['Scope / tool',*[str(t)+' workers' for t in [1,2,4,6,12]]],
      [[k,*[f"{100*r[str(t)]['efficiency']:.2f}%" for t in [1,2,4,6,12]]] for k,r in result['competitive_scaling'].items()]),
      '## Sterimol cost and error, twelve-molecule in-memory batch\n',
      'The fixed target is 0.01 Å on every L/B1/B5. Morfeus and DBSTEP each miss one B1 target; their timings are cost/error observations, excluded from equal-accuracy speed ratios. Setup/containers/imports are outside the timer; alignment and result construction are inside.\n',
      table(['Method','Batch median µs','MAD µs','Max residual Å','Target passed'],
      [[k,f"{v['twelve_molecule_batch_ns']['median']/1e3:.3f}",f"{v['twelve_molecule_batch_ns']['mad']/1e3:.3f}",f"{v['max_absolute_error_angstrom']:.9g}",v['target_pass']] for k,v in result['sterimol_competitive'].items()])]
    (HERE/'FINAL_TABLES.md').write_text('\n'.join(sections))
    print(json.dumps({'complete':True,'cli_trials':len(cli['attempts']),'competitive_trials':len(comp['runs']),
                     'sterimol_trials':len(st['runs']),'accuracy':result['accuracy_counts'],'profile_overhead_ratio':result['profile_overhead_ratio']}))

if __name__=='__main__':main()
