"""Derive complete distributions and per-input accuracy frontiers from admitted runs."""
import argparse, collections, json, statistics
from pathlib import Path
from campaign import ROOT, EVIDENCE, DOCS, sha, save

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--window',action='append',required=True,help='Successful stage receipts; multiple retained windows allowed')
a=p.parse_args()
bindings={}; stages={}
def read(path):
    bindings[str(path.relative_to(ROOT))]=sha(path)
    return json.loads(path.read_text())
def quiet(value):
    if isinstance(value,dict):
        if 'quiet' in value:assert value['quiet'] is True
        for v in value.values():quiet(v)
    elif isinstance(value,list):
        for v in value:quiet(v)
def stats(values):
    values=list(values);m=statistics.median(values)
    return dict(values=values,n=len(values),median=m,mad=statistics.median(abs(v-m) for v in values),min=min(values),max=max(values))
def metrics(rows):
    return {key:stats(r['metrics'][key] for r in rows) for key in
            ['wall_ns','cpu_user_s','cpu_system_s','peak_rss_bytes','minor_faults','major_faults','voluntary_context_switches','involuntary_context_switches']}
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|',*['| '+' | '.join(map(str,r))+' |' for r in rows]])+'\n'
for name in a.window:
    directory=EVIDENCE/name;receipt=read(directory/'receipt.json')
    for row in receipt['runs']:
        if row['status']=='completed' and row.get('returncode')==0 and row.get('host_preflight',{}).get('quiet') and row.get('host_postflight',{}).get('quiet'):
            assert row['name'] not in stages,'Do not select repeated successful stages after seeing timings'
            stages[row['name']]=directory/row['name']
assert set(stages)=={'cli_ab','kernel_ab','chord_ab','kernel_profile','cli_profile','accuracy','competitive','radius_competitive'},stages
result={'complete':True,'scope':'All samples in prospectively gated final stages; no outlier removal. Ratios of medians, unscaled MAD, ranges, and complete distributions. Quiet endpoints cannot guarantee no interruption inside a stage.', 'bindings':bindings,'windows':a.window}
cli=read(stages['cli_ab']/'summary.json');assert cli['complete'];quiet(cli)
assert len(cli['attempts'])==630 and all(r['status']=='passed' for r in cli['attempts'])
result['cli']=cli['results']
for key,row in result['cli'].items():
    workload,t=key.rsplit('/t',1)
    for side in ['baseline','candidate']:
        rows=[r for r in cli['attempts'] if r['workload']==workload and r['threads']==int(t) and r['side']==side and not r['warmup']]
        assert len(rows)==8
        row[side]['process_metrics']=metrics(rows)
result['scaling']={}
for t in [1,2,4,6,12]:
    row=result['cli'][f'descriptors_10000/t{t}'];entry={'speedup':row['speedup']}
    for side in ['baseline','candidate']:
        ns=row[side]['median_ns'];speed=result['cli']['descriptors_10000/t1'][side]['median_ns']/ns
        entry[side]={'files_per_second':1e13/ns,'speedup_vs_one_thread':speed,'parallel_efficiency':speed/t}
    result['scaling'][str(t)]=entry
kernel=read(stages['kernel_ab']/'summary.json');assert kernel['complete']
result['kernel_ab']=kernel['summary']
chord=read(stages['chord_ab']/'summary.json');assert chord['complete']
result['chord_ab']=chord
kp=read(stages['kernel_profile']/'journal.json');assert kp['complete'] and len(kp['runs'])==71
grouped=collections.defaultdict(list)
for row in kp['runs']:grouped[row['name'].rsplit('/',1)[0]].append(row)
result['kernel_profile']={}
for name,rows in grouped.items():
    native=[r for r in rows if r['kernel']['allocations'] is None]
    result['kernel_profile'][name]={'ns_per_call':stats(r['kernel']['elapsed_ns']/r['kernel']['iterations'] for r in native), 'process_metrics':metrics(native),'allocations':[r for r in rows if r['kernel']['allocations'] is not None]}
cp=read(stages['cli_profile']/'journal.json');assert cp['complete'] and len(cp['runs'])==175
result['cli_profile']=cp
for stage in ['competitive','radius_competitive']:
    data=read(stages[stage]/'summary.json');assert data['complete'];quiet(data)
    assert len(data['runs'])==(432 if stage=='competitive' else 162)
    grouped=collections.defaultdict(list)
    for r in data['runs']:
        assert r.get('exact',r.get('status')=='passed')
        if not r['warmup']:grouped[(r.get('size',r.get('count')),r['threads'],r['scope'],r['tool'])].append(r)
    result[stage]={}
    for (n,t,scope,tool),rows in grouped.items():
        assert len(rows)==8
        result[stage][f'n{n}/t{t}/{scope}/{tool}']={'ns':stats(r['metrics']['wall_ns'] if scope=='fresh_process' else statistics.mean(r['warm_elapsed_ns']) for r in rows),'process_metrics':metrics(rows),'inner_samples':[r['warm_elapsed_ns'] for r in rows]}
accuracy=read(stages['accuracy']/'summary.json');assert accuracy['complete'];quiet(accuracy)
assert len(accuracy['runs'])==1152
memory=read(EVIDENCE/'accuracy_memory_v2/summary.json');assert memory['complete']
alloc={r['id']:r for r in memory['runs']};grouped=collections.defaultdict(list)
for row in accuracy['runs']:grouped[row['id']].append(row)
pareto=[]
for identity,rows in grouped.items():
    assert len(rows)==8
    r=rows[0];bound=r['result'].get('absolute_error_bound_angstrom3')
    pareto.append({'id':identity,'case':r['source_id'],'method':r['configuration'],
        'ns':stats(v for row in rows for v in row['elapsed_ns']), 'process_metrics':metrics(rows),
        'certified_error':bound,'observed_error_interval':r['observed_error_interval'],
        'reference_scope':r['reference_scope'],'result':r['result'],
        'allocations':alloc[identity]['allocations'],'diagnostic_process_metrics':alloc[identity]['metrics']})
for row in pareto:
    peers=[r for r in pareto if r['case']==row['case']]
    for mode in ['certified','observed']:
        error=lambda r:r['certified_error'] if mode=='certified' else r['observed_error_interval'][1]
        own=error(row)
        row[mode+'_pareto']=None if own is None else not any(error(r) is not None and r['ns']['median']<=row['ns']['median'] and error(r)<=own and (r['ns']['median']<row['ns']['median'] or error(r)<own) for r in peers)
frontiers=[]
for case in sorted({r['case'] for r in pareto}):
    for tolerance in [.1,.05,.01]:
        for mode in ['certified','observed']:
            admitted=[r for r in pareto if r['case']==case and (r['certified_error'] if mode=='certified' else r['observed_error_interval'][1]) is not None and (r['certified_error'] if mode=='certified' else r['observed_error_interval'][1])<=tolerance]
            fastest=min(admitted,key=lambda r:r['ns']['median']) if admitted else None
            frontiers.append({'case':case,'tolerance':tolerance,'accuracy_class':mode,'fastest_tested_method':fastest['method'] if fastest else None,'median_ns':fastest['ns']['median'] if fastest else None})
save(DOCS/'accuracy_pareto.json',{'scope':'Per-input tested frontier. Certified intrinsic bounds and observed error upper bounds are separate classes. Molecular references use a separate same-algorithm 0.001 certificate; analytical references are independent. No untested-method optimality claim.','rows':pareto,'fastest_by_target':frontiers,'bindings':bindings})
result['accuracy_rows']=pareto;result['fastest_by_target']=frontiers
save(DOCS/'final_measurements.json',result)
sections=['# Final controlled measurements\n','All trials, CPU/wall time, RSS, allocation diagnostics, MAD and ranges are in [final_measurements.json](final_measurements.json). Warmups remain in the raw receipts. No outliers were removed. Fresh processes use warmed file caches.\n',
    '## 10,000-file throughput\n',table(['Threads','Baseline files/s','Final files/s','Speedup','Final efficiency'],[[t,f"{r['baseline']['files_per_second']:,.1f}",f"{r['candidate']['files_per_second']:,.1f}",f"{r['speedup']:.4f}×",f"{100*r['candidate']['parallel_efficiency']:.2f}%"] for t,r in result['scaling'].items()]),
    '## Paired kernel latency\n',table(['Mode','Baseline µs','Final µs','Final MAD µs','Speedup'],[[k,f"{r['baseline']['median_ns']/1000:.3f}",f"{r['candidate']['median_ns']/1000:.3f}",f"{r['candidate']['mad_ns']/1000:.3f}",f"{r['speedup']:.4f}×"] for k,r in result['kernel_ab'].items()]),
    '## Whole-process workload distributions\n',table(['Workload/threads','Baseline ms','Final ms','Final MAD ms','Final range ms','Speedup'],[[k,f"{r['baseline']['median_ns']/1e6:.4f}",f"{r['candidate']['median_ns']/1e6:.4f}",f"{r['candidate']['MAD_ns']/1e6:.4f}",f"{r['candidate']['minimum_ns']/1e6:.4f}–{r['candidate']['maximum_ns']/1e6:.4f}",f"{r['speedup']:.4f}×"] for k,r in result['cli'].items()])]
for stage in ['competitive','radius_competitive']:
    sections += ['## '+stage.replace('_',' ')+'\n',table(['Size/workers/scope/tool','Median ms','MAD ms','Range ms'],[[k,f"{r['ns']['median']/1e6:.5f}",f"{r['ns']['mad']/1e6:.5f}",f"{r['ns']['min']/1e6:.5f}–{r['ns']['max']/1e6:.5f}"] for k,r in result[stage].items()])]
sections += ['## Accuracy and computational cost\n','Each row identifies its reference class in the linked JSON. An observed error against a fine same-algorithm molecular enclosure is not independent exact truth. Midpoint methods have no intrinsic certificate.\n',table(['Case/method','Median ms','MAD ms','Certified error Å³','Observed error upper Å³','Target met'],[[r['id'],f"{r['ns']['median']/1e6:.4f}",f"{r['ns']['mad']/1e6:.4f}",f"{r['certified_error']:.6g}" if r['certified_error'] is not None else 'uncertified',f"{r['observed_error_interval'][1]:.6g}",r['result'].get('target_met','not certified')] for r in pareto])]
(DOCS/'FINAL_TABLES.md').write_text('\n'.join(sections))
print(json.dumps({'complete':True,'scaling':result['scaling'],'accuracy_configurations':len(pareto)},indent=2))
