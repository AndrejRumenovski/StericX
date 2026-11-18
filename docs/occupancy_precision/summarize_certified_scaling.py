"""Derive fixed-accuracy batch throughput, CPU, memory and worker efficiency."""
import argparse,json,statistics
from pathlib import Path
from campaign import ROOT,DOCS,EVIDENCE,sha,save
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',required=True);a=p.parse_args()
directory=EVIDENCE/a.input;raw=json.loads((directory/'summary.json').read_text());assert raw['complete'] and not raw['preflight_only']
assert raw['host_preflight']['quiet'] and raw['host_postflight']['quiet']
assert len(raw['runs'])==45 and all(r['exact'] for r in raw['runs'])
def stats(values):
    values=list(values);m=statistics.median(values)
    return {'values':values,'median':m,'mad':statistics.median(abs(v-m) for v in values),'min':min(values),'max':max(values),'n':len(values)}
rows={}
for threads in [1,2,4,6,12]:
    selected=[r for r in raw['runs'] if r['threads']==threads and not r['warmup']];assert len(selected)==8
    rows[str(threads)]={'batch_ns':stats(r['elapsed_ns'] for r in selected),'process_metrics':{k:stats(r['metrics'][k] for r in selected) for k in ['wall_ns','cpu_user_s','cpu_system_s','peak_rss_bytes','minor_faults','major_faults','voluntary_context_switches','involuntary_context_switches']}}
for threads,row in rows.items():
    row['speedup']=rows['1']['batch_ns']['median']/row['batch_ns']['median'];row['efficiency']=row['speedup']/int(threads);row['molecules_per_second']=12e9/row['batch_ns']['median']
save(DOCS/'certified_scaling.json',{'complete':True,'scope':raw['scope'],'absolute_tolerance':.05,'molecules_per_batch':12,'rows':rows,'source_receipt':{'path':str(directory/'summary.json'),'sha256':sha(directory/'summary.json')},'controller_sha256':sha(Path(__file__))})
lines=['# Certified integration batch scaling','',raw['scope'],'','All 45 processes preserve every individual certificate and output order. Values below are medians of eight measured processes; full ranges, MAD, CPU and RSS distributions are in [certified_scaling.json](certified_scaling.json).','', '| Workers | Batch ms | MAD ms | Full range ms | Molecules/s | Speedup | Efficiency | Median process peak RSS MiB |','|---|---|---|---|---|---|---|---|']
for threads,row in rows.items():
    s=row['batch_ns'];lines.append(f"| {threads} | {s['median']/1e6:.3f} | {s['mad']/1e6:.3f} | {s['min']/1e6:.3f}–{s['max']/1e6:.3f} | {row['molecules_per_second']:.2f} | {row['speedup']:.4f}× | {100*row['efficiency']:.2f}% | {row['process_metrics']['peak_rss_bytes']['median']/1048576:.3f} |")
lines += ['', 'This is batch-level concurrency for independent molecules. Each certified integration remains serial, with its own candidates, heap and scratch buffers. The pool is created explicitly by the benchmark caller, so no nested pool or new default concurrency policy was added to production. Peak RSS includes the pool and warmup; more simultaneous adaptive heaps can increase memory. Twelve heterogeneous tasks impose a load-balance limit, and the result does not establish scaling of one large integration.','']
(DOCS/'CERTIFIED_SCALING.md').write_text('\n'.join(lines));print(json.dumps(rows,indent=2))
