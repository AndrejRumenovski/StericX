#!/usr/bin/env python3
"""Gated matched-workflow preflight, timing, or separate process-tree memory sampling."""
import argparse,gzip,json,os,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import profile_stericx as h
from host_activity import measure as host_measure
BASE=ROOT/'.stericx/ultimate_engine'
AFFINITY={1:[2],2:[2,3],4:list(range(4)),6:list(range(6)),12:list(range(12))}

def snapshot(pid):
    todo=[pid];seen=set();rows=[]
    while todo:
        pid=todo.pop()
        if pid in seen:continue
        seen.add(pid)
        try:
            fields=dict(s.split(':',1) for s in Path(f'/proc/{pid}/status').read_text().splitlines() if ':' in s)
            rows.append({'pid':pid,'name':fields['Name'].strip(),'rss_bytes':int(fields.get('VmRSS','0 kB').split()[0])*1024})
            for task in Path(f'/proc/{pid}/task').iterdir():
                try:todo.extend(map(int,(task/'children').read_text().split()))
                except OSError:pass
        except OSError:pass
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['preflight','timing','memory'])
    p.add_argument('--gate',type=Path,default=BASE/'competitive_gate/summary.json')
    p.add_argument('--preflight',type=Path,default=BASE/'competitive_preflight/summary.json')
    p.add_argument('--native',type=Path,default=BASE/'competitive_warm/target/release/volume_workflow')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    gate=a.gate.resolve();g=json.loads(gate.read_text());assert g['passed'] and g['complete']
    pre=json.loads(a.preflight.read_text());assert pre['complete']
    native=a.native.resolve()
    worker=Path(__file__).with_name('competitive_worker.py')
    binaries=[native,worker,gate,a.preflight,Path(__file__)]
    shutil.copy2(Path(__file__),out/Path(__file__).name);shutil.copy2(worker,out/worker.name)
    receipt={'phase':a.phase,'scope':'56 unique geometries; repeated paths form the 1000-calculation batch. Matched finite-grid volume only.',
        'bindings':{str(p):h.digest(p) for p in binaries},'runs':[],'complete':False}
    def save():(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if a.phase=='timing':
        receipt['host_preflight']=host_measure();save()
        if not receipt['host_preflight']['quiet']:raise SystemExit('No timing: host requires quiet window.')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    env.pop('STERICX_PROFILE_PATH',None)
    launcher=h.build_launcher(out/'launcher') if a.phase=='timing' else None
    configs=[(56,t) for t in AFFINITY] if a.phase=='preflight' else [(1000,t) for t in AFFINITY] if a.phase=='memory' else [(1,1),*[(56,t) for t in AFFINITY],(1000,1),(1000,6)]
    for size,threads in configs:
        source=[g['rows'][i%56] for i in range(size)];paths=[r['path'] for r in source]
        for scope in (['fresh_process','warm_workflow'] if a.phase=='timing' else ['preflight'] if a.phase=='preflight' else ['memory']):
            for pair in (range(-1,8) if a.phase=='timing' else [0]):
                order=['native','morfeus','dbstep']
                if pair%2:order.reverse()
                for tool in order:
                    iterations=(100 if size==1 else 3 if size==56 else 1) if scope=='warm_workflow' else 1 if a.phase=='preflight' else 0
                    command=[str(native),str(iterations),*paths] if tool=='native' else [str(BASE.parent/'paton/reference/venv/bin/python'),str(worker),'--tool',tool,'--workers',str(threads),'--warm-iterations',str(iterations),*paths]
                    folder=out/f'n{size}_t{threads}_{scope}_{pair}_{tool}'
                    row={'size':size,'threads':threads,'scope':scope,'pair':pair,'warmup':pair<0,'tool':tool,'status':'started'}
                    receipt['runs'].append(row);save();child_env=dict(env,RAYON_NUM_THREADS=str(threads))
                    if a.phase=='timing':
                        row['metrics']=h.measure(command,child_env,set(AFFINITY[threads]),folder,launcher)
                        assert row['metrics'].get('exit_code',row['metrics'].get('returncode',0))==0
                    else:
                        folder.mkdir();samples=[]
                        launched=['taskset','-c',','.join(map(str,AFFINITY[threads])),*command]
                        with (folder/'stdout.txt').open('w') as stream,(folder/'stderr.txt').open('w') as err:
                            proc=subprocess.Popen(launched,env=child_env,stdout=stream,stderr=err)
                            while proc.poll() is None:
                                if a.phase=='memory':
                                    start=time.monotonic_ns();processes=snapshot(proc.pid)
                                    samples.append({'time_ns':start,'processes':processes,'sum_rss_bytes':sum(p['rss_bytes'] for p in processes),'sampling_ns':time.monotonic_ns()-start})
                                time.sleep(.01)
                        assert proc.returncode==0,(tool,size,threads)
                        (folder/'command.json').write_text(json.dumps(launched))
                        if samples:
                            with gzip.open(folder/'memory_samples.json.gz','wt') as f:json.dump(samples,f)
                            row['sampled_peak_sum_rss_bytes']=max(s['sum_rss_bytes'] for s in samples)
                            row['maximum_sample_spacing_seconds']=max((b['time_ns']-a['time_ns'])/1e9 for a,b in zip(samples,samples[1:])) if len(samples)>1 else None
                    value=json.loads((folder/'stdout.txt').read_text())
                    assert value['rows']==[r[tool] for r in source],(tool,size,threads)
                    row.update(status='passed',warm_elapsed_ns=value['warm_elapsed_ns'],stdout_sha256=h.digest(folder/'stdout.txt'))
                    save()
    for path,digest in receipt['bindings'].items():assert h.digest(Path(path))==digest,path
    if a.phase=='timing':receipt['host_postflight']=host_measure()
    receipt['memory_caveat']='Summed concurrent RSS double-counts shared pages and sampling can miss peaks. wait4 max RSS is not total process-tree memory. Memory sampling is never used for speed claims.'
    receipt['complete']=True;save()

if __name__=='__main__':main()
