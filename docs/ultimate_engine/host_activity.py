#!/usr/bin/env python3
"""Measure interval host activity before the final timing window; never stop other apps."""
import argparse,json,os,time
from pathlib import Path
from benchmark_artifacts import drain_artifacts
def snapshot():
    cpu={}
    for line in Path('/proc/stat').read_text().splitlines():
        key,*values=line.split()
        if key.startswith('cpu'):cpu[key]=list(map(int,values))
    processes={}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            raw=(p/'stat').read_text();name=raw[raw.index('(')+1:raw.rindex(')')]
            fields=raw[raw.rindex(')')+2:].split()
            processes[p.name]=(name,int(fields[11])+int(fields[12]),fields[19])
        except (OSError,ValueError,IndexError):continue
    return cpu,processes
def measure(seconds=10.):
    a,pa=snapshot();start=time.monotonic();time.sleep(seconds);b,pb=snapshot();elapsed=time.monotonic()-start
    hz=os.sysconf('SC_CLK_TCK');d=[v-u for u,v in zip(a['cpu'],b['cpu'])]
    # Guest time is already included in user/nice; do not count it twice.
    busy=sum(d[:8])-d[3]-d[4]
    rows=[]
    for pid,(name,ticks,identity) in pb.items():
        if pid in pa and pa[pid][2]==identity and ticks>pa[pid][1]:
            rows.append({'pid':int(pid),'name':name,'core_equivalents':(ticks-pa[pid][1])/hz/elapsed})
    rows.sort(key=lambda r:r['core_equivalents'],reverse=True)
    result={'interval_seconds':elapsed,'busy_core_equivalents':busy/hz/elapsed,
        'io_wait_core_equivalents':d[4]/hz/elapsed,'processes':rows[:20],
        'pressure_cpu':Path('/proc/pressure/cpu').read_text(),'pressure_io':Path('/proc/pressure/io').read_text(),
        'criteria':'At most 0.15 busy CPU cores and 0.05 cores of I/O wait averaged over this interval. Conservative local campaign gate, not a hardware universal.'}
    result['quiet']=result['busy_core_equivalents']<=.15 and result['io_wait_core_equivalents']<=.05
    return result
def postflight(directory):
    flushed=drain_artifacts(directory)
    result=measure()
    result['artifact_flush']=flushed
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--seconds',type=float,default=10)
    a=p.parse_args();assert 1<=a.seconds<=30
    value=measure(a.seconds)
    with a.output.open('x') as f:json.dump(value,f,indent=2);f.write('\n')
    print(json.dumps(value,indent=2))
