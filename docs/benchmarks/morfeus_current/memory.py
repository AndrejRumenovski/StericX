"""Separate process-tree RSS sampling; not used for headline wall times.

RSS is summed at each snapshot, never by adding individual peaks. Shared pages
are counted once per process, so this is not unique memory/PSS. Sampling can miss
short peaks. Every snapshot, process name, PID, RSS and HWM is retained.
"""
import gzip
import json
import os
from pathlib import Path
import subprocess
import time
from benchmark import HERE,AFFINITIES,ENV,argv,read,write,verify_freeze,sha,now


def snapshot(root_pid):
    todo=[root_pid]
    seen=set()
    records=[]
    while todo:
        pid=todo.pop()
        if pid in seen: continue
        seen.add(pid)
        try:
            status=Path(f"/proc/{pid}/status").read_text()
            fields=dict(line.split(":",1) for line in status.splitlines() if ":" in line)
            records.append(dict(pid=pid,name=fields["Name"].strip(),rss_bytes=int(fields.get("VmRSS","0 kB").split()[0])*1024,hwm_bytes=int(fields.get("VmHWM","0 kB").split()[0])*1024))
            for task in Path(f"/proc/{pid}/task").iterdir():
                try: todo.extend(map(int,(task/"children").read_text().split()))
                except FileNotFoundError: pass
        except (FileNotFoundError,ProcessLookupError): pass
    return records


def main():
    verify_freeze()
    assert read(HERE/"raw_timings.json")["status"]=="complete"
    preflight=read(HERE/"preflight.json")
    for p,h in preflight["identities"].items(): assert sha(HERE/p)==h
    outdir=HERE/"memory"
    outdir.mkdir(exist_ok=False)
    paths=read(HERE/"workloads.json")["10000"]
    fields=read(HERE/"scientific_gate.json")["descriptors"]
    results=[]
    for workers,cpus in AFFINITIES.items():
        for tool in ["stericx","morfeus"]:
            stem=outdir/f"{tool}-w{workers}"
            command=["taskset","-c",",".join(map(str,cpus)),*argv(tool,workers,paths)]
            write(stem.with_suffix(".command.json"),dict(argv=command,started_utc=now(),sample_interval_target_s=0.01))
            samples=[]
            with stem.with_suffix(".stdout").open("wb") as out,stem.with_suffix(".stderr").open("wb") as err:
                started=time.monotonic_ns()
                p=subprocess.Popen(command,env=os.environ|ENV|{"RAYON_NUM_THREADS":str(workers)},stdout=out,stderr=err)
                while p.poll() is None:
                    before=time.monotonic_ns()
                    processes=snapshot(p.pid)
                    samples.append(dict(elapsed_ns=before-started,processes=processes,sum_rss_bytes=sum(x["rss_bytes"] for x in processes),sample_duration_ns=time.monotonic_ns()-before))
                    time.sleep(0.01)
                assert p.returncode==0
            content=stem.with_suffix(".stdout").read_bytes()
            rows=json.loads(content)
            assert len(rows)==len(paths)
            for i,r in enumerate(rows):
                assert r["file"]==paths[i]
                assert {k:r[k] for k in fields}==preflight["expected"][tool][i%56]
            stem.with_suffix(".stdout.gz").write_bytes(gzip.compress(content,mtime=0))
            stem.with_suffix(".stdout").unlink()
            with gzip.open(stem.with_suffix(".samples.json.gz"),"wt") as f: json.dump(samples,f)
            peak=max(samples,key=lambda s:s["sum_rss_bytes"])
            result=dict(tool=tool,workers=workers,structures=10000,sampled_peak_sum_rss_bytes=peak["sum_rss_bytes"],peak_snapshot=peak,samples=len(samples),max_sample_spacing_s=max((b["elapsed_ns"]-a["elapsed_ns"])/1e9 for a,b in zip(samples,samples[1:])),stdout_sha256=__import__("hashlib").sha256(content).hexdigest())
            results.append(result)
            write(HERE/"memory_results.json",dict(results=results,complete=False))
            print(f"{tool} w{workers}: {peak['sum_rss_bytes']/2**20:.2f} MiB summed RSS",flush=True)
    write(HERE/"memory_results.json",dict(results=results,complete=True,scope="Separate untimed-for-speed 10000-file launches. Sum concurrent process RSS including parent, workers and resource tracker; double counts shared pages; sampling can miss short peaks. Never sum independent high-water marks."))


if __name__=="__main__": main()
