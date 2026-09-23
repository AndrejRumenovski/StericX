"""Gated native E2E timing: one warmup, four alternating AB/BA pairs.

Every stdout/stderr and measurement is preserved. No outlier removal.
Multiprocessing RSS is explicitly not represented by wait4's ru_maxrss.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
from datetime import datetime, timezone

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
AFFINITIES={1:[2],2:[2,3],4:[0,1,2,3],6:[0,1,2,3,4,5]}
ENV={k:"1" for k in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"]}


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text())
def write(p,obj): Path(p).write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
def now(): return datetime.now(timezone.utc).isoformat()


def verify_freeze():
    frozen=read(HERE/"freeze.json")
    for name,digest in frozen["source_sha256"].items(): assert sha(ROOT/name)==digest,name
    assert sha(HERE/"build/stericx")==frozen["executable_sha256"]
    gate=read(HERE/"scientific_gate.json")
    assert gate["passed"] and not gate["failures"] and gate["numerical_tolerance"] is None
    return frozen,gate


def argv(tool,workers,paths):
    if tool=="stericx": return [str(HERE/"adapter/target/release/volume-batch"),*paths]
    return [str(HERE/"environment/venv/bin/python"),str(HERE/"morfeus_batch.py"),"--workers",str(workers),*paths]


def canonical(rows,fields):
    return [{k:r[k] for k in fields} for r in rows]


def prepare():
    frozen,gate=verify_freeze()
    assert not (HERE/"preflight.json").exists()
    manifest=read(HERE/"input_manifest.json")
    paths=[str(HERE/r["filename"]) for r in manifest["corpus"]]
    for path,row in zip(paths,manifest["corpus"],strict=True): assert sha(path)==row["sha256"]
    # Same exact original repeated workload, materialized with the original names.
    batch=HERE/"inputs/repeated_10000"
    batch.mkdir(exist_ok=True)
    repeated=[]
    for row in manifest["repeated_10000"]:
        p=batch/row["filename"]
        shutil.copyfile(paths[row["corpus_index"]],p)
        assert sha(p)==row["sha256"]
        repeated.append(str(p))
    fields=gate["descriptors"]
    observations=[json.loads(line) for line in (HERE/"raw/stericx-observer.stdout").read_text().splitlines()]
    refs=read(HERE/"raw/morfeus-reference.json")["source_float64"]
    expected={}
    checks=[]
    # Compare the actual timed adapters over all 56 geometries at every worker count.
    for workers in AFFINITIES:
        for tool in ["stericx","morfeus"]:
            command=["taskset","-c",",".join(map(str,AFFINITIES[workers])),*argv(tool,workers,paths)]
            p=subprocess.run(command,env=os.environ|ENV|{"RAYON_NUM_THREADS":str(workers)},capture_output=True)
            stem=HERE/"raw"/f"preflight-{tool}-{workers}"
            stem.with_suffix(".stdout").write_bytes(p.stdout)
            stem.with_suffix(".stderr").write_bytes(p.stderr)
            write(stem.with_suffix(".command.json"),dict(argv=command,returncode=p.returncode,timing=False))
            p.check_returncode()
            rows=json.loads(p.stdout)
            assert len(rows)==56 and [r["file"] for r in rows]==paths
            values=canonical(rows,fields)
            if tool not in expected:
                expected[tool]=values
                if tool=="stericx":
                    import struct
                    for r,o in zip(rows,observations,strict=True):
                        for k in fields: assert struct.pack("f",r[k])==struct.pack("f",o[k]),(tool,k)
                else:
                    for r,ref in zip(rows,refs,strict=True):
                        for k in fields:
                            v=ref["values"][k]
                            if k=="percent_buried_volume":
                                import math
                                v=100*ref["values"]["buried_volume"]/(4*math.pi*3.5**3/3)
                            assert r[k]==v,(tool,k,r[k],v)
            assert values==expected[tool]
            checks.append(dict(tool=tool,workers=workers,rows=56,passed=True))
    launcher=HERE/"build/profile_child.c"
    shutil.copyfile(ROOT/"scripts/profile_child.c",launcher)
    p=subprocess.run(["cc","-O2","-o",str(HERE/"build/measure-child"),str(launcher)],capture_output=True)
    (HERE/"raw/launcher-build.stdout").write_bytes(p.stdout)
    (HERE/"raw/launcher-build.stderr").write_bytes(p.stderr)
    p.check_returncode()
    write(HERE/"workloads.json",{"56":paths,"1000":repeated[:1000],"10000":repeated})
    identities={str(p.relative_to(HERE)):sha(p) for p in [HERE/"adapter/target/release/volume-batch",HERE/"morfeus_batch.py",HERE/"build/measure-child",HERE/"scientific_gate.json",HERE/"benchmark.py"]}
    write(HERE/"preflight.json",dict(passed=True,created_utc=now(),checks=checks,expected=expected,identities=identities,worker_counts=list(AFFINITIES),pairs=4,warmups_per_tool=1))


def distribution(values):
    median=statistics.median(values)
    return dict(median=median,MAD=statistics.median(abs(x-median) for x in values),minimum=min(values),maximum=max(values),samples=values)


def run():
    verify_freeze()
    preflight=read(HERE/"preflight.json")
    assert preflight["passed"]
    for p,digest in preflight["identities"].items(): assert sha(HERE/p)==digest,p
    directory=HERE/"timings"
    directory.mkdir(exist_ok=False)
    workloads=read(HERE/"workloads.json")
    fields=read(HERE/"scientific_gate.json")["descriptors"]
    journal=dict(started_utc=now(),status="running",samples=[],pairs=4,warmups_per_tool=1,exclusions=[])
    write(HERE/"raw_timings.json",journal)
    try:
        for size,paths in workloads.items():
            for workers,cpus in AFFINITIES.items():
                config=directory/f"n{size}-w{workers}"
                config.mkdir()
                commands={tool:["taskset","-c",",".join(map(str,cpus)),str(HERE/"build/measure-child"),"METRICS",*argv(tool,workers,paths)] for tool in ["stericx","morfeus"]}
                with gzip.open(config/"commands.json.gz","wt") as f: json.dump(commands,f)
                for pair in [-1,0,1,2,3]:
                    order=["stericx","morfeus"] if pair<0 or pair%2==0 else ["morfeus","stericx"]
                    for tool in order:
                        phase="warmup" if pair<0 else f"pair{pair}"
                        stem=config/f"{phase}-{tool}"
                        command=list(commands[tool]);command[4]=str(stem.with_suffix(".metrics.json"))
                        attempt=dict(tool=tool,workers=workers,structures=int(size),phase=phase,pair=pair,affinity=cpus,started_utc=now(),order=order,status="running")
                        journal["samples"].append(attempt);write(HERE/"raw_timings.json",journal)
                        with stem.with_suffix(".stdout").open("wb") as out,stem.with_suffix(".stderr").open("wb") as err:
                            p=subprocess.run(command,env=os.environ|ENV|{"RAYON_NUM_THREADS":str(workers)},stdout=out,stderr=err)
                        metrics=read(stem.with_suffix(".metrics.json"))
                        attempt.update(metrics=metrics,status="measured",stdout_sha256=sha(stem.with_suffix(".stdout")),stderr_sha256=sha(stem.with_suffix(".stderr")))
                        write(HERE/"raw_timings.json",journal)
                        p.check_returncode()
                        rows=read(stem.with_suffix(".stdout"))
                        assert len(rows)==int(size)
                        for i,row in enumerate(rows):
                            assert row["file"]==paths[i]
                            assert {k:row[k] for k in fields}==preflight["expected"][tool][i%56]
                        content=stem.with_suffix(".stdout").read_bytes()
                        stem.with_suffix(".stdout.gz").write_bytes(gzip.compress(content,mtime=0))
                        assert gzip.decompress(stem.with_suffix(".stdout.gz").read_bytes())==content
                        stem.with_suffix(".stdout").unlink()
                        attempt.update(status="passed",raw_prefix=str(stem.relative_to(HERE)))
                        write(HERE/"raw_timings.json",journal)
                        print(f"n{size} w{workers} {phase} {tool}: {metrics['wall_ns']/1e9:.4f} s",flush=True)
        verify_freeze()
        for p,digest in preflight["identities"].items(): assert sha(HERE/p)==digest,p
        journal.update(status="complete",finished_utc=now())
    except BaseException as e:
        journal.update(status="failed",error=str(e));raise
    finally: write(HERE/"raw_timings.json",journal)


def summarize():
    raw=read(HERE/"raw_timings.json")
    assert raw["status"]=="complete"
    summaries=[]
    for size in [56,1000,10000]:
        for workers in AFFINITIES:
            rows=[s for s in raw["samples"] if s["structures"]==size and s["workers"]==workers and s["pair"]>=0]
            tools={}
            for tool in ["stericx","morfeus"]:
                samples=[s["metrics"] for s in rows if s["tool"]==tool]
                tools[tool]={k:distribution([s[k] for s in samples]) for k in ["wall_ns","cpu_user_s","cpu_system_s","peak_rss_bytes"]}
                tools[tool]["cpu_total_s"]=distribution([s["cpu_user_s"]+s["cpu_system_s"] for s in samples])
                tools[tool]["molecules_per_second"]=size/(tools[tool]["wall_ns"]["median"]/1e9)
            ratios=[]
            for pair in range(4):
                sample={s["tool"]:s["metrics"]["wall_ns"] for s in rows if s["pair"]==pair}
                ratios.append(sample["morfeus"]/sample["stericx"])
            summaries.append(dict(structures=size,workers=workers,tools=tools,paired_speedup=distribution(ratios),ratio_of_wall_medians=tools["morfeus"]["wall_ns"]["median"]/tools["stericx"]["wall_ns"]["median"]))
    write(HERE/"benchmark_results.json",dict(status="complete",scientific_gate=read(HERE/"scientific_gate.json"),freeze=read(HERE/"freeze.json"),summaries=summaries,raw_timings="raw_timings.json",memory_caveat="wait4 peak_rss_bytes is a single-process/maximum-descendant high-water mark, NOT total multiprocessing memory. See separate process-tree samples.",measurement_gate=all(x["paired_speedup"]["minimum"]>1 for x in summaries)))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("phase",choices=["prepare","run","summarize"])
    args=p.parse_args();globals()[args.phase]()
