"""Separate, exact-rounding-gated pyramidalization end-to-end experiment.

Run only AFTER volume measurements and builds have finished. Reuses the same
measurement launcher/protocol, and writes a separate raw sample collection.
"""
import gzip
import json
import os
import struct
import subprocess
from pathlib import Path
import benchmark as protocol
from benchmark import HERE,AFFINITIES,ENV,read,write,sha,now


def command(tool,workers,paths):
    if tool=="stericx":return [str(HERE/"adapter/target/release/pyramid-batch"),*paths]
    return [str(HERE/"environment/venv/bin/python"),str(HERE/"morfeus_pyramid.py"),"--workers",str(workers),*paths]


def main():
    protocol.verify_freeze()
    assert read(HERE/"raw_timings.json")["status"]=="complete"
    out=HERE/"pyramidalization"
    out.mkdir(exist_ok=False)
    workloads=read(HERE/"workloads.json")
    fields=["pyr_p","pyr_alpha"]
    observed=[json.loads(x) for x in (HERE/"raw/stericx-observer.stdout").read_text().splitlines()]
    expected={}
    checks=[]
    for workers,cpus in AFFINITIES.items():
        for tool in ["stericx","morfeus"]:
            cmd=["taskset","-c",",".join(map(str,cpus)),*command(tool,workers,workloads["56"])]
            p=subprocess.run(cmd,env=os.environ|ENV|{"RAYON_NUM_THREADS":str(workers)},capture_output=True)
            stem=out/f"preflight-{tool}-w{workers}"
            stem.with_suffix(".stdout").write_bytes(p.stdout);stem.with_suffix(".stderr").write_bytes(p.stderr)
            write(stem.with_suffix(".command.json"),dict(argv=cmd,returncode=p.returncode,timing=False))
            p.check_returncode()
            rows=json.loads(p.stdout)
            assert len(rows)==56 and [r["file"] for r in rows]==workloads["56"]
            vals=[{k:r[k] for k in fields} for r in rows]
            for a,b in zip(vals,observed,strict=True):
                for k in fields:assert struct.pack("f",a[k])==struct.pack("f",b[k]),(tool,k)
            if tool not in expected:expected[tool]=vals
            assert expected[tool]==vals
            checks.append(dict(tool=tool,workers=workers,passed=True))
    identities={str(p.relative_to(HERE)):sha(p) for p in [HERE/"adapter/target/release/pyramid-batch",HERE/"morfeus_pyramid.py",Path(__file__),HERE/"build/measure-child"]}
    gate=dict(passed=True,failures=[],numerical_tolerance=None,descriptors=fields,scope="All 56 conformers, identical f32-represented XYZ coordinates and donor neighbors; morfeus P and alpha rounded to f32 exactly equal every native output; zero tolerance.",checks=checks,identities=identities)
    write(out/"scientific_gate.json",gate)
    write(out/"freeze.json",read(HERE/"freeze.json"))
    journal=dict(status="running",started_utc=now(),samples=[],pairs=4,warmups_per_tool=1,exclusions=[])
    write(out/"raw_timings.json",journal)
    try:
        for size,paths in workloads.items():
            for workers,cpus in AFFINITIES.items():
                config=out/f"n{size}-w{workers}";config.mkdir()
                commands={tool:["taskset","-c",",".join(map(str,cpus)),str(HERE/"build/measure-child"),"METRICS",*command(tool,workers,paths)] for tool in ["stericx","morfeus"]}
                with gzip.open(config/"commands.json.gz","wt") as f:json.dump(commands,f)
                for pair in [-1,0,1,2,3]:
                    order=["stericx","morfeus"] if pair<0 or pair%2==0 else ["morfeus","stericx"]
                    for tool in order:
                        phase="warmup" if pair<0 else f"pair{pair}"
                        stem=config/f"{phase}-{tool}"
                        cmd=list(commands[tool]);cmd[4]=str(stem.with_suffix(".metrics.json"))
                        sample=dict(tool=tool,workers=workers,structures=int(size),phase=phase,pair=pair,affinity=cpus,started_utc=now(),order=order,status="running")
                        journal["samples"].append(sample);write(out/"raw_timings.json",journal)
                        with stem.with_suffix(".stdout").open("wb") as stdout,stem.with_suffix(".stderr").open("wb") as stderr:
                            p=subprocess.run(cmd,env=os.environ|ENV|{"RAYON_NUM_THREADS":str(workers)},stdout=stdout,stderr=stderr)
                        sample.update(metrics=read(stem.with_suffix(".metrics.json")),status="measured",stdout_sha256=sha(stem.with_suffix(".stdout")),stderr_sha256=sha(stem.with_suffix(".stderr")))
                        write(out/"raw_timings.json",journal);p.check_returncode()
                        content=stem.with_suffix(".stdout").read_bytes();rows=json.loads(content)
                        assert len(rows)==int(size)
                        for i,r in enumerate(rows):assert r["file"]==paths[i] and {k:r[k] for k in fields}==expected[tool][i%56]
                        stem.with_suffix(".stdout.gz").write_bytes(gzip.compress(content,mtime=0));stem.with_suffix(".stdout").unlink()
                        sample.update(status="passed",raw_prefix=str(stem.relative_to(HERE)))
                        write(out/"raw_timings.json",journal)
                        print(f"pyramid n{size} w{workers} {phase} {tool}: {sample['metrics']['wall_ns']/1e9:.4f} s",flush=True)
        protocol.verify_freeze()
        for p,h in identities.items():assert sha(HERE/p)==h
        journal.update(status="complete",finished_utc=now())
    except BaseException as e:journal.update(status="failed",error=str(e));raise
    finally:write(out/"raw_timings.json",journal)
    # The identical summary function consumes the same sample schema.
    protocol.HERE=out
    protocol.summarize()


if __name__=="__main__":main()
