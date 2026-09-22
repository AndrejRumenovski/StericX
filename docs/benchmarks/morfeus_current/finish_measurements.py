raise SystemExit("Archived pre-policy campaign: publication claim withdrawn. See BENCHMARK.md and EQUIVALENCE_POLICY.md.")
"""Sequential continuation; never build/reference-check alongside native timing."""
import json
from pathlib import Path
import subprocess
import time
from benchmark import HERE,write


def main():
    while True:
        try:state=json.loads((HERE/"raw_timings.json").read_text())["status"]
        except json.JSONDecodeError:time.sleep(1);continue
        if state=="complete":break
        if state=="failed":raise RuntimeError("Primary volume timing failed")
        time.sleep(1)
    python=str(HERE/"environment/venv/bin/python")
    commands=[
        ("volume-summary",[python,str(HERE/"benchmark.py"),"summarize"]),
        ("pyramid-build",["cargo","build","--offline","--locked","--release","--manifest-path",str(HERE/"adapter/Cargo.toml"),"--bin","pyramid-batch"]),
        ("protocol-review",[python,str(HERE/"protocol_review.py")]),
        ("sterimol-diagnosis",[python,str(HERE/"diagnose_sterimol.py")]),
        ("pyramid-run",[python,str(HERE/"pyramid_benchmark.py")]),
        ("memory-run",[python,str(HERE/"memory.py")]),
    ]
    for name,cmd in commands:
        stem=HERE/"raw"/name
        with stem.with_suffix(".stdout").open("wb") as stdout,stem.with_suffix(".stderr").open("wb") as stderr:
            p=subprocess.run(cmd,stdout=stdout,stderr=stderr)
        write(stem.with_suffix(".command.json"),dict(argv=cmd,returncode=p.returncode))
        p.check_returncode()
        print(name+" complete",flush=True)


if __name__=="__main__":main()
