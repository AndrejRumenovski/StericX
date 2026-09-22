"""Create a NEW sibling evidence directory, using only retained published inputs.

Run from a checkout whose production sources match freeze.json. This checks all
production source hashes; source.tar.gz also contains the exact source snapshot.
No old evidence is overwritten and no hidden .stericx input is needed.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
from datetime import datetime,timezone
from freeze import HERE,ROOT,sha,write


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--python",default="python3.12")
    p.add_argument("--prepare-only",action="store_true",help="Rebuild and run all scientific/preflight gates, then stop before timing")
    args=p.parse_args()
    out=args.output.resolve()
    assert out.parent==HERE.parent,"Use a new sibling under docs/benchmarks for the relative Cargo dependency"
    frozen=json.loads((HERE/"freeze.json").read_text())
    for name,digest in frozen["source_sha256"].items(): assert sha(ROOT/name)==digest,name
    assert sha(ROOT/"scripts/profile_child.c")==sha(HERE/"build/profile_child.c"),"Changed measurement launcher source"
    out.mkdir(exist_ok=False)
    for directory in ["build","environment","raw","inputs","adapter"]: (out/directory).mkdir()
    for src in HERE.glob("*.py"): shutil.copyfile(src,out/src.name)
    shutil.copyfile(HERE/".gitignore",out/".gitignore")
    for src in (HERE/"adapter").rglob("*"):
        if src.is_file() and "target" not in src.relative_to(HERE/"adapter").parts:
            dst=out/src.relative_to(HERE);dst.parent.mkdir(exist_ok=True,parents=True);shutil.copyfile(src,dst)
    shutil.copyfile(HERE/"input_manifest.json",out/"input_manifest.json")
    manifest=json.loads((out/"input_manifest.json").read_text())
    for row in manifest["corpus"]:
        src=HERE/row["filename"];assert sha(src)==row["sha256"]
        dst=out/row["filename"];dst.parent.mkdir(exist_ok=True,parents=True);shutil.copyfile(src,dst)
    shutil.copyfile(HERE/"environment/requirements.txt",out/"environment/requirements.txt")
    shutil.copyfile(HERE/"environment/reference-sources.stdout",out/"environment/reference-sources.stdout")
    env=os.environ.copy()
    for k in list(env):
        if k.startswith("CARGO_PROFILE_") or k in ["RUSTFLAGS","CARGO_ENCODED_RUSTFLAGS","RUSTC_WRAPPER","RUSTC_WORKSPACE_WRAPPER","CARGO_BUILD_TARGET"]: del env[k]
    def run(name,cmd):
        with (out/"raw"/(name+".stdout")).open("wb") as stdout,(out/"raw"/(name+".stderr")).open("wb") as stderr:
            result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=stdout,stderr=stderr)
        write(out/"raw"/(name+".command.json"),dict(argv=cmd,returncode=result.returncode))
        result.check_returncode()
    python=str(out/"environment/venv/bin/python")
    run("venv",["uv","venv","--python",args.python,str(out/"environment/venv")])
    run("install",["uv","pip","install","--python",python,"--index-url","https://pypi.org/simple","-r",str(out/"environment/requirements.txt")])
    build=["cargo","build","--locked","--release","--bin","stericx","--target-dir",str(out/"build/target")]
    run("build",build)
    shutil.copy2(out/"build/target/release/stericx",out/"build/stericx")
    run("adapter-build",["cargo","build","--locked","--release","--manifest-path",str(out/"adapter/Cargo.toml")])
    frozen.update(created_utc=datetime.now(timezone.utc).isoformat(),build_command=build,executable_sha256=sha(out/"build/stericx"),reproduced_from=str(HERE),original_freeze_sha256=sha(HERE/"freeze.json"),affinity_inherited=sorted(os.sched_getaffinity(0)))
    write(out/"freeze.json",frozen)
    for name,cmd in {"cpu":["lscpu","-J"],"kernel":["uname","-a"],"rustc":["rustc","-Vv"],"python":[python,"-VV"],"packages":["uv","pip","freeze","--python",python]}.items():run(name,cmd)
    for name,cmd in [
        ("compare",[python,str(out/"compare.py")]),
        ("exact-gate",[python,str(out/"exact_volume_gate.py")]),
        ("prepare",[python,str(out/"benchmark.py"),"prepare"]),
        ("timing",[python,str(out/"benchmark.py"),"run"]),
        ("summary",[python,str(out/"benchmark.py"),"summarize"]),
        ("sterimol-diagnosis",[python,str(out/"diagnose_sterimol.py")]),
        ("pyramid",[python,str(out/"pyramid_benchmark.py")]),
        ("memory",[python,str(out/"memory.py")]),
        ("protocol-review",[python,str(out/"protocol_review.py")]),
    ]:
        run(name,cmd)
        if args.prepare_only and name=="prepare":break
    print(out)


if __name__=="__main__":main()
