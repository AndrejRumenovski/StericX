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
    p.add_argument("output",type=Path)
    p.add_argument("--python",default="3.12.13")
    p.add_argument("--prepare-only",action="store_true",help="Rebuild and run all scientific/preflight gates, then stop before timing")
    args=p.parse_args()
    out=args.output.resolve()
    assert out.parent==HERE.parent,"Use a new sibling under docs/benchmarks for the relative Cargo dependency"
    frozen=json.loads((HERE/"freeze.json").read_text())
    for name,digest in frozen["source_sha256"].items(): assert sha(ROOT/name)==digest,name
    assert sha(ROOT/"scripts/profile_child.c")==sha(HERE/"build/profile_child.c"),"Changed measurement launcher source"
    out.mkdir(exist_ok=False)
    for directory in ["build","environment","raw","inputs","adapter","revision_policy"]: (out/directory).mkdir()
    for src in HERE.glob("*.py"): shutil.copyfile(src,out/src.name)
    shutil.copyfile(HERE/".gitignore",out/".gitignore")
    shutil.copyfile(HERE/"EQUIVALENCE_POLICY.md",out/"EQUIVALENCE_POLICY.md")
    for name in ["reference.py","evaluate.py","selftest.py","grid_append.rs","run_admitted.py","cpython-3.12.13-bltinmodule.c","reference_revision.json"]:
        shutil.copyfile(HERE/"revision_policy"/name,out/"revision_policy"/name)
    write(out/"revision_policy/policy_freeze.json",dict(policy="EQUIVALENCE_POLICY.md",sha256=sha(out/"EQUIVALENCE_POLICY.md"),frozen_utc=datetime.now(timezone.utc).isoformat(),prior_residual_exposure=True,prior_timings_admissible=False,classification_not_yet_run=True,reproduced_policy=True))
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
    run("interpreter-check",[python,"-c","import sys; assert sys.version_info[:3] == (3,12,13), sys.version"])
    run("reference-selftest",[python,str(out/"revision_policy/selftest.py")])
    shutil.copyfile(out/"raw/reference-selftest.stdout",out/"revision_policy/selftest.json")
    run("scientific-evaluation",[python,str(out/"revision_policy/evaluate.py")])
    gate=json.loads((out/"revision_policy/evaluation_2/scientific_gate.json").read_text())
    assert gate["passed"] and not gate["failures"], "Scientific policy did not admit this reproduction; no timing authorized"
    if not args.prepare_only:
        run("admitted-campaign",[python,str(out/"revision_policy/run_admitted.py")])
    print(out)


if __name__=="__main__":main()
