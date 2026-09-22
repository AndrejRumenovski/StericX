"""Freeze current source, executable, reference environment and existing workload."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p, obj):
    Path(p).write_text(json.dumps(obj, indent=2) + "\n")


def command(name, argv, env=None):
    p = subprocess.run(argv, cwd=ROOT, capture_output=True, env=env)
    base = HERE / "environment" / name
    base.with_suffix(".stdout").write_bytes(p.stdout)
    base.with_suffix(".stderr").write_bytes(p.stderr)
    write(base.with_suffix(".command.json"), dict(argv=argv, returncode=p.returncode))
    p.check_returncode()
    return p.stdout.decode()


def main():
    for name in ["environment", "build", "inputs", "raw"]:
        (HERE / name).mkdir(exist_ok=True)
    assert not (HERE / "freeze.json").exists(), "Freeze already exists"
    commit = command("git-commit", ["git", "rev-parse", "HEAD"]).strip()
    command("git-status", ["git", "status", "--short"])
    names = subprocess.check_output(["git", "ls-files", "src", "Cargo.toml", "Cargo.lock"], cwd=ROOT, text=True).splitlines()
    source = {n: sha(ROOT / n) for n in names}
    with tarfile.open(HERE / "build/source.tar.gz", "w:gz") as tar:
        for n in names:
            tar.add(ROOT / n, arcname=n)
    env = os.environ.copy()
    for k in list(env):
        if k.startswith("CARGO_PROFILE_") or k in ["RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "RUSTC_WRAPPER", "RUSTC_WORKSPACE_WRAPPER", "CARGO_BUILD_TARGET"]:
            del env[k]
    build = ["cargo", "build", "--locked", "--offline", "--release", "--bin", "stericx", "--target-dir", str(HERE / "build/target")]
    command("build", build, env)
    shutil.copy2(HERE / "build/target/release/stericx", HERE / "build/stericx")
    for name, argv in {
        "rustc": ["rustc", "-Vv"], "cargo": ["cargo", "-V"],
        "rust-target-cfg": ["rustc", "--print", "cfg", "-C", "opt-level=3"],
        "cpu": ["lscpu", "-J"], "cpu-topology": ["lscpu", "-e=CPU,CORE,SOCKET,ONLINE"],
        "kernel": ["uname", "-a"], "uv": ["uv", "--version"],
    }.items():
        command(name, argv)
    for name in ["/proc/meminfo", "/proc/cpuinfo", "/etc/os-release"]:
        shutil.copyfile(name, HERE / "environment" / Path(name).name)
    metadata = urllib.request.urlopen("https://pypi.org/pypi/morfeus-ml/json").read()
    (HERE / "environment/pypi-morfeus.json").write_bytes(metadata)
    version = json.loads(metadata)["info"]["version"]
    python = str(HERE / "environment/venv/bin/python")
    command("venv-create", ["uv", "venv", "--python", str(ROOT / ".venv/bin/python"), str(HERE / "environment/venv")])
    command("install", ["uv", "pip", "install", "--python", python, "--index-url", "https://pypi.org/simple", "morfeus-ml==" + version])
    requirements = command("pip-freeze", ["uv", "pip", "freeze", "--python", python])
    (HERE / "environment/requirements.txt").write_text(requirements)
    command("python", [python, "-VV"])
    command("numpy-config", [python, "-c", "import numpy; numpy.show_config()"])
    command("reference-sources", [python, "-c", "import morfeus,pathlib,hashlib,json; p=pathlib.Path(morfeus.__file__).parent; print(json.dumps({str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(p.rglob('*.py'))},indent=2))"])
    wpath = ROOT / ".stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2/workloads/manifest.json"
    assert sha(wpath) == "a1041cc07fc5583c33d2a7cf584dbc638866c64bdd0b64e3720ef83cc85e1e6a"
    workloads = json.loads(wpath.read_text())["workloads"]
    paths = [Path(p) for p in workloads["conformers_56"]["argv"] if p.endswith(".xyz")]
    assert len(paths) == 56
    entries = []
    for p in paths:
        dest = HERE / "inputs" / p.parent.name / p.name
        dest.parent.mkdir(exist_ok=True)
        shutil.copyfile(p, dest)
        lines = p.read_text().splitlines()
        atoms = [line.split()[0] for line in lines[2:]]
        assert len(atoms) == int(lines[0]) and atoms.count("P") == 1
        entries.append(dict(filename=str(dest.relative_to(HERE)), source=str(p.relative_to(ROOT)), ligand=p.parent.name, conformer=p.stem, sha256=sha(p), atom_count=len(atoms), donor_element="P", donor_index_zero_based=atoms.index("P")))
    # Preserve all original 10,000 filenames/hashes/order, not 10,000 new chemistries.
    repeated = []
    for i, name in enumerate(p for p in workloads["descriptors_10000"]["argv"] if p.endswith(".xyz")):
        p = Path(name)
        assert sha(p) == entries[i % 56]["sha256"]
        repeated.append(dict(filename=p.name, source=str(p.relative_to(ROOT)), corpus_index=i % 56, sha256=sha(p), atom_count=entries[i % 56]["atom_count"], donor_index_zero_based=entries[i % 56]["donor_index_zero_based"]))
    assert len(repeated) == 10000
    write(HERE / "input_manifest.json", dict(source_manifest_sha256=sha(wpath), corpus=entries, repeated_10000=repeated, structures_1000="First 1000 entries of repeated_10000; optional timing workload", config=dict(sphere_radius=3.5, density=0.01, center_distance=2.28, radii_scale=1.17, include_hydrogens_volume=False, include_hydrogens_sterimol=True, sterimol_axis="coordination", L_correction=0.4, volume_planes=3)))
    write(HERE / "freeze.json", dict(created_utc=datetime.now(timezone.utc).isoformat(), stericx_commit=commit, source_sha256=source, executable_sha256=sha(HERE / "build/stericx"), build_command=build, release=dict(opt_level=3, lto="thin", codegen_units=1, cargo_features=[], target_cpu="rustc default x86-64", rustflags=None), morfeus_version=version, installation="Fresh uv virtual environment; PyPI morfeus-ml latest version resolved then pinned; environment/install.command.json", affinity_inherited=sorted(os.sched_getaffinity(0)), planned_affinities={1:[2],2:[2,3],4:[0,1,2,3],6:[0,1,2,3,4,5]}, input_manifest_sha256=sha(HERE / "input_manifest.json")))


if __name__ == "__main__":
    main()
