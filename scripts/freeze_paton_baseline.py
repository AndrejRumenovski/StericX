"""Seal task baseline evidence without changing scientific inputs or outputs.

Run with the existing research environment after the source-bound build.
The archived complete observations are reused only after exact scientific-source
identity and output bytes are verified. A fresh replay is recorded separately.
"""

import hashlib
import importlib.metadata
import json
import platform
import re
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def record(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def command(args):
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    return {
        "argv": args,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def main():
    destination = ROOT / "docs/paton_recommendations/baseline/manifest.json"
    if destination.exists():
        raise FileExistsError("baseline seal already exists; never overwrite it")
    local = ROOT / ".stericx/paton/baseline"
    build = local / "build/manifest.json"
    provider = ROOT / "docs/kraken_residual_forensics/candidates/angular_sampling"
    observations = provider / "full_observations/manifest.json"
    identical_sources = []
    for path in [
        *sorted((ROOT / "src").rglob("*.rs")),
        ROOT / "Cargo.toml",
        ROOT / "Cargo.lock",
    ]:
        old = provider / "validated_build/source" / path.relative_to(ROOT)
        if path.read_bytes() != old.read_bytes():
            raise ValueError(f"cannot reuse observations: different source {path}")
        identical_sources.append(record(path))
    # Require completed provider observations and a fresh source-bound build.
    for path in (build, observations):
        if not path.is_file():
            raise FileNotFoundError(path)
    tracked = (
        subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT)
        .decode()
        .split("\0")
    )
    sources = []
    reports = []
    for name in sorted(filter(None, tracked)):
        path = ROOT / name
        if name.startswith(
            ("src/", "tests/", "scripts/", "studies/", ".github/")
        ) or name in {
            "Cargo.toml",
            "Cargo.lock",
            "pyproject.toml",
            "uv.lock",
        }:
            sources.append(record(path))
        if name.startswith("docs/") and path.suffix == ".md" and "/frozen/" not in name:
            reports.append(record(path))
    receipts = [record(build), record(observations)]
    provider_data = json.loads(observations.read_text())
    for name, lane in provider_data["lanes"].items():
        path = observations.parent / name / "stdout.jsonl"
        if record(path)["sha256"] != lane["stdout"]["sha256"]:
            raise ValueError(f"changed archived observations in {name}")
    for path in sorted(observations.parent.rglob("*")):
        if path.is_file():
            receipts.append(record(path))
    for path in (local / "build/bin").iterdir():
        receipts.append(record(path))
    receipts.append(record(local / "stericx"))
    receipts.append(record(Path(__file__)))
    for path in sorted(destination.parent.glob("*.log")):
        if path.name != "observe.log":  # fresh replay is still live, sealed later
            receipts.append(record(path))
    rust_log = (destination.parent / "rust_tests.log").read_text()
    python_log = (destination.parent / "python_tests.log").read_text()
    rust_counts = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed", rust_log)
    python_count = re.search(r"Ran (\d+) tests", python_log)
    if (
        not rust_counts
        or any(int(failed) for _, failed in rust_counts)
        or "test result: FAILED" in rust_log
        or python_count is None
        or "\nOK\n" not in python_log
    ):
        raise ValueError("engineering logs do not establish passing baseline tests")
    manifest = {
        "kind": "paton_task_baseline",
        "created_utc": datetime.now(UTC).isoformat(),
        "git_commit": command(["git", "rev-parse", "HEAD"]),
        "git_tree": command(["git", "rev-parse", "HEAD^{tree}"]),
        "tracked_changes": command(["git", "diff", "--stat"]),
        "rustc": command(["rustc", "-Vv"]),
        "cargo": command(["cargo", "-V"]),
        "python": {"version": sys.version, "executable": sys.executable},
        "platform": platform.platform(),
        "dependencies": sorted(
            [
                {"name": d.metadata["Name"], "version": d.version}
                for d in importlib.metadata.distributions()
            ],
            key=lambda d: d["name"].lower(),
        ),
        "sources": sources,
        "report_inventory": reports,
        "identical_scientific_sources_with_observation_provider": identical_sources,
        "observation_provider": record(observations),
        "fresh_replay_status_at_freeze": "running separately on frozen source",
        "execution_artifacts": receipts,
        "engineering": {
            "rust_tests": sum(int(passed) for passed, _ in rust_counts),
            "python_tests": int(python_count[1]),
            "failures": 0,
        },
        "scientific_global_pass": False,
        "scope": (
            "Frozen baseline observations, including scientific errors; "
            "no new validation claim"
        ),
        "preserved_limitations": [
            "Missing historical Kraken energies and exact conformer identity",
            "B1 finite angular phase; finite-grid/f32 occupancy and argmin sensitivity",
            "Radii, center and regional-normalization convention differences",
            "XYZ inferred connectivity and ambiguous planar centers",
            "Force-field/electronic weights are not verified equilibrium populations",
            "Failed folds, poor predictive validation and Ni-hDA 2064 provenance",
            "No DBSTEP or electron-density validation at baseline; no Python extension",
        ],
    }
    destination.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(record(destination)))


if __name__ == "__main__":
    main()
