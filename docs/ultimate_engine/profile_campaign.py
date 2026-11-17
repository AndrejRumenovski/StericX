"""Fresh baseline measurements. Every attempt is retained; no production edits.

Run from the repository root with .venv/bin/python. This is exploratory profiling,
not the final quiet-window A/B acceptance benchmark.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import sys
from benchmark_artifacts import compact_metrics, retain_fingerprint
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import profile_stericx as harness  # noqa: E402


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def host():
    return {name: Path("/proc", name).read_text() for name in
            ["loadavg", "pressure/cpu", "pressure/io", "meminfo"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--phase", choices=["cli", "kernels", "valgrind"], required=True)
    args = parser.parse_args()
    build, out = args.build.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, LC_ALL="C", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    env.pop("STERICX_PROFILE_PATH", None)
    launcher = harness.build_launcher(out / "launcher")
    journal = {"purpose": __doc__, "argv": sys.argv, "host_before": host(),
               "controller_sha256": harness.digest(Path(__file__)),
               "artifact_helper_sha256": harness.digest(Path(__file__).with_name("benchmark_artifacts.py")),
               "runs": []}

    def run(name, argv, threads, affinity, instrumented=False):
        directory = out / name
        child_env = dict(env, RAYON_NUM_THREADS=str(threads))
        if instrumented:
            child_env["STERICX_PROFILE_PATH"] = str(directory / "profile.json")
        argv = [str(a).replace("{run_dir}", str(directory)) for a in argv]
        metrics = harness.measure(argv, child_env, set(affinity), directory, launcher)
        row = {"name": name, "threads": threads,
               "affinity": affinity, "instrumented": instrumented,
               "metrics": compact_metrics(metrics, directory)}
        journal["runs"].append(row)
        save(out / "journal.json", journal)
        print(name, metrics.get("exit_code", metrics.get("returncode")), flush=True)
        if metrics.get("exit_code", metrics.get("returncode", 0)) != 0:
            raise RuntimeError(f"failed measurement retained: {directory}")
        return directory, row

    if args.phase == "cli":
        path = ROOT / ".stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2/workloads/manifest.json"
        manifest = json.loads(path.read_text())
        harness.verify_manifest(manifest)
        journal["workload_manifest"] = {"path": str(path), "sha256": harness.digest(path)}
        native = build / "native-target/release/stericx"
        diagnostic = build / "diagnostic-target/release/stericx"
        journal["binaries"] = {str(p): harness.digest(p) for p in [native, diagnostic]}
        expected = {}
        for threads, affinity in [(1, [2]), (2, [2, 3]), (4, [0, 1, 2, 3]),
                                  (6, list(range(6))), (12, list(range(12)))]:
            for name in ["single_small", "single_large", "conformers_56", "ensemble_sdf",
                         "descriptors_10000", "parse_ensembles_1000", "db_build_ensembles"]:
                work = manifest["workloads"][name]
                for rep in range(4):
                    directory, row = run(f"t{threads}/{name}/native_{rep}",
                                         [native, *work["argv"]], threads, affinity)
                    row["warmup"] = rep == 0
                    scientific = harness.check_result(work, directory)
                    expected.setdefault(name, scientific)
                    if expected[name] != scientific:
                        raise RuntimeError(f"nonidentical outputs: {name}")
                    row["scientific"] = retain_fingerprint(scientific, out / "fingerprints" / (name + ".json"))
                directory, row = run(f"t{threads}/{name}/diagnostic", [diagnostic, *work["argv"]],
                                     threads, affinity, True)
                scientific = harness.check_result(work, directory)
                if expected[name] != scientific:
                    raise RuntimeError(f"instrumentation changed outputs: {name}")
                row["scientific"] = retain_fingerprint(scientific, out / "fingerprints" / (name + ".json"))
        harness.verify_manifest(manifest)
    elif args.phase == "kernels":
        molecule = ROOT / "python/examples/data/ligand_401_ensemble.xyz"
        options = ROOT / "python/examples/data/ligand_401_options.json"
        cube = ROOT / "python/tests/data/Ne_medium.cube"
        journal["inputs"] = {str(p): harness.digest(p) for p in [molecule, options, cube]}
        native, diagnostic = build / "bin/observer", build / "bin/observer_alloc"
        journal["binaries"] = {str(p): harness.digest(p) for p in [native, diagnostic]}
        modes = {"sterimol": 2000, "volume": 100, "regions": 100,
                 "pyramidalization": 10000, "center": 10000, "snapshot": 2000,
                 "radius_profile": 20, "ensemble": 20, "statistics": 1000,
                 "serialize": 1000, "parse": 1000, "density": 100,
                 "density_profile": 30, "density_parse": 5}
        for kind, iterations in modes.items():
            input_path = cube if kind.startswith("density") else molecule
            expected = None
            for rep in range(4):
                binary = diagnostic if rep == 3 else native
                directory, row = run(f"{kind}/{rep}", [binary, kind, input_path, options, iterations], 1, [2])
                result = json.loads((directory / "stdout.txt").read_text())
                row["kernel"] = {k: v for k, v in result.items() if k != "result"}
                row["result_sha256"] = fingerprint(result["result"])
                expected = expected or row["result_sha256"]
                if expected != row["result_sha256"]:
                    raise RuntimeError(f"kernel result changed: {kind}")
        for threads, affinity in [(1, [2]), (2, [2, 3]), (4, [0, 1, 2, 3]),
                                  (6, list(range(6))), (12, list(range(12)))]:
            for rep in range(3):
                directory, row = run(f"batch/t{threads}/{rep}",
                    [native, "batch", molecule, options, 40, threads], threads, affinity)
                result = json.loads((directory / "stdout.txt").read_text())
                row["kernel"] = {k: v for k, v in result.items() if k != "result"}
                row["result_sha256"] = fingerprint(result["result"])
    else:
        prefix = Path("/tmp/stericx-ultimate-valgrind")
        env["VALGRIND_LIB"] = str(prefix / "usr/libexec/valgrind")
        valgrind = prefix / "usr/bin/valgrind.bin"
        binary = build / "bin/observer"
        command = [binary, "ensemble", ROOT / "python/examples/data/ligand_401_ensemble.xyz",
                   ROOT / "python/examples/data/ligand_401_options.json", "5"]
        journal["binaries"] = {str(p): harness.digest(p) for p in [valgrind, binary]}
        for tool, flags in [
            ("callgrind", ["--callgrind-out-file={run_dir}/callgrind.out"]),
            ("cachegrind", ["--cache-sim=yes", "--branch-sim=yes", "--I1=32768,8,64",
                "--D1=32768,8,64", "--LL=16777216,16,64", "--cachegrind-out-file={run_dir}/cachegrind.out"]),
            ("dhat", ["--dhat-out-file={run_dir}/dhat.json"]),
            ("memcheck", ["--error-exitcode=99", "--leak-check=full", "--errors-for-leak-kinds=definite"]),
        ]:
            directory, row = run(tool, [valgrind, f"--tool={tool}", *flags, *command], 1, [2])
            result = json.loads((directory / "stdout.txt").read_text())
            row["result_sha256"] = fingerprint(result["result"])
        if len({r["result_sha256"] for r in journal["runs"]}) != 1:
            raise RuntimeError("profiler outputs differ")
    journal["host_after"] = host()
    journal["complete"] = True
    save(out / "journal.json", journal)


if __name__ == "__main__":
    main()
