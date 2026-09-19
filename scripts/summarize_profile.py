"""Summarize native timings and optional probes without mixing worker and wall time."""

from __future__ import annotations

import argparse
import gzip
import json
import statistics
from collections import defaultdict
from pathlib import Path


def median(values: list[float]) -> float:
    return statistics.median(values)


def amdahl(fraction: float, improvement: float | None = None) -> float | None:
    """An infinite local improvement removes that fraction of measured wall time."""
    if not 0 <= fraction <= 1:
        raise ValueError(f"invalid wall-time fraction: {fraction}")
    denominator = 1 - fraction
    if improvement is not None:
        if improvement < 1:
            raise ValueError("local improvement must be >= 1")
        denominator += fraction / improvement
    return 1 / denominator if denominator else None


def summarize(root: Path, baseline: str, diagnostic: str) -> dict:
    native_path = root / "runs" / baseline / "summary.json"
    probe_path = root / "runs" / diagnostic / "summary.json"
    native = json.loads(native_path.read_text())
    probes = json.loads(probe_path.read_text())
    if native["status"] != "complete" or probes["status"] != "complete":
        raise ValueError("only completed measurement series can be summarized")
    if native["manifest_sha256"] != probes["manifest_sha256"]:
        raise ValueError("native and diagnostic inputs differ")
    if (
        native["affinity"] != probes["affinity"]
        or native["environment"] != probes["environment"]
    ):
        raise ValueError("native and diagnostic execution settings differ")
    results = {}
    for name, baseline_result in native["results"].items():
        probe_result = probes["results"][name]
        if baseline_result["fingerprints"] != probe_result["fingerprints"]:
            raise ValueError(f"scientific output changed: {name}")
        files = sorted(
            (root / "runs" / diagnostic / name).glob("measured_*/stages.json")
        )
        if len(files) != probe_result["measured_runs"]:
            raise ValueError(f"missing stage reports: {name}")
        stage_samples = defaultdict(list)
        function_samples = defaultdict(list)
        allocation_samples = defaultdict(list)
        residuals = []
        thread_counts = []
        worker_samples = defaultdict(list)
        worker_function_samples = defaultdict(list)
        for sample_index, file in enumerate(files):
            report = json.loads(file.read_text())
            metrics = json.loads(file.with_name("metrics.json").read_text())
            if report["dropped_scopes"]:
                raise ValueError(f"scope capacity exceeded: {file}")
            functions = report["functions"]
            command_roots = [f for f in functions if f["function"] == "cli::run"]
            if len(command_roots) != 1:
                raise ValueError(f"expected exactly one CLI wall root: {file}")
            command_root = command_roots[0]
            main_thread = command_root["thread_index"]
            main_functions = [f for f in functions if f["thread_index"] == main_thread]
            wall_ns = metrics["wall_ns"]
            root_ns = command_root["inclusive_ns"]
            if sum(f["exclusive_ns"] for f in main_functions) != root_ns:
                raise ValueError(
                    f"serial scopes do not partition CLI wall time: {file}"
                )
            if root_ns > wall_ns:
                raise ValueError(
                    f"internal scope exceeds measured process wall: {file}"
                )
            residuals.append(wall_ns - root_ns)
            thread_counts.append(report["observed_threads"])
            by_stage = defaultdict(int)
            seen_functions = set()
            for function in main_functions:
                by_stage[function["stage"]] += function["exclusive_ns"]
                key = (function["stage"], function["function"])
                seen_functions.add(key)
                if key not in function_samples:
                    function_samples[key] = [
                        dict(
                            exclusive_ns=0,
                            inclusive_ns=0,
                            calls=0,
                            fraction=0,
                            inclusive_fraction=0,
                        )
                        for _ in range(sample_index)
                    ]
                function_samples[key].append(
                    {
                        "exclusive_ns": function["exclusive_ns"],
                        "inclusive_ns": function["inclusive_ns"],
                        "calls": function["calls"],
                        "fraction": function["exclusive_ns"] / wall_ns,
                        "inclusive_fraction": function["inclusive_ns"] / wall_ns,
                    }
                )
            for key in function_samples.keys() - seen_functions:
                function_samples[key].append(
                    dict(
                        exclusive_ns=0,
                        inclusive_ns=0,
                        calls=0,
                        fraction=0,
                        inclusive_fraction=0,
                    )
                )
            by_stage["outside_cli_scope"] = wall_ns - root_ns
            for stage in stage_samples.keys() | by_stage.keys():
                if stage not in stage_samples:
                    stage_samples[stage] = [(0, 0)] * sample_index
                duration = by_stage.get(stage, 0)
                stage_samples[stage].append((duration, duration / wall_ns))
            workers = defaultdict(int)
            worker_functions = defaultdict(
                lambda: dict(exclusive_ns=0, inclusive_ns=0, calls=0)
            )
            for function in functions:
                if function["thread_index"] != main_thread:
                    workers[function["stage"]] += function["exclusive_ns"]
                    key = (function["stage"], function["function"])
                    for field in ("exclusive_ns", "inclusive_ns", "calls"):
                        worker_functions[key][field] += function[field]
            for stage in worker_samples.keys() | workers.keys():
                if stage not in worker_samples:
                    worker_samples[stage] = [0] * sample_index
                worker_samples[stage].append(workers.get(stage, 0))
            for key in worker_function_samples.keys() | worker_functions.keys():
                zero = dict(exclusive_ns=0, inclusive_ns=0, calls=0)
                if key not in worker_function_samples:
                    worker_function_samples[key] = [zero] * sample_index
                worker_function_samples[key].append(worker_functions.get(key, zero))
            for key, value in report["allocations"].items():
                allocation_samples[key].append(value)

        stages = {}
        for stage, samples in stage_samples.items():
            fraction = median([v[1] for v in samples])
            stages[stage] = {
                "median_exclusive_ns": median([v[0] for v in samples]),
                "median_diagnostic_wall_fraction": fraction,
                "maximum_end_to_end_speedup": amdahl(fraction),
                "end_to_end_speedup_if_stage_twice_as_fast": amdahl(fraction, 2),
            }
        ranked_functions = []
        for (stage, function), samples in function_samples.items():
            fraction = median([v["fraction"] for v in samples])
            inclusive_fraction = median([v["inclusive_fraction"] for v in samples])
            ranked_functions.append(
                {
                    "stage": stage,
                    "function": function,
                    "median_calls": median([v["calls"] for v in samples]),
                    "median_exclusive_ns": median([v["exclusive_ns"] for v in samples]),
                    "median_inclusive_ns": median([v["inclusive_ns"] for v in samples]),
                    "median_diagnostic_wall_fraction": fraction,
                    "maximum_end_to_end_speedup_removing_self_work": amdahl(fraction),
                    "median_inclusive_diagnostic_wall_fraction": inclusive_fraction,
                    "maximum_end_to_end_speedup_removing_function_and_callees": amdahl(
                        inclusive_fraction
                    ),
                }
            )
        ranked_functions.sort(key=lambda f: f["median_exclusive_ns"], reverse=True)
        native_wall = baseline_result["summary"]["wall_ns"]["median"]
        diagnostic_wall = probe_result["summary"]["wall_ns"]["median"]
        results[name] = {
            "native": baseline_result,
            "diagnostic_process_wall_ns": diagnostic_wall,
            "diagnostic_to_native_wall_ratio": diagnostic_wall / native_wall,
            "main_thread_stage_wall": stages,
            "main_thread_functions": ranked_functions,
            "worker_summed_elapsed_ns_by_stage": {
                k: median(v) for k, v in worker_samples.items()
            },
            "worker_functions_summed_elapsed": sorted(
                [
                    {
                        "stage": stage,
                        "function": function,
                        **{
                            f"median_{field}": median([v[field] for v in samples])
                            for field in ("exclusive_ns", "inclusive_ns", "calls")
                        },
                    }
                    for (stage, function), samples in worker_function_samples.items()
                ],
                key=lambda f: f["median_exclusive_ns"],
                reverse=True,
            ),
            "observed_threads": sorted(set(thread_counts)),
            "requested_heap": {k: median(v) for k, v in allocation_samples.items()},
            "outside_cli_scope_median_ns": median(residuals),
        }
    return {
        "schema_version": 1,
        "baseline": baseline,
        "diagnostic": diagnostic,
        "manifest_sha256": native["manifest_sha256"],
        "baseline_binary_sha256": native["binary_sha256"],
        "diagnostic_binary_sha256": probes["binary_sha256"],
        "timing_method": (
            "Native timing is authoritative. Function and stage times are measured "
            "in a separate diagnostic build. Its main-thread exclusive times "
            "partition the CLI root; worker totals are shown separately and never "
            "added to wall shares. Missing scope samples count as zero."
        ),
        "amdahl_method": (
            "Per-run diagnostic fraction f = exclusive stage/function ns divided "
            "by external process wall ns. Report the median f; infinite local "
            "optimization has ceiling 1/(1-f), and 2x local improvement gives "
            "1/(1-f+f/2). These are diagnostic-workload bounds, not measured native "
            "speedups. Instrumentation overhead is reported per workload. Function "
            "self-work and inclusive ceilings are separate; nested inclusive "
            "ceilings must never be added. Outside-CLI residual includes loading, "
            "exit, launcher wait and diagnostic report writing, and is not a "
            "single optimizable function."
        ),
        "allocation_method": (
            "Requested System heap counters include process startup and diagnostic "
            "session metadata and exclude report construction. Requested bytes "
            "include new sizes of reallocations. Heap highwater excludes mappings, "
            "stacks, static arrays and allocator metadata; native RSS is measured "
            "separately."
        ),
        "results": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(".stericx/profiling"))
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--diagnostic", required=True)
    parser.add_argument("--output", type=Path, default=Path("docs/profiling"))
    args = parser.parse_args()
    report = summarize(args.root.resolve(), args.baseline, args.diagnostic)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "profile_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    manifest = args.root / "workloads/manifest.json"
    (args.output / "workloads.json.gz").write_bytes(
        gzip.compress(manifest.read_bytes(), mtime=0)
    )
    for name, result in report["results"].items():
        print(
            name, round(result["native"]["summary"]["wall_ns"]["median"] / 1e6, 3), "ms"
        )


if __name__ == "__main__":
    main()
