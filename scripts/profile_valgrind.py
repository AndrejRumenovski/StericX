#!/usr/bin/env python3
"""Run diagnostic Valgrind profiles against the frozen native workload oracle.

Cachegrind's cache/branch counts are simulations, not native hardware events;
its function costs are instruction counts, not elapsed time. DHAT heap sizes
are not RSS and its lifetime units are instructions. Run separately from native
timing experiments. Never scale or change a selected manifest workload.
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import statistics
import subprocess
import time
from pathlib import Path

from profile_stericx import (
    CONTROLLED_ENV,
    DEFAULT_ROOT,
    REPO,
    check_result,
    digest,
    now,
    verify_manifest,
)

CACHE_EVENTS = (
    "Ir",
    "I1mr",
    "ILmr",
    "Dr",
    "D1mr",
    "DLmr",
    "Dw",
    "D1mw",
    "DLmw",
    "Bc",
    "Bcm",
    "Bi",
    "Bim",
)
SIMULATION = {
    "native_hardware_counters": False,
    "function_cost_units": "executed instructions; not time",
    "cache_model": (
        "Split I1 and D1 with unified LL; write allocate, bit-indexed sets, "
        "approximately inclusive LL. Exact simulated configuration is retained "
        "from events.out desc lines. This is not a complete modern CPU model."
    ),
    "branch_model": (
        "Conditional: 16384 two-bit counters indexed by branch address/history. "
        "Indirect: 512-entry last-target predictor; returns assumed perfect. "
        "These simulated misses are not measurements of this CPU's predictor."
    ),
    "instruction_caveat": (
        "REP-prefixed instruction iterations are counted separately, unlike "
        "the corresponding native hardware instruction event."
    ),
    "documentation": "https://valgrind.org/docs/manual/cg-manual.html",
}


def save(path: Path, value: object) -> None:
    """Keep incremental status readable even if the controller is interrupted."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def fraction(numerator: int, denominator: int) -> dict:
    return {
        "numerator": numerator,
        "denominator": denominator,
        "fraction": numerator / denominator if denominator else None,
    }


def cache_geometry(value: str) -> dict:
    size, ways, line = map(int, value.split(","))
    if min(size, ways, line) <= 0 or size % (ways * line):
        raise ValueError("cache geometry requires size,ways,line with integral sets")
    sets = size // (ways * line)
    if sets & (sets - 1) or line & (line - 1):
        raise ValueError("cache sets and line size must be powers of two")
    return {"size_bytes": size, "ways": ways, "line_bytes": line}


def sysfs_caches(affinity: set[int]) -> dict:
    result = {}
    fields = (
        "level",
        "type",
        "size",
        "ways_of_associativity",
        "coherency_line_size",
        "number_of_sets",
        "shared_cpu_list",
    )
    for cpu in sorted(affinity):
        entries = []
        for directory in sorted(
            Path(f"/sys/devices/system/cpu/cpu{cpu}/cache").glob("index*")
        ):
            entry = {"sysfs_path": str(directory)}
            for field in fields:
                path = directory / field
                entry[field] = path.read_text().strip() if path.is_file() else None
            entries.append(entry)
        result[str(cpu)] = entries
    return result


def parse_cachegrind(path: Path) -> dict:
    descriptions, events, totals = [], None, None
    with path.open() as source:
        for line in source:
            if line.startswith("desc:"):
                descriptions.append(line.removeprefix("desc:").strip())
            elif line.startswith("events:"):
                if events is not None:
                    raise RuntimeError(f"duplicate events header: {path}")
                events = line.split()[1:]
            elif line.startswith("summary:"):
                if totals is not None:
                    raise RuntimeError(f"duplicate summary: {path}")
                totals = [int(value) for value in line.split()[1:]]
    if (
        events is None
        or totals is None
        or len(events) != len(totals)
        or len(set(events)) != len(events)
        or set(events) != set(CACHE_EVENTS)
        or any(value < 0 for value in totals)
    ):
        raise RuntimeError(f"incomplete or unsupported Cachegrind summary: {path}")
    counts = dict(zip(events, totals, strict=True))
    d_refs = counts["Dr"] + counts["Dw"]
    d1_misses = counts["D1mr"] + counts["D1mw"]
    ll_data_misses = counts["DLmr"] + counts["DLmw"]
    ll_requests = counts["I1mr"] + d1_misses
    ll_misses = counts["ILmr"] + ll_data_misses
    branch_refs = counts["Bc"] + counts["Bi"]
    branch_misses = counts["Bcm"] + counts["Bim"]
    return {
        "events_in_file_order": events,
        "counts": counts,
        "cache_configuration": descriptions,
        "derived_counts": {
            "data_references": d_refs,
            "D1_misses": d1_misses,
            "LL_requests_after_L1_misses": ll_requests,
            "LL_misses": ll_misses,
            "branches": branch_refs,
            "branch_mispredictions": branch_misses,
        },
        "simulated_rates": {
            "I1_misses_per_instruction": fraction(counts["I1mr"], counts["Ir"]),
            "D1_misses_per_data_reference": fraction(d1_misses, d_refs),
            "LL_misses_per_LL_request_local": fraction(ll_misses, ll_requests),
            "LL_misses_per_all_reference_global": fraction(
                ll_misses, counts["Ir"] + d_refs
            ),
            "LL_instruction_misses_per_instruction": fraction(
                counts["ILmr"], counts["Ir"]
            ),
            "LL_data_misses_per_data_reference": fraction(ll_data_misses, d_refs),
            "branch_mispredictions_per_branch": fraction(branch_misses, branch_refs),
            "conditional_mispredictions_per_conditional_branch": fraction(
                counts["Bcm"], counts["Bc"]
            ),
            "indirect_mispredictions_per_indirect_branch": fraction(
                counts["Bim"], counts["Bi"]
            ),
        },
    }


def parse_dhat(path: Path) -> dict:
    profile = json.loads(path.read_text())
    if (
        profile.get("dhatFileVersion") != 2
        or profile.get("mode") != "heap"
        or profile.get("tu") != "instrs"
        or not isinstance(profile.get("pps"), list)
    ):
        raise RuntimeError(f"incomplete or unsupported DHAT profile: {path}")
    points = profile["pps"]
    keys = ("tb", "tbk", "gb", "gbk", "eb", "ebk", "rb", "wb", "tl")
    totals = {key: sum(point[key] for point in points) for key in keys}
    if any(value < 0 for value in totals.values()):
        raise RuntimeError(f"negative DHAT counter: {path}")
    top = []
    for point in sorted(points, key=lambda item: item["tb"], reverse=True)[:30]:
        top.append(
            {
                "total_bytes": point["tb"],
                "total_blocks": point["tbk"],
                "bytes_at_global_heap_peak": point["gb"],
                "average_lifetime_instructions": (
                    point["tl"] / point["tbk"] if point["tbk"] else None
                ),
                "allocation_stack": [profile["ftbl"][frame] for frame in point["fs"]],
            }
        )
    return {
        "total_requested_bytes": totals["tb"],
        "total_allocation_blocks_including_realloc": totals["tbk"],
        "peak_live_heap_bytes": totals["gb"],
        "blocks_at_global_heap_peak": totals["gbk"],
        "end_live_heap_bytes": totals["eb"],
        "end_live_heap_blocks": totals["ebk"],
        "heap_read_bytes": totals["rb"],
        "heap_written_bytes": totals["wb"],
        "average_lifetime_instructions": (
            totals["tl"] / totals["tbk"] if totals["tbk"] else None
        ),
        "end_instructions": profile["te"],
        "global_heap_peak_at_instruction": profile["tg"],
        "program_points": len(points),
        "top_allocation_sites_by_requested_bytes": top,
        "semantics": (
            "Whole-process intercepted heap operations, not RSS. realloc adds a "
            "block and its entire requested new size; attribution remains at the "
            "original allocation site. Lifetimes/progress use instructions, not "
            "elapsed time. Site-local peak sizes are not summed for global peak."
        ),
        "documentation": "https://valgrind.org/docs/manual/dh-manual.html",
    }


def allocation_comparison(root: Path, label: str, name: str, dhat: dict) -> dict:
    """Compare scopes explicitly; libc and profiler metadata prevent equality."""
    directory = root / "runs" / label
    summary_path = directory / "summary.json"
    summary = json.loads(summary_path.read_text())
    if summary.get("status") != "complete":
        raise RuntimeError(f"allocator reference is incomplete: {summary_path}")
    profiles = sorted((directory / name).glob("measured_*/stages.json"))
    if not profiles:
        raise RuntimeError(f"no measured allocator reports: {directory / name}")
    rows = []
    for path in profiles:
        profile = json.loads(path.read_text())
        if profile.get("dropped_scopes") != 0:
            raise RuntimeError(f"incomplete allocator reference: {path}")
        rows.append(profile["allocations"])
    blocks = statistics.median(
        row["allocation_calls"] + row["reallocation_calls"] for row in rows
    )
    requested = statistics.median(row["requested_bytes"] for row in rows)
    peak = statistics.median(row["peak_live_requested_bytes"] for row in rows)
    return {
        "reference_label": label,
        "reference_manifest_sha256": summary["manifest_sha256"],
        "reference_summary_sha256": digest(summary_path),
        "reports": [{"path": str(path), "sha256": digest(path)} for path in profiles],
        "rust_allocator_median_allocation_plus_reallocation_calls": blocks,
        "rust_allocator_median_requested_bytes": requested,
        "rust_allocator_median_peak_live_requested_bytes": peak,
        "dhat_minus_rust_block_count": dhat["total_allocation_blocks_including_realloc"]
        - blocks,
        "dhat_minus_rust_requested_bytes": dhat["total_requested_bytes"] - requested,
        "dhat_minus_rust_peak_bytes": dhat["peak_live_heap_bytes"] - peak,
        "interpretation": (
            "Cross-check, not an equality gate: DHAT also intercepts libc heap "
            "allocations and includes shutdown; the Rust allocator probe covers "
            "System requests up to its pre-report snapshot and adds session metadata."
        ),
    }


def check_against_native(workload: dict, directory: Path, expected: dict) -> dict:
    """Reuse the native scientific oracle without tolerances or profiler noise."""
    fingerprints = check_result(workload, directory)
    if fingerprints != expected:
        raise RuntimeError(f"scientific output differs from native oracle: {directory}")
    return fingerprints


def terminate_group(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def run_child(
    argv: list[str],
    env: dict,
    affinity: set[int],
    stdout: Path,
    stderr: Path,
    started_callback,
    timeout: float | None,
) -> dict:
    with stdout.open("wb") as out, stderr.open("wb") as err:
        previous = os.sched_getaffinity(0)
        started = time.perf_counter_ns()
        try:
            os.sched_setaffinity(0, affinity)
            process = subprocess.Popen(
                argv,
                cwd=REPO,
                env=env,
                stdin=subprocess.DEVNULL,
                stdout=out,
                stderr=err,
                start_new_session=True,
            )
        finally:
            os.sched_setaffinity(0, previous)
        try:
            started_callback(process.pid)
            returncode = process.wait(timeout=timeout)
        except BaseException:
            terminate_group(process)
            raise
    return {
        "returncode": returncode,
        "profiler_controller_elapsed_ns_not_native_time": time.perf_counter_ns()
        - started,
    }


def run(args: argparse.Namespace) -> None:
    root = args.root.resolve()
    manifest_path = root / "workloads/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    verify_manifest(manifest)
    for value in (args.label, args.compare_to, args.allocator_reference):
        if value is not None and (
            not value or Path(value).name != value or value in (".", "..")
        ):
            raise RuntimeError("run labels must be single directory names")
    baseline_path = root / "runs" / args.compare_to / "summary.json"
    baseline = json.loads(baseline_path.read_text())
    manifest_hash = digest(manifest_path)
    if (
        baseline.get("status") != "complete"
        or baseline["manifest_sha256"] != manifest_hash
    ):
        raise RuntimeError(
            "native oracle must be complete and use the same frozen manifest"
        )
    default_cases = (
        list(manifest["workloads"])
        if args.tool == "cachegrind"
        else ["single_small", "conformers_56"]
    )
    selected = args.workloads.split(",") if args.workloads else default_cases
    if len(selected) != len(set(selected)) or any(
        name not in manifest["workloads"] or name not in baseline["results"]
        for name in selected
    ):
        raise RuntimeError("unknown, duplicate, or unbaselined workload")
    affinity = set(map(int, args.affinity.split(",")))
    if not affinity or not affinity <= os.sched_getaffinity(0) or args.threads < 1:
        raise RuntimeError("invalid affinity or thread count")
    if baseline["environment"]["RAYON_NUM_THREADS"] != str(args.threads):
        raise RuntimeError("oracle and profiler must use the same Rayon thread count")
    configured_caches = {
        "I1": cache_geometry(args.i1),
        "D1": cache_geometry(args.d1),
        "LL": cache_geometry(args.ll),
    }
    binary = (args.binary or root / "baseline/symbol-target/release/stericx").resolve(
        strict=True
    )
    capabilities_path = args.capabilities or root / "probes/capabilities.json"
    capabilities = json.loads(capabilities_path.read_text())
    valgrind = Path(capabilities["valgrind_binary"])
    library = Path(capabilities["VALGRIND_LIB"])
    if (
        not valgrind.is_file()
        or not library.is_dir()
        or any(c.isspace() for c in str(library))
    ):
        raise RuntimeError("Valgrind and its no-space VALGRIND_LIB path must exist")
    env = dict(os.environ)
    env.update(
        CONTROLLED_ENV, RAYON_NUM_THREADS=str(args.threads), VALGRIND_LIB=str(library)
    )
    env.pop("STERICX_PROFILE_PATH", None)
    version = subprocess.check_output(
        [str(valgrind), "--version"], env=env, text=True
    ).strip()
    target = root / "valgrind" / args.label
    target.mkdir(parents=True, exist_ok=False)
    summary = {
        "schema_version": 1,
        "created_utc": now(),
        "status": "running",
        "tool": args.tool,
        "tool_version": version,
        "valgrind_binary": str(valgrind),
        "valgrind_binary_sha256": digest(valgrind),
        "capabilities": capabilities,
        "capabilities_sha256": digest(capabilities_path),
        "binary": str(binary),
        "binary_sha256": digest(binary),
        "manifest": str(manifest_path),
        "manifest_sha256": manifest_hash,
        "input_fingerprints": manifest["files"],
        "script_sha256": digest(Path(__file__)),
        "oracle_helper_sha256": digest(REPO / "scripts/profile_stericx.py"),
        "baseline_label": args.compare_to,
        "baseline_summary_sha256": digest(baseline_path),
        "environment": {
            **CONTROLLED_ENV,
            "RAYON_NUM_THREADS": str(args.threads),
            "VALGRIND_LIB": str(library),
        },
        "affinity": sorted(affinity),
        "cache_configuration": configured_caches if args.tool == "cachegrind" else None,
        "sysfs_cache_geometry": sysfs_caches(affinity),
        "cache_configuration_policy": (
            "Explicit I1/D1/LL geometry avoids incorrect CPUID auto-detection. "
            "Defaults match this baseline host's cpu2 sysfs geometry; override "
            "for another host. Intermediate L2 is not modeled by Cachegrind. "
            "Matching sizes/ways does not make the basic simulation a native "
            "hardware-counter measurement."
        ),
        "per_process_timeout_seconds": args.timeout,
        "workloads": selected,
        "results": {},
        "simulation": SIMULATION if args.tool == "cachegrind" else None,
        "timing_policy": (
            "Profiler elapsed time retained only as run metadata; never use for "
            "native wall-time or Amdahl estimates."
        ),
    }
    save(target / "summary.json", summary)
    try:
        for name in selected:
            workload = manifest["workloads"][name]
            directory = target / name
            directory.mkdir()
            app_argv = [
                str(binary),
                *(
                    value.replace("{run_dir}", str(directory))
                    for value in workload["argv"]
                ),
            ]
            events = directory / (
                "events.out" if args.tool == "cachegrind" else "dhat.json"
            )
            argv = [
                str(valgrind),
                f"--tool={args.tool}",
                f"--log-file={directory / 'profiler.log'}",
            ]
            if args.tool == "cachegrind":
                argv.extend(
                    [
                        "--cache-sim=yes",
                        "--branch-sim=yes",
                        f"--I1={args.i1}",
                        f"--D1={args.d1}",
                        f"--LL={args.ll}",
                        f"--cachegrind-out-file={events}",
                    ]
                )
            else:
                argv.extend(
                    ["--mode=heap", "--num-callers=20", f"--dhat-out-file={events}"]
                )
            argv.extend(app_argv)
            state = {
                "status": "starting",
                "created_utc": now(),
                "argv": argv,
                "application_argv": app_argv,
            }
            save(directory / "run.json", state)
            summary["active_workload"] = name
            save(target / "summary.json", summary)

            def started(pid: int, state=state, directory=directory) -> None:
                state.update(status="running", pid=pid)
                save(directory / "run.json", state)

            print(
                f"{args.tool} {name}: starting unchanged manifest workload", flush=True
            )
            try:
                state.update(
                    run_child(
                        argv,
                        env,
                        affinity,
                        directory / "stdout.txt",
                        directory / "stderr.txt",
                        started,
                        args.timeout,
                    )
                )
                if state["returncode"] != 0:
                    raise RuntimeError(f"profiler/target failed: {directory}")
                loader_errors = (directory / "stderr.txt").read_text() + (
                    directory / "profiler.log"
                ).read_text()
                if (
                    "cannot be preloaded" in loader_errors
                    or "Fatal error" in loader_errors
                ):
                    raise RuntimeError(f"invalid profiler runtime: {directory}")
                fingerprints = check_against_native(
                    workload, directory, baseline["results"][name]["fingerprints"]
                )
                result = (
                    parse_cachegrind(events)
                    if args.tool == "cachegrind"
                    else parse_dhat(events)
                )
                if args.tool == "cachegrind":
                    annotate = [
                        capabilities["cg_annotate"],
                        "--no-annotate",
                        "--threshold=0.1",
                        "--show=" + ",".join(CACHE_EVENTS),
                        "--sort=Ir",
                        str(events),
                    ]
                    state["annotation_argv"] = annotate
                    annotation = run_child(
                        annotate,
                        env,
                        affinity,
                        directory / "functions.txt",
                        directory / "annotation.stderr.txt",
                        lambda pid: None,
                        args.timeout,
                    )
                    if annotation["returncode"] != 0:
                        raise RuntimeError(f"cg_annotate failed: {directory}")
                elif args.allocator_reference:
                    comparison = allocation_comparison(
                        root, args.allocator_reference, name, result
                    )
                    if comparison["reference_manifest_sha256"] != manifest_hash:
                        raise RuntimeError(
                            "allocator reference used different frozen inputs"
                        )
                    result["rust_allocator_crosscheck"] = comparison
                state.update(
                    status="complete",
                    fingerprints=fingerprints,
                    exact_scientific_output_match=True,
                )
                result.update(
                    fingerprints=fingerprints,
                    exact_scientific_output_match=True,
                    exactness_scope=workload.get(
                        "exactness_scope", "all emitted scientific values byte-for-byte"
                    ),
                    artifacts={
                        path.name: digest(path)
                        for path in directory.iterdir()
                        if path.is_file() and path.name != "run.json"
                    },
                )
                save(directory / "summary.json", result)
                summary["results"][name] = result
                print(
                    f"{args.tool} {name}: complete; "
                    "scientific output matches native oracle",
                    flush=True,
                )
            except BaseException as error:
                state.update(status="failed", error=f"{type(error).__name__}: {error}")
                raise
            finally:
                state["finished_utc"] = now()
                save(directory / "run.json", state)
                save(target / "summary.json", summary)
        verify_manifest(manifest)
        if digest(binary) != summary["binary_sha256"]:
            raise RuntimeError("profiled binary changed during run")
        summary.pop("active_workload", None)
        summary["status"] = "complete"
    except BaseException as error:
        summary.update(status="failed", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        summary["finished_utc"] = now()
        save(target / "summary.json", summary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--label", required=True)
    parser.add_argument(
        "--compare-to", required=True, help="completed native run label"
    )
    parser.add_argument(
        "--binary", type=Path, help="default frozen baseline symbol build"
    )
    parser.add_argument("--capabilities", type=Path)
    parser.add_argument("--tool", choices=("cachegrind", "dhat"), default="cachegrind")
    parser.add_argument(
        "--workloads",
        help=(
            "comma-separated unchanged manifest names; default all for Cachegrind, "
            "single_small/conformers_56 for DHAT"
        ),
    )
    parser.add_argument(
        "--allocator-reference",
        help="completed feature-instrumented run label for DHAT count cross-check",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--affinity", default="2")
    parser.add_argument("--i1", default="32768,8,64", help="I1 size,ways,line bytes")
    parser.add_argument("--d1", default="32768,8,64", help="D1 size,ways,line bytes")
    parser.add_argument(
        "--ll", default="16777216,16,64", help="LL size,ways,line bytes"
    )
    parser.add_argument(
        "--timeout",
        type=float,
        help=(
            "per-process timeout in seconds (default: none); expiration fails the "
            "run and stops its process group, never reduces the workload"
        ),
    )
    args = parser.parse_args()
    if args.timeout is not None and args.timeout <= 0:
        parser.error("--timeout must be positive")

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"received signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    run(args)


if __name__ == "__main__":
    main()
