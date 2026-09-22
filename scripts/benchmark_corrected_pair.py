#!/usr/bin/env python3
"""Run admitted native A/B pairs; see benchmark_corrected_pair.md for the gate.

No candidate process, including a warmup, starts before every scientific gate
passes. Native timings include startup and output, using the existing wait4
launcher. Every attempted sample survives; no outliers are removed. Prediction
CLI output is only a preview/aggregate: its complete f32 stream is a separate
mandatory admission check, never claimed to be emitted by each timed process.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import statistics
import sys
from pathlib import Path

import compare_corrected_cli as cli_compare
import profile_stericx as profile
import scientifically_exact_optimization as exact

oracle = exact.oracle
replay = exact.replay
THREADS = (1, 2, 4, 6)
PAIRS = 8
WORKLOADS = {
    "single_small",
    "single_large",
    "ensemble_sdf",
    "conformers_56",
    "descriptors_10000",
    "parse_ensembles_1000",
    "predict_1000000",
    "search_database",
    "screen_1000",
    "db_build_ensembles",
}
ENGINEERING = {
    "rust-tests",
    "python-tests",
    "clippy",
    "rust-format",
    "rustdoc",
    "ruff-check",
    "ruff-format",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_records(value):
    """Rehash linked raw records, without treating their presence as success."""
    if isinstance(value, dict):
        if {"path", "bytes", "sha256"} <= value.keys():
            oracle.verify(value)
        else:
            for child in value.values():
                verify_records(child)
    elif isinstance(value, list):
        for child in value:
            verify_records(child)


def receipt(record, kind=None):
    path = oracle.verify(record)
    result = replay.load_manifest(path)
    if kind:
        require(result.get("kind") == kind, f"Expected {kind}: {path}")
    verify_records(result)
    return result


def bound_builds(value, baseline, candidate):
    require(value.get("baseline_build") == baseline, "Wrong baseline build")
    require(value.get("candidate_build") == candidate, "Wrong candidate build")


def contains_record(value, wanted):
    if isinstance(value, dict):
        return value == wanted or any(
            contains_record(child, wanted) for child in value.values()
        )
    if isinstance(value, list):
        return any(contains_record(child, wanted) for child in value)
    return False


def verify_result(record):
    """Bind legacy path->SHA inventories as well as modern file records."""
    path = oracle.verify(record)
    result = oracle.read(path)
    verify_records(result)
    for field in ("inputs", "artifacts"):
        values = result.get(field, {})
        if not isinstance(values, dict):
            continue
        for name, expected in values.items():
            if (
                isinstance(expected, str)
                and len(expected) == 64
                and all(c in "0123456789abcdef" for c in expected)
            ):
                source = Path(name)
                if not source.is_absolute():
                    source = path.parent / source
                require(
                    source.is_file() and oracle.sha(source) == expected,
                    f"Legacy reference raw evidence changed: {source}",
                )
    return result


def verify_observation_comparison(record, frozen, baseline, candidate):
    value = receipt(record, "exact_corrected_baseline_comparison")
    require(
        value.get("passed") is True and value.get("tolerance") == 0,
        "Observation comparison did not pass exactly",
    )
    require(
        value["executing_helper"]["sha256"] == oracle.sha(Path(exact.__file__)),
        "Unknown observation comparison helper",
    )
    require(
        value["baseline"] == frozen["admitted_observations"],
        "Observation baseline is not the admitted frozen baseline",
    )
    observed = [
        receipt(value[side], "scientific_remediation_observations")
        for side in ("baseline", "candidate")
    ]
    for item, build in zip(observed, (baseline, candidate), strict=True):
        require(
            item.get("complete_capture") is True and item["build_manifest"] == build,
            "Incomplete or wrong-build observations",
        )
        built = oracle.read(oracle.verify(build))
        require(
            item["binary"] == built["binaries"]["observer"], "Wrong observer binary"
        )
    old, new = observed
    require(
        len(old["lanes"]) == 13
        and set(old["lanes"]) == set(new["lanes"]) == set(value["lanes"]),
        "All thirteen observation lanes required",
    )
    for name, left in old["lanes"].items():
        right, checked = new["lanes"][name], value["lanes"][name]
        require(
            left["requests"]["sha256"] == right["requests"]["sha256"],
            f"Changed observation requests: {name}",
        )
        require(
            left["rows"] > 0
            and left["rows"] == right["rows"] == checked["rows"]
            and left["returncode"] == right["returncode"] == 0
            and checked["passed"] is True
            and checked["mismatches"] == 0
            and checked["stderr_equal"] is True
            and left["stderr"]["sha256"] == right["stderr"]["sha256"],
            f"Incomplete or failed observation lane: {name}",
        )
        require(
            left["stdout"]["sha256"] == right["stdout"]["sha256"],
            f"Observation stdout bytes differ despite comparison claim: {name}",
        )


def verify_cli_comparisons(records, baseline, candidate):
    require(
        set(records) == {str(t) for t in THREADS},
        "CLI comparisons need all four threads",
    )
    for thread, record in records.items():
        value = receipt(record, "corrected_cli_exact_comparison")
        require(
            value.get("complete") is True
            and value.get("exact_equivalence") is True
            and value.get("differences") == []
            and value.get("cases") == 93
            and value.get("mode") == "candidate"
            and value.get("baseline_threads")
            == value.get("candidate_threads")
            == int(thread),
            "Incomplete or failed candidate CLI comparison",
        )
        helpers = {Path(item["path"]).name: item for item in value["executing_helpers"]}
        require(
            helpers[Path(cli_compare.__file__).name]["sha256"]
            == oracle.sha(Path(cli_compare.__file__)),
            "Unknown CLI comparison helper",
        )
        captures = [
            receipt(value[side], "corrected_cli_capture")
            for side in ("baseline", "candidate")
        ]
        plan = cli_compare.verify_plan(Path(value["plan"]["path"]).parent)
        for capture, build in zip(captures, (baseline, candidate), strict=True):
            require(
                capture.get("complete_capture") is True
                and capture["build_manifest"] == build
                and capture["threads"] == int(thread)
                and capture["plan"] == value["plan"]
                and len(capture["cases"]) == 93,
                "CLI capture does not bind all cases to the requested build",
            )
            built = oracle.read(oracle.verify(build))
            require(
                capture["binary"] == built["binaries"]["native"], "Wrong CLI binary"
            )
        for side, capture in zip(("baseline", "candidate"), captures, strict=True):
            capture_root = Path(value[side]["path"]).parent
            require(
                set(capture["cases"]) == {case["id"] for case in plan["cases"]},
                "CLI capture membership differs from the canonical plan",
            )
            for case in plan["cases"]:
                directory = capture_root / case["id"]
                observed = cli_compare.cli.fingerprint(directory)
                require(
                    observed == capture["cases"][case["id"]]["fingerprint"],
                    "CLI stored fingerprint does not match raw output",
                )
                cli_compare.required_outputs(case, directory, observed["returncode"])
        require(captures[0]["cases"] == captures[1]["cases"], "CLI fingerprints differ")


def verify_reference_review(record, baseline, candidate):
    """A reviewed conclusion must point to actual passing reference result fields."""
    value = receipt(record, "corrected_candidate_independent_reference_gate")
    bound_builds(value, baseline, candidate)
    require(
        value.get("complete") is True and value.get("passed") is True,
        "Independent reference review is incomplete",
    )
    require(
        set(value["domains"]) == {"geometry", "thermodynamics", "models"},
        "All three independent reference domains are required",
    )
    for name, domain in value["domains"].items():
        require(
            domain.get("passed") is True
            and domain.get("checks")
            and domain.get("evidence"),
            f"Reference evidence/checks missing: {name}",
        )
        for evidence in domain["evidence"]:
            # Reference receipts frequently inventory their own raw observations
            # and unchanged equation programs. Bind those bytes as well.
            verify_result(evidence)
        for check in domain["checks"]:
            require(
                check["evidence"] in domain["evidence"],
                "Reference check is not inventoried",
            )
            current = verify_result(check["evidence"])
            require(
                contains_record(current, candidate),
                f"Checked reference result lacks candidate build provenance: {name}",
            )
            require(check.get("field"), "Reference check field is empty")
            for key in check["field"]:
                current = current[key]
            # Success, zero failures, or an empty failure list; no arbitrary
            # string assertion can stand in for a scientific review result.
            expected = check["equals"]
            require(
                expected is True
                or (type(expected) is int and expected == 0)
                or expected == [],
                "Reference check must assert success or zero failures",
            )
            require(
                type(current) is type(expected) and current == expected,
                f"Independent reference check failed: {name}/{check['field']}",
            )


def verify_engineering(record, baseline, candidate, built):
    value = receipt(record, "corrected_candidate_engineering_gate")
    bound_builds(value, baseline, candidate)
    require(
        value.get("complete") is True and value.get("passed") is True,
        "Engineering gate is incomplete",
    )
    require(
        value["sources"] == built["sources"] + built["python_sources"],
        "Engineering gate does not bind the complete candidate source snapshot",
    )
    commands = value["commands"]
    require(
        ENGINEERING <= {item["name"] for item in commands}, "Missing engineering checks"
    )
    require(
        len({item["name"] for item in commands}) == len(commands),
        "Duplicate engineering check",
    )
    for item in commands:
        require(
            item.get("returncode") == 0 and item.get("argv") and item.get("log"),
            "Engineering command failed or lacks raw evidence",
        )


def verify_prediction_exporters(value, baseline, candidate):
    require(
        set(value.get("exporter_builds", {})) == {"baseline", "candidate"},
        "Both prediction exporter build receipts are required",
    )
    exporter_sources = []
    for side, build_record in (("baseline", baseline), ("candidate", candidate)):
        proof = receipt(
            value["exporter_builds"][side], "corrected_frozen_diagnostic_builds"
        )
        build = oracle.read(oracle.verify(build_record))
        require(
            proof.get("complete") is True
            and proof.get("native_identical_to_admitted_executable") is True,
            "Prediction exporter build is incomplete",
        )
        native = proof["binaries"]["native"]
        require(
            native["stericx"]["sha256"]
            == build["binaries"]["native"]["sha256"]
            == proof["native_timing_executable"]["sha256"],
            "Prediction exporter is not linked to the admitted native build",
        )
        commands = [
            command for command in proof["commands"] if command["kind"] == "native"
        ]
        require(
            len(commands) == 1 and commands[0]["returncode"] == 0,
            "Native exporter compilation failed",
        )
        command = commands[0]
        require(
            "--release" in command["argv"]
            and "--locked" in command["argv"]
            and "--features" not in command["argv"]
            and "--all-features" not in command["argv"]
            and "profile_predictions" in command["argv"],
            "Wrong exporter compilation settings",
        )
        source = Path(command["cwd"])
        for item in build["sources"]:
            relative = Path(item["path"]).relative_to(build["source"])
            copied = oracle.file_record(source / relative)
            require(
                copied["sha256"] == item["sha256"] and copied in proof["bindings"],
                "Prediction exporter compiled different or unbound source",
            )
        example = oracle.file_record(source / "examples/profile_predictions.rs")
        require(example in proof["bindings"], "Prediction export adapter is not bound")
        exporter_sources.append(example["sha256"])
        exports = [item for item in proof["exports"] if item["kind"] == "native"]
        require(
            len(exports) == 4 and {item["threads"] for item in exports} == set(THREADS),
            "All four native prediction captures are required",
        )
        for capture in exports:
            thread = str(capture["threads"])
            stream = value["threads"][thread][side]
            require(
                capture["returncode"] == 0
                and capture["vector"] == stream
                and capture["affinity"] == exact.inventory.AFFINITIES[int(thread)]
                and capture["argv"][-4:]
                == [
                    native["profile_predictions"]["path"],
                    value["data"]["path"],
                    value["weights"]["path"],
                    stream["path"],
                ],
                "Prediction capture does not bind its stream to executable/settings",
            )
    require(len(set(exporter_sources)) == 1, "Prediction exporter adapter changed")


def verify_prediction(record, baseline, candidate, workloads):
    require(
        sys.byteorder == "little",
        "The frozen prediction stream contract requires little-endian host bytes",
    )
    value = receipt(record, "corrected_prediction_bitstream_comparison")
    bound_builds(value, baseline, candidate)
    require(
        value.get("complete") is True
        and value.get("exact_equivalence") is True
        and value.get("encoding") == "ieee754-f32-le"
        and value.get("records") == 1_000_000,
        "Complete million-prediction f32 streams are required",
    )
    argv = workloads["workloads"]["predict_1000000"]["argv"]
    for field, flag in (("data", "--data"), ("weights", "--weights")):
        require(
            value[field] == oracle.file_record(Path(argv[argv.index(flag) + 1])),
            "Prediction gate uses different frozen inputs",
        )
    require(
        set(value["threads"]) == {str(t) for t in THREADS},
        "Prediction gate needs all four threads",
    )
    for streams in value["threads"].values():
        require(
            set(streams) == {"baseline", "candidate"},
            "Both prediction streams required",
        )
        left, right = streams["baseline"], streams["candidate"]
        require(
            left["bytes"] == right["bytes"] == 4_000_000
            and left["sha256"] == right["sha256"],
            "Full prediction bitstreams differ or are truncated",
        )
    verify_prediction_exporters(value, baseline, candidate)


def full_fingerprint(workload, directory):
    return {
        **profile.check_result(workload, directory),
        "stderr_sha256": profile.digest(directory / "stderr.txt"),
    }


def verify_baseline_measurements(records, frozen, workloads, binary):
    require(
        set(records) == {str(t) for t in THREADS},
        "Native baselines need all four threads",
    )
    fingerprints = {}
    for thread, record in records.items():
        path = oracle.verify(record)
        value = oracle.read(path)
        require(
            value.get("status") == "complete"
            and value["binary_sha256"] == binary["sha256"]
            and value["manifest_sha256"] == frozen["workload_manifest"]["sha256"]
            and value["environment"]["RAYON_NUM_THREADS"] == thread
            and value["affinity"]
            == [int(cpu) for cpu in frozen["thread_affinities"][thread].split(",")]
            and value["warmups"] == 1
            and value["repetitions"] >= 7
            and set(value["results"]) == WORKLOADS,
            "Wrong or incomplete native baseline measurements",
        )
        fingerprints[thread] = {}
        for name, result in value["results"].items():
            require(
                result["measured_runs"] == value["repetitions"],
                "Incomplete baseline sample count",
            )
            observed = None
            for phase, count in (("warmup", 1), ("measured", value["repetitions"])):
                for index in range(count):
                    directory = path.parent / name / f"{phase}_{index:02d}"
                    current = full_fingerprint(workloads["workloads"][name], directory)
                    metrics = oracle.read(directory / "metrics.json")
                    require(
                        metrics["returncode"] == 0
                        and metrics["phase"] == phase
                        and metrics["fingerprints"] == result["fingerprints"]
                        and metrics["raw_stdout_sha256"]
                        == profile.digest(directory / "stdout.txt")
                        and metrics["raw_stderr_sha256"] == current["stderr_sha256"]
                        and {k: v for k, v in current.items() if k != "stderr_sha256"}
                        == result["fingerprints"],
                        "Baseline raw output differs from its receipt",
                    )
                    if observed is not None:
                        require(
                            current == observed,
                            "Baseline output/stderr is nondeterministic",
                        )
                    observed = current
            fingerprints[thread][name] = observed
    return fingerprints


def verify_gate(snapshot, candidate_build, gate_path):
    # Called only at phase boundaries, never for individual process launches.
    frozen = exact.verify_frozen(snapshot / "manifest.json")
    gate_identity = oracle.file_record(gate_path)
    gate = receipt(gate_identity, "corrected_candidate_performance_gate")
    require(gate.get("complete") is True, "Scientific gate is incomplete")
    require(
        gate["baseline_snapshot"] == oracle.file_record(snapshot / "manifest.json"),
        "Wrong baseline snapshot",
    )
    copies = {
        Path(item["snapshot"]["path"]).relative_to(snapshot).as_posix(): item
        for item in frozen["copies"]
    }
    baseline = copies["receipts/build.json"]["original"]
    candidate = oracle.file_record(candidate_build)
    bound_builds(gate, baseline, candidate)
    built, _ = replay.verify_build(candidate_build)
    baseline_binary = copies["bin/native"]["snapshot"]
    require(
        baseline_binary["path"] == str(snapshot / "bin/native"),
        "Wrong frozen native path",
    )
    workloads = oracle.read(oracle.verify(frozen["workload_manifest"]))
    require(
        workloads.get("kind") == "corrected_science_profile_inputs"
        and workloads["build_manifest"] == baseline
        and set(workloads["workloads"]) == WORKLOADS,
        "Expected all ten corrected frozen workloads",
    )
    require(
        frozen["thread_affinities"]
        == {str(t): exact.inventory.AFFINITIES[t] for t in THREADS},
        "Unexpected frozen thread affinities",
    )
    verify_observation_comparison(
        gate["observation_comparison"], frozen, baseline, candidate
    )
    verify_cli_comparisons(gate["cli_comparisons"], baseline, candidate)
    verify_reference_review(gate["independent_references"], baseline, candidate)
    verify_engineering(gate["engineering"], baseline, candidate, built)
    verify_prediction(gate["prediction_comparison"], baseline, candidate, workloads)
    expected = verify_baseline_measurements(
        gate["baseline_measurements"], frozen, workloads, baseline_binary
    )
    oracle.verify(gate_identity)
    return dict(
        frozen=frozen,
        workloads=workloads,
        expected=expected,
        gate=gate_identity,
        builds={"baseline": baseline, "candidate": candidate},
        binaries={
            "baseline": baseline_binary,
            "candidate": built["binaries"]["native"],
        },
    )


def schedule():
    for side in ("baseline", "candidate"):
        yield "warmup", 0, side
    for pair in range(PAIRS):
        order = (
            ("baseline", "candidate") if pair % 2 == 0 else ("candidate", "baseline")
        )
        for side in order:
            yield "measured", pair, side


def settings(thread, affinities):
    require(thread in THREADS, "Unsupported thread count")
    affinity = set(map(int, affinities[str(thread)].split(",")))
    require(
        affinity and affinity <= os.sched_getaffinity(0), "Requested CPUs unavailable"
    )
    env = dict(os.environ)
    env.update(profile.CONTROLLED_ENV, RAYON_NUM_THREADS=str(thread))
    env.pop("STERICX_PROFILE_PATH", None)
    return env, affinity


def distribution(values):
    require(
        values and all(math.isfinite(value) for value in values),
        "Nonfinite sample distribution",
    )
    median = statistics.median(values)
    return dict(
        values=values,
        median=median,
        mad=statistics.median(abs(v - median) for v in values),
        minimum=min(values),
        maximum=max(values),
    )


def paired_summary(samples):
    measured = [row for row in samples if row["phase"] == "measured"]
    by_side = {
        side: [r for r in measured if r["side"] == side]
        for side in ("baseline", "candidate")
    }
    require(
        len(measured) == PAIRS * 2
        and all(len(rows) == PAIRS for rows in by_side.values()),
        "Incomplete measured pairs cannot be summarized as complete",
    )
    pairs = []
    for index in range(PAIRS):
        pair = {row["side"]: row for row in measured if row["pair"] == index}
        require(
            set(pair) == {"baseline", "candidate"}, "Missing or duplicate pair member"
        )
        a, b = pair["baseline"], pair["candidate"]
        require(a["wall_ns"] > 0 and b["wall_ns"] > 0, "Invalid native duration")
        pairs.append(a["wall_ns"] / b["wall_ns"])
    summaries = {
        side: {
            **profile.summarize(rows),
            "cpu_total_s": distribution(
                [row["cpu_user_s"] + row["cpu_system_s"] for row in rows]
            ),
        }
        for side, rows in by_side.items()
    }
    return dict(
        **summaries,
        paired_wall_speedup=distribution(pairs),
        speedup_definition="baseline wall_ns / candidate wall_ns within the same pair",
        pairs=PAIRS,
        exclusions=[],
    )


def measure_configuration(name, workload, thread, context, output, launcher, journal):
    env, affinity = settings(thread, context["frozen"]["thread_affinities"])
    expected = context["expected"][str(thread)][name]
    samples = []
    for phase, pair, side in schedule():
        directory = output / f"t{thread}" / name / f"{phase}_{pair:02d}" / side
        argv = [
            context["binaries"][side]["path"],
            *(arg.replace("{run_dir}", str(directory)) for arg in workload["argv"]),
        ]
        attempt = dict(
            workload=name,
            threads=thread,
            phase=phase,
            pair=pair,
            side=side,
            directory=str(directory),
            argv=argv,
        )
        journal["attempts"].append(attempt)
        profile.write_json(output / "summary.json", journal)
        try:
            metrics = profile.measure(argv, env, affinity, directory, launcher)
            fingerprint = full_fingerprint(workload, directory)
            metrics.update(phase=phase, pair=pair, side=side, fingerprints=fingerprint)
            profile.write_json(directory / "metrics.json", metrics)
            # Every warmup and measured sample is checked against the admitted
            # baseline measurement, including exact successful stderr bytes.
            require(fingerprint == expected, f"Exact output mismatch: {directory}")
            attempt.update(
                status="complete",
                metrics=oracle.file_record(directory / "metrics.json"),
            )
            samples.append(metrics)
        except BaseException as error:
            attempt.update(status="failed", error=str(error))
            raise
        finally:
            profile.write_json(output / "summary.json", journal)
        print(
            f"t{thread} {name} {phase} {pair + 1} {side}: "
            f"{metrics['wall_ns'] / 1e6:.3f} ms",
            flush=True,
        )
    return dict(
        summary=paired_summary(samples),
        fingerprints=expected,
        exactness_scope=workload.get(
            "exactness_scope", "all emitted scientific values byte-for-byte"
        ),
    )


def executing_helpers():
    modules = (
        sys.modules[__name__],
        profile,
        exact,
        exact.inventory,
        replay,
        oracle,
        replay.refs,
        cli_compare,
        cli_compare.capture,
        cli_compare.cli,
    )
    return sorted(
        {Path(module.__file__).resolve() for module in modules}
        | {Path(profile.__file__).with_name("profile_child.c")}
    )


def run(args):
    require(sys.platform == "linux", "Native wait4 measurements require Linux")
    snapshot, build, gate = (
        args.snapshot.resolve(),
        args.candidate_build.resolve(),
        args.gate.resolve(),
    )
    destination = args.output.resolve()
    require(
        not destination.is_relative_to(snapshot)
        and not destination.is_relative_to(build.parent),
        "Benchmark output must be outside frozen snapshot and candidate build",
    )
    output = replay.new_output(args.output)
    helpers = [oracle.file_record(path) for path in executing_helpers()]
    journal = dict(
        kind="corrected_native_paired_benchmark",
        status="admitting",
        created_utc=profile.now(),
        helpers=helpers,
        attempts=[],
        results={},
        pairs_per_configuration=PAIRS,
        warmups_per_binary=1,
        threads=list(THREADS),
        exclusions=[],
    )
    profile.write_json(output / "summary.json", journal)
    context = None
    try:
        context = verify_gate(snapshot, build, gate)
        journal.update(
            status="running",
            gate=context["gate"],
            builds=context["builds"],
            binaries=context["binaries"],
            workload_manifest=context["frozen"]["workload_manifest"],
            baseline_snapshot=oracle.file_record(snapshot / "manifest.json"),
            cache_policy=context["workloads"]["cache_policy"],
            environment={
                key: value
                for key, value in os.environ.items()
                if key in exact.inventory.ENV_KEYS
            },
            controlled_environment=profile.CONTROLLED_ENV,
            thread_affinities=context["frozen"]["thread_affinities"],
            expected_launches=len(WORKLOADS) * len(THREADS) * (2 + 2 * PAIRS),
            preflight_verified_utc=profile.now(),
        )
        for item in [*helpers, context["gate"]]:
            source = oracle.verify(item)
            target = output / "helpers" / source.name
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(source, target)
            require(
                oracle.sha(target) == item["sha256"],
                "Archived helper/gate copy changed",
            )
        journal["archived_helpers"] = [
            oracle.file_record(path) for path in sorted((output / "helpers").iterdir())
        ]
        launcher = profile.build_launcher(output / "launcher")
        launcher_records = [
            oracle.file_record(path)
            for path in sorted(launcher.parent.iterdir())
            if path.is_file()
        ]
        journal["launcher"] = launcher_records
        # Validate every affinity before any process starts.
        for thread in THREADS:
            settings(thread, context["frozen"]["thread_affinities"])
        for thread in THREADS:
            for name, workload in context["workloads"]["workloads"].items():
                journal["results"][f"t{thread}/{name}"] = measure_configuration(
                    name, workload, thread, context, output, launcher, journal
                )
        verify_records(launcher_records)
        journal["status"] = "complete"
    except BaseException as error:
        journal.update(status="failed", error=str(error))
        raise
    finally:
        try:
            verify_records(helpers)
            verify_records(journal.get("archived_helpers", []))
            verify_records(journal.get("launcher", []))
            if context is not None:
                oracle.verify(context["gate"])
                after = verify_gate(snapshot, build, gate)
                require(
                    after == context,
                    "Scientific inputs or admission changed during measurements",
                )
                for attempt in journal["attempts"]:
                    if attempt.get("status") != "complete":
                        continue
                    metrics = oracle.read(oracle.verify(attempt["metrics"]))
                    directory = Path(attempt["directory"])
                    require(
                        full_fingerprint(
                            context["workloads"]["workloads"][attempt["workload"]],
                            directory,
                        )
                        == context["expected"][str(attempt["threads"])][
                            attempt["workload"]
                        ]
                        and profile.digest(directory / "stdout.txt")
                        == metrics["raw_stdout_sha256"],
                        "Retained sample output changed after measurement",
                    )
                journal["postflight_verified_utc"] = profile.now()
        except BaseException as error:
            journal.update(status="failed", postflight_error=str(error))
        journal["finished_utc"] = profile.now()
        profile.write_json(output / "summary.json", journal)
        # Inventory every attempted sample, even if admission or a process failed.
        replay.manifest(
            output / "manifest.json",
            dict(
                kind="corrected_native_paired_benchmark_inventory",
                complete=journal["status"] == "complete",
                files=[
                    oracle.file_record(path)
                    for path in sorted(output.rglob("*"))
                    if path.is_file()
                ],
            ),
        )
    require(
        journal["status"] == "complete",
        journal.get("postflight_error", "Benchmark incomplete"),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "candidate-build", "gate", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
