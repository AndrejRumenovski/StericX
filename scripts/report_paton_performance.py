#!/usr/bin/env python3
"""Recheck retained Paton timings and produce the confirmed-window report tables.

Reads existing evidence only; it never launches a benchmark or drops a sample.
Run from the repository root. Original /tmp paths are resolved to retained copies.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import benchmark_corrected_pair as paired
import benchmark_paton_extensions as benchmark
import profile_stericx as profile

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / ".stericx/paton"
CAMPAIGNS = (
    "performance_v1",
    "performance_v2_tmpfs",
    "performance_idle_v1",
    "performance_idle_tmpfs_v1",
    "performance_confirmed_20261102",
)


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_at(item, path):
    current = benchmark.record(path)
    require(
        all(current[key] == item[key] for key in ("sha256", "bytes")),
        f"Evidence changed: {path}",
    )


def review_campaign(name, workloads):
    directory = EVIDENCE / name
    summary = read(directory / "summary.json")
    require(summary["complete"], f"Incomplete campaign: {name}")
    require(len(summary["attempts"]) == 504, f"Missing samples: {name}")
    for item in summary["bindings"] + list(summary["binaries"].values()):
        benchmark.verify(item)
    original = Path(summary["attempts"][0]["directory"]).parents[3]
    groups = {}
    for attempt in summary["attempts"]:
        require(attempt["status"] == "complete", f"Failed sample: {attempt}")
        retained = directory / Path(attempt["directory"]).relative_to(original)
        verify_at(attempt["metrics"], retained / "metrics.json")
        metrics = read(retained / "metrics.json")
        fingerprint = paired.full_fingerprint(workloads[attempt["workload"]], retained)
        require(metrics["fingerprints"] == fingerprint, "Output fingerprint changed")
        key = (str(attempt["threads"]), attempt["workload"])
        groups.setdefault(key, []).append(metrics)
    rows, triggers = [], []
    require(len(groups) == 28, f"Missing configurations: {name}")
    for (thread, workload), samples in sorted(groups.items()):
        require(len(samples) == 18, "Warmup/measured coverage changed")
        expected = summary["results"][thread][workload]
        require(
            all(s["fingerprints"] == expected["fingerprints"] for s in samples),
            "Baseline/candidate scientific outputs differ",
        )
        actual = paired.paired_summary(samples)
        actual["fingerprints"] = expected["fingerprints"]
        for side in ("baseline", "candidate"):
            actual[side]["records_per_second"] = (
                workloads[workload]["expected_records"]
                * 1e9
                / actual[side]["wall_ns"]["median"]
            )
        require(actual == expected, f"Summary changed: {name}/{thread}/{workload}")
        wall_ratio = (
            actual["candidate"]["wall_ns"]["median"]
            / actual["baseline"]["wall_ns"]["median"]
        )
        paired_ratio = actual["paired_wall_speedup"]["median"]
        row = dict(
            threads=int(thread),
            workload=workload,
            candidate_over_baseline_medians=wall_ratio,
            paired_baseline_over_candidate=paired_ratio,
            metrics=actual,
        )
        rows.append(row)
        if wall_ratio > 1.1 or paired_ratio < 1 / 1.1:
            triggers.append(
                {key: value for key, value in row.items() if key != "metrics"}
            )
    return dict(
        campaign=name,
        source=benchmark.record(directory / "summary.json"),
        attempts=len(summary["attempts"]),
        confirmed_quiet_window=name == CAMPAIGNS[-1],
        all_raw_scientific_fingerprints_exact=True,
        all_summaries_recomputed_exactly=True,
        regression_triggers=triggers,
        maximum_median_slowdown_percent=100
        * (max(r["candidate_over_baseline_medians"] for r in rows) - 1),
        rows=rows,
    )


def timing_table(review):
    lines = [
        "# Confirmed-window native performance, 2026-11-02",
        "",
        "Eight AB/BA pairs per configuration; one warmup per binary. No exclusions.",
        "Times include process startup and regular-file output on the original disk.",
        "MAD is median absolute deviation. Speed ratio is the median of within-pair",
        "baseline/final wall times; it need not equal the ratio of side medians.",
        "",
        "| Threads | Workload | Baseline median ± MAD (ms) | Final median ± MAD (ms) | "
        "Baseline range (ms) | Final range (ms) | Final/baseline medians | "
        "Paired speed ratio |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in review["rows"]:
        b, c = (row["metrics"][side]["wall_ns"] for side in ("baseline", "candidate"))
        lines.append(
            f"| {row['threads']} | {row['workload']} | "
            f"{b['median'] / 1e6:.3f} ± {b['mad'] / 1e6:.3f} | "
            f"{c['median'] / 1e6:.3f} ± {c['mad'] / 1e6:.3f} | "
            f"{b['minimum'] / 1e6:.3f}-{b['maximum'] / 1e6:.3f} | "
            f"{c['minimum'] / 1e6:.3f}-{c['maximum'] / 1e6:.3f} | "
            f"{row['candidate_over_baseline_medians']:.5f} | "
            f"{row['paired_baseline_over_candidate']:.5f} |"
        )
    lines += [
        "",
        "## CPU, memory and throughput",
        "",
        "RSS is per-child peak resident memory; entries are medians in MiB.",
        "CPU is summed user + system time across threads. Throughput uses completed",
        "output records: one for the 86-conformer SDF and 11 for database builds.",
        "",
        "| Threads | Workload | Baseline / final CPU (ms) | "
        "Baseline / final RSS (MiB) | Baseline / final records/s |",
        "|---:|---|---:|---:|---:|",
    ]
    for row in review["rows"]:
        b, c = (row["metrics"][side] for side in ("baseline", "candidate"))
        lines.append(
            f"| {row['threads']} | {row['workload']} | "
            f"{b['cpu_total_s']['median'] * 1000:.3f} / "
            f"{c['cpu_total_s']['median'] * 1000:.3f} | "
            f"{b['peak_rss_bytes']['median'] / 2**20:.3f} / "
            f"{c['peak_rss_bytes']['median'] / 2**20:.3f} | "
            f"{b['records_per_second']:.2f} / {c['records_per_second']:.2f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    confirmed = read(EVIDENCE / CAMPAIGNS[-1] / "summary.json")
    workload_path = next(
        Path(item["path"])
        for item in confirmed["bindings"]
        if item["path"].endswith("workloads/manifest.json")
    )
    workloads = read(workload_path)["workloads"]
    campaigns = []
    for name in CAMPAIGNS:
        campaigns.append(review_campaign(name, workloads))
        print(f"Rechecked all 504 raw samples and summaries: {name}", flush=True)
    review = dict(
        kind="confirmed_paton_performance_review",
        created_utc=profile.now(),
        controller=benchmark.record(__file__),
        total_retained_attempts=sum(r["attempts"] for r in campaigns),
        campaigns=campaigns,
        confirmed_window_performance_preserved=not campaigns[-1]["regression_triggers"],
        scientific_global_pass=False,
        release_admitted=False,
        interpretation=(
            "No >10% median regression in the confirmed quiet window. "
            "All earlier campaigns and their triggers remain separate; "
            "no sample exclusions or production numerical changes. "
            "This finite performance result does not admit a release or "
            "override the known failed strict DBSTEP volume gate."
        ),
    )
    profile.write_json(args.output / "confirmed_review.json", review)
    (args.output / "CONFIRMED_TIMING_TABLE.md").write_text(timing_table(campaigns[-1]))


if __name__ == "__main__":
    main()
