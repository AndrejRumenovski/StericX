#!/usr/bin/env python3
"""Measure same-input native, JSON-bridge and typed-Python ensemble calls."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

import stericx as sx


def distribution(values):
    median = statistics.median(values)
    return {
        "values_ns_per_call": values,
        "median_ns": median,
        "mad_ns": statistics.median(abs(x - median) for x in values),
        "minimum_ns": min(values),
        "maximum_ns": max(values),
    }


def record(path):
    path = Path(path).resolve()
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--options", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--iterations", default=100, type=int)
    args = parser.parse_args()
    if args.iterations < 1:
        raise ValueError("iterations must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, {2})
    os.environ.update(RAYON_NUM_THREADS="1", OMP_NUM_THREADS="1")
    ensemble = sx.Ensemble.from_file(args.input)
    values = json.loads(args.options.read_text())
    values["volume"] = sx.VolumeOptions(**values.get("volume", {}))
    options = sx.DescriptorOptions(**values)
    options_json = json.dumps(asdict(options))

    def bridge():
        return ensemble._engine.analyze(
            options_json, '{"missing":"require_complete"}', "null", None
        )

    def public():
        return ensemble.descriptors(options)

    def native(iterations):
        command = [
            str(args.native.resolve()),
            str(args.input),
            str(args.options),
            str(iterations),
        ]
        raw = subprocess.check_output(command)
        return raw, json.loads(raw)

    raw, checked = native(1)
    (args.output / "native_preflight.json").write_bytes(raw)
    expected = json.loads(bridge())
    if expected != checked["result"]:
        raise ValueError("Native and Python-bridge scientific results differ")
    typed = public()
    if typed.to_dict() != expected:
        raise ValueError("Typed public result differs from the native transport")
    (args.output / "bridge_preflight.json").write_text(json.dumps(expected))
    for method in [bridge, public]:
        method()
    samples = {key: [] for key in ["native", "bridge", "public"]}
    attempts = []
    result = {
        "command": sys.argv,
        "python": sys.version,
        "platform": platform.platform(),
        "stericx": sx.__version__,
        "bindings": [
            record(path)
            for path in [
                __file__,
                args.input,
                args.options,
                args.native,
                sx.__file__,
                sx._native.__file__,
                Path(sx.__file__).with_name("_types.py"),
            ]
        ],
        "iterations_per_block": args.iterations,
        "affinity": [2],
        "warmups": 1,
        "pairs": 8,
        "scientific_results_exact": True,
        "attempts": attempts,
        "complete": False,
    }
    receipt = args.output / "summary.json"
    for pair in range(8):
        order = ["native", "bridge", "public"]
        if pair % 2:
            order.reverse()
        for name in order:
            attempt = {"pair": pair, "method": name, "status": "started"}
            attempts.append(attempt)
            receipt.write_text(json.dumps(result, indent=2) + "\n")
            if name == "native":
                raw, observed = native(args.iterations)
                (args.output / f"native_{pair}.json").write_bytes(raw)
                if observed["result"] != expected:
                    raise ValueError("Native timed output changed")
                elapsed = observed["elapsed_ns"]
            else:
                method = bridge if name == "bridge" else public
                start = time.perf_counter_ns()
                for _ in range(args.iterations):
                    value = method()
                    del value
                elapsed = time.perf_counter_ns() - start
            samples[name].append(elapsed / args.iterations)
            attempt.update(status="complete", elapsed_ns=elapsed)
            receipt.write_text(json.dumps(result, indent=2) + "\n")
    result.update(
        complete=True,
        results={name: distribution(values) for name, values in samples.items()},
        interpretation=(
            "Native calculates descriptors and statistics with parsed input. "
            "Bridge adds options decoding, GIL transition and native JSON encoding. "
            "Public adds Python options encoding, JSON parsing and typed objects. "
            "Native executable and portable wheel have different build targets; "
            "their difference is not a pure FFI cost. Public minus bridge isolates "
            "Python conversion in the same installed wheel."
        ),
    )
    for item in result["bindings"]:
        if record(item["path"]) != item:
            raise ValueError("Bound input or executable changed during timing")
    receipt.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
