#!/usr/bin/env python3
"""Replay the unchanged admitted 17-command model demonstration on a new build."""

from __future__ import annotations

import argparse
import importlib.util
import shutil
from pathlib import Path

import replay_scientific_remediation as replay

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / ".stericx/profiling/scientifically_validated_optimization"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.build = args.build.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    built, _ = replay.verify_build(args.build / "manifest.json")
    controller = FROZEN / "candidate_reference_controller_v1/controller.py"
    spec = importlib.util.spec_from_file_location(
        "frozen_candidate_reference", controller
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    frozen_path = FROZEN / "accurate_baseline_v1/manifest.json"
    frozen = replay.load_manifest(frozen_path)
    baseline = module.Baseline(frozen)
    journal = {
        "kind": "paton_fixed_model_demonstration",
        "commands": [],
        "helper": replay.oracle.file_record(Path(__file__)),
        "frozen_controller": replay.oracle.file_record(controller),
        "frozen_inventory": replay.oracle.file_record(frozen_path),
        "build": replay.oracle.file_record(args.build / "manifest.json"),
    }
    shutil.copyfile(controller, args.output / "frozen_controller.py")
    module.demo(args, built, baseline, journal)
    replay.verify_build(args.build / "manifest.json")
    module.verify_all(journal)
    replay.manifest(args.output / "receipt.json", journal)


if __name__ == "__main__":
    main()
