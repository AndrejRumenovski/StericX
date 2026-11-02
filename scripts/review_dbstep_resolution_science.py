#!/usr/bin/env python3
"""Review fresh legacy replay with unchanged admitted comparison equations.

This checks legacy scientific preservation only. The new DBSTEP discrepancy
analysis and release decision require separate reviewed evidence.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import replay_scientific_remediation as replay

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / ".stericx/profiling/scientifically_validated_optimization"


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("build", "science", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    args.build, args.science, args.output = (
        p.resolve() for p in (args.build, args.science, args.output)
    )
    args.output.mkdir(parents=True, exist_ok=False)
    controller = FROZEN / "candidate_reference_controller_v1/controller.py"
    spec = importlib.util.spec_from_file_location(
        "retained_reference_review", controller
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    frozen = replay.load_manifest(FROZEN / "accurate_baseline_v1/manifest.json")
    baseline = module.Baseline(frozen)
    built, _ = replay.verify_build(args.build / "manifest.json")
    observed, observed_build = replay.verify_observations(
        args.science / "observations/manifest.json"
    )
    require(observed_build == built, "Wrong observation build")
    require(
        len(observed["lanes"]) == 13
        and sum(lane["rows"] for lane in observed["lanes"].values()) == 128021,
        "Incomplete legacy observations",
    )
    evidence = [replay.oracle.file_record(controller)]
    changes, comparisons = [], []
    for label, components in [
        ("inferred", ["geometry", "kraken"]),
        ("topology", ["kraken"]),
    ]:
        old = module.BASE / f"references_v4_{label}"
        new = args.science / f"references_{label}"
        baseline.check(old / "result.json")
        receipt = replay.load_manifest(new / "result.json")
        module.verify_all(receipt)
        started = read(Path(receipt["started"]["path"]))
        module.verify_all(started)
        require(
            started["observations"]
            == replay.oracle.file_record(args.science / "observations/manifest.json")
            and started["components"] == components
            and started["kraken_observation_lane"]
            == ("kraken" if label == "inferred" else "kraken_topology"),
            "Reference replay is not bound to the requested fresh observations",
        )
        require(
            receipt["commands_completed"]
            and all(
                r["returncode"] == 0 and not r["timed_out"] for r in receipt["runs"]
            ),
            "Reference commands incomplete",
        )

        class SealedBaseline:
            audit = old

            def check(self, relative):
                return baseline.check(self.audit / relative)

        destination = args.output / label
        destination.mkdir()
        result = replay.refs.compare_artifacts(
            SealedBaseline(), new, components, destination
        )
        require(
            {p.relative_to(old) for p in replay.refs.scientific_files(old, components)}
            == {
                p.relative_to(new)
                for p in replay.refs.scientific_files(new, components)
            },
            "Reference artifact coverage changed",
        )
        replay.manifest(destination / "comparison.json", result)
        comparisons.append(result)
        evidence += [
            replay.oracle.file_record(p)
            for p in (new / "result.json", destination / "comparison.json")
        ]

    volume = read(args.science / "volume_reference/complete.json")
    require(
        volume["complete"]
        and volume["limit"] is None
        and volume["requests"] == volume["compared"] == volume["source_rows"] == 31721
        and volume["reference_errors"] == 0
        and volume["complete_source_stream_validated"],
        "Incomplete full volume reference",
    )
    bins = read(args.science / "all_bins/complete.json")
    require(
        (bins["rows"], bins["frames"], bins["bins"]) == (31721, 95163, 1141956),
        "Incomplete bins",
    )
    union = read(args.science / "fixed_union/complete.json")
    old_union = read(
        baseline.check(module.BASE / "fixed_reference_union_v4/complete.json")
    )
    require(
        union["complete"]
        and union["union_count"] == 2517
        and union["memberships"] == old_union["memberships"],
        "Fixed union changed",
    )
    pairs = [
        (
            module.BASE / "buried_volume_reference_v4",
            args.science / "volume_reference",
            ["reference_results.jsonl"],
        ),
        (
            module.GEOMETRY / "all_bins_reference_v4",
            args.science / "all_bins",
            ["per_id_residuals.jsonl", "metrics.json", "top20.json"],
        ),
        (
            module.BASE / "thermodynamics/final_v4",
            args.science / "thermodynamics",
            [
                "summary.json",
                "kinetics_comparison.json",
                "aggregation_comparison.json",
                "python.json",
            ],
        ),
    ]
    for lane in ("inferred", "topology"):
        require(
            union["lanes"][lane]["rows"] == 2517
            and union["lanes"][lane]["coverage"]
            == old_union["lanes"][lane]["coverage"],
            "Fixed variant coverage changed",
        )
        pairs.append(
            (
                module.BASE / "fixed_reference_union_v4" / lane,
                args.science / "fixed_union" / lane,
                [
                    "selected_inputs.jsonl",
                    "reference_raw.jsonl",
                    "mean_recombination.jsonl",
                ],
            )
        )
    for old, new, names in pairs:
        for name in names:
            evidence.append(
                module.compare_file(old / name, new / name, baseline, changes, old, new)
            )

    thermo = args.science / "thermodynamics"
    for name, key in (("stericx", "native"), ("observer", "observer")):
        require(
            replay.oracle.sha(thermo / name) == built["binaries"][key]["sha256"],
            "Wrong thermodynamic binary",
        )
    sources = {
        Path(r["path"]).relative_to(built["source"]).as_posix(): r
        for r in built["sources"] + built["python_sources"]
    }
    for name in module.THERMO_SOURCE:
        require(
            replay.oracle.sha(thermo / "source" / name) == sources[name]["sha256"],
            "Wrong thermodynamic source",
        )
    for directory, count, failures in (
        ("models", 7175, "failed_comparisons"),
        ("model_supplement", 449, "failed"),
    ):
        path = args.science / directory / "manifest.json"
        value = replay.load_manifest(path)
        require(
            value["complete"]
            and value["numerical_contract_checks_passed"]
            and value["comparisons"] == count
            and not value[failures],
            "Independent model equations failed",
        )
        module.verify_all(value)
        evidence.append(replay.oracle.file_record(path))
    for threads in (2, 4, 6):
        path = args.science / f"cli_comparison_{threads}/manifest.json"
        value = replay.load_manifest(path)
        require(
            value["complete"]
            and value["exact_equivalence"]
            and value["cases"] == 93
            and not value["differences"],
            "CLI comparison failed",
        )
        module.verify_all(value)
        evidence.append(replay.oracle.file_record(path))
    legacy_path = args.science / "legacy_comparison.json"
    legacy = read(legacy_path)
    require(
        legacy["total_rows"] == 128021
        and legacy["exact_legacy_replay_pass"]
        and not legacy["failures"],
        "Exact replay failed",
    )
    for path, digest in legacy["verified_files"].items():
        require(replay.oracle.sha(Path(path)) == digest, "Raw legacy file changed")
    engineering_path = args.science / "execution/engineering/journal.json"
    engineering = read(engineering_path)
    expected_checks = [
        "rust_format",
        "clippy",
        "rustdoc",
        "rust_tests",
        "screening_regression",
        "python_tests",
        "dbstep_upstream_tests",
        "archived_morfeus_evidence",
        "retrospective_screening",
    ]
    require(
        engineering["complete"]
        and [r["name"] for r in engineering["runs"]] == expected_checks
        and engineering["build"]
        == replay.oracle.file_record(args.build / "manifest.json")
        and all(r["returncode"] == 0 for r in engineering["runs"]),
        "Engineering checks incomplete",
    )
    archived_controller = engineering_path.with_name("controller.py")
    require(
        replay.oracle.sha(archived_controller) == engineering["controller"]["sha256"],
        "Executed engineering controller archive changed",
    )
    evidence.append(replay.oracle.file_record(archived_controller))
    for row in engineering["runs"]:
        for stream in ("stdout", "stderr"):
            replay.oracle.verify(row[stream])
    replay.manifest(args.output / "differences.json", changes)
    failures = len(changes) + sum(c["changed_values_or_structure"] for c in comparisons)
    evidence = list(module.records(evidence))
    evidence += [
        replay.oracle.file_record(p)
        for p in (
            Path(__file__),
            args.build / "manifest.json",
            legacy_path,
            engineering_path,
            args.output / "differences.json",
            args.science / "observations/manifest.json",
            args.science / "volume_reference/complete.json",
            args.science / "all_bins/complete.json",
            args.science / "fixed_union/complete.json",
        )
    ]
    result = dict(
        kind="dbstep_resolution_legacy_scientific_review",
        build=replay.oracle.file_record(args.build / "manifest.json"),
        evidence=evidence,
        complete_reference_coverage=True,
        changed_legacy_reference_values=failures,
        legacy_reference_preservation_pass=failures == 0,
        scientific_global_pass=False,
        release_admitted=False,
        interpretation=(
            "Legacy preservation only; new DBSTEP analyses and release policy are "
            "separate gates. No negative historical result is upgraded."
        ),
    )
    module.verify_all(result)
    replay.manifest(args.output / "review.json", result)
    print(
        json.dumps(
            {
                "changed_legacy_reference_values": failures,
                "legacy_reference_preservation_pass": failures == 0,
            }
        )
    )
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
