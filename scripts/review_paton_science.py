#!/usr/bin/env python3
"""Compare complete fresh references with the independently admitted old evidence.

This admits only descriptor-preservation timing. Known DBSTEP convention and
boundary differences remain failures of their original strict comparison gates.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import replay_scientific_remediation as replay

ROOT = Path(__file__).resolve().parents[1]
PATON = ROOT / ".stericx/paton"
FROZEN = ROOT / ".stericx/profiling/scientifically_validated_optimization"


def read(path):
    return json.loads(Path(path).read_text())


def require(value, message):
    if not value:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    controller = FROZEN / "candidate_reference_controller_v1/controller.py"
    spec = importlib.util.spec_from_file_location(
        "paton_reference_comparator", controller
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    frozen = replay.load_manifest(FROZEN / "accurate_baseline_v1/manifest.json")
    baseline = module.Baseline(frozen)
    current_build, _ = replay.verify_build(PATON / "final_build/manifest.json")
    evidence = [replay.oracle.file_record(controller)]
    changes, comparisons = [], []
    for label, components in [
        ("inferred", ["geometry", "kraken"]),
        ("topology", ["kraken"]),
    ]:
        old = module.BASE / f"references_v4_{label}"
        new = PATON / f"final_references_{label}_v2"
        baseline.check(old / "result.json")
        receipt = replay.load_manifest(new / "result.json")
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

        destination = output / label
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
        evidence.extend(
            replay.oracle.file_record(p)
            for p in [new / "result.json", destination / "comparison.json"]
        )

    volume = read(PATON / "final_volume_reference_v2/complete.json")
    require(
        volume["complete"]
        and volume["limit"] is None
        and volume["requests"] == volume["compared"] == volume["source_rows"] == 31721
        and volume["reference_errors"] == 0
        and volume["complete_source_stream_validated"],
        "Incomplete full volume reference",
    )
    bins = read(PATON / "final_all_bins/complete.json")
    require(
        (bins["rows"], bins["frames"], bins["bins"]) == (31721, 95163, 1141956),
        "Incomplete regional bins",
    )
    union = read(PATON / "final_fixed_union/complete.json")
    old_union = read(
        baseline.check(module.BASE / "fixed_reference_union_v4/complete.json")
    )
    require(
        union["complete"]
        and union["union_count"] == 2517
        and union["memberships"] == old_union["memberships"],
        "Fixed union changed",
    )
    file_pairs = [
        (
            module.BASE / "buried_volume_reference_v4",
            PATON / "final_volume_reference_v2",
            ["reference_results.jsonl"],
        ),
        (
            module.GEOMETRY / "all_bins_reference_v4",
            PATON / "final_all_bins",
            ["per_id_residuals.jsonl", "metrics.json", "top20.json"],
        ),
        (
            module.BASE / "thermodynamics/final_v4",
            PATON / "final_thermodynamics",
            [
                "summary.json",
                "kinetics_comparison.json",
                "aggregation_comparison.json",
                "python.json",
            ],
        ),
    ]
    for lane in ["inferred", "topology"]:
        require(
            union["lanes"][lane]["rows"] == 2517
            and union["lanes"][lane]["coverage"]
            == old_union["lanes"][lane]["coverage"],
            "Fixed variant coverage changed",
        )
        file_pairs.append(
            (
                module.BASE / "fixed_reference_union_v4" / lane,
                PATON / "final_fixed_union" / lane,
                [
                    "selected_inputs.jsonl",
                    "reference_raw.jsonl",
                    "mean_recombination.jsonl",
                ],
            )
        )
    for old, new, names in file_pairs:
        for name in names:
            evidence.append(
                module.compare_file(old / name, new / name, baseline, changes, old, new)
            )

    # Bind copied thermodynamic executables and all eight copied source files.
    thermo = PATON / "final_thermodynamics"
    for name, key in [("stericx", "native"), ("observer", "observer")]:
        require(
            replay.oracle.sha(thermo / name)
            == current_build["binaries"][key]["sha256"],
            "Wrong thermodynamic executable",
        )
    source = {
        Path(r["path"]).relative_to(current_build["source"]).as_posix(): r
        for r in current_build["sources"] + current_build["python_sources"]
    }
    for name in module.THERMO_SOURCE:
        require(
            replay.oracle.sha(thermo / "source" / name) == source[name]["sha256"],
            "Wrong thermodynamic source",
        )

    for directory, count, failures in [
        ("final_models", 7175, "failed_comparisons"),
        ("final_models_supplement_v2", 449, "failed"),
    ]:
        value = replay.load_manifest(PATON / directory / "manifest.json")
        require(
            value["complete"]
            and value["numerical_contract_checks_passed"]
            and value["comparisons"] == count
            and not value[failures],
            "Independent model equations failed",
        )
        module.verify_all(value)
        evidence.append(replay.oracle.file_record(PATON / directory / "manifest.json"))

    legacy = read(
        ROOT / "docs/paton_recommendations/validation/final_replay_comparison.json"
    )
    require(
        legacy["total_rows"] == 128021
        and legacy["exact_legacy_replay_pass"]
        and not legacy["failures"],
        "Legacy exact replay failed",
    )
    for path, digest in legacy["verified_files"].items():
        require(replay.oracle.sha(Path(path)) == digest, "Legacy raw evidence changed")
    for thread in [2, 4, 6]:
        path = PATON / f"cli_comparison_{thread}/manifest.json"
        value = replay.load_manifest(path)
        require(
            value["complete"]
            and value["exact_equivalence"]
            and value["cases"] == 93
            and not value["differences"],
            "CLI thread comparison failed",
        )
        module.verify_all(value)
        evidence.append(replay.oracle.file_record(path))

    density = read(PATON / "final_density/summary.json")
    sterimol = read(PATON / "final_dbstep_sterimol/summary.json")
    dbvolume = read(PATON / "final_dbstep_volume/summary.json")
    profiles = read(PATON / "final_profiles/summary.json")
    populations = read(PATON / "final_populations/report.json")
    require(
        density["admitted"] and density["completed"] == 12 and not density["failures"],
        "Density reference failed",
    )
    require(
        sterimol["successful_cases"] == 61
        and sterimol["conditional_numerical_gate_pass"]
        and not sterimol["errors"],
        "Sterimol conditional reference failed",
    )
    require(
        dbvolume["cases"] == 69
        and dbvolume == read(ROOT / "docs/validation/dbstep/volume_v2/summary.json"),
        "DBSTEP volume result or failure set changed",
    )
    require(
        profiles["profile_preservation_pass"]
        and profiles["positions"] == 20
        and profiles["reference_comparisons"] == 16
        and not profiles["failures"],
        "Profile reference failed",
    )
    require(not populations["failures"], "Population reference failed")
    ensemble = read(PATON / "final_ensemble_reference/summary.json")
    require(
        ensemble["conformers"] == 56
        and ensemble["ligands"] == 11
        and not ensemble["reference_failures"]
        and not ensemble["conformer_investigations"]
        and not ensemble["ensemble_investigations"],
        "Ensemble independent reference failed",
    )
    for directory, name in [
        ("final_density", "summary.json"),
        ("final_dbstep_sterimol", "summary.json"),
        ("final_dbstep_volume", "summary.json"),
        ("final_profiles", "summary.json"),
        ("final_populations", "report.json"),
        ("final_ensemble_reference", "summary.json"),
        ("final_volume_reference_v2", "complete.json"),
        ("final_all_bins", "complete.json"),
        ("final_fixed_union", "complete.json"),
    ]:
        evidence.append(replay.oracle.file_record(PATON / directory / name))
    replay.manifest(output / "differences.json", changes)
    failures = len(changes) + sum(c["changed_values_or_structure"] for c in comparisons)
    # Flatten comparison records to the same verifiable file-record contract.
    evidence = list(module.records(evidence))
    evidence.extend(
        replay.oracle.file_record(p)
        for p in [
            Path(__file__),
            output / "differences.json",
            ROOT / "docs/paton_recommendations/validation/final_replay_comparison.json",
        ]
    )
    result = {
        "kind": "paton_independent_scientific_review",
        "baseline_build": replay.oracle.file_record(
            PATON / "baseline/build/manifest.json"
        ),
        "candidate_build": replay.oracle.file_record(
            PATON / "final_build/manifest.json"
        ),
        "evidence": evidence,
        "complete_reference_coverage": True,
        "changed_legacy_reference_values": failures,
        "legacy_reference_preservation_pass": failures == 0,
        "dbstep_strict_volume_gate": dbvolume["finite_lattice_gate_pass"],
        "density_admitted": True,
        "descriptor_performance_admitted": failures == 0,
        "scientific_global_pass": False,
        "interpretation": (
            "Preserved admitted legacy results and new scoped validations. "
            "Known f32/f64 boundary, radii, angular, historical-provenance and "
            "predictive limitations remain. Strict DBSTEP volume gate remains "
            "false; this is not a release admission."
        ),
    }
    module.verify_all(result)
    replay.manifest(output / "review.json", result)
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "changed_legacy_reference_values",
                    "descriptor_performance_admitted",
                    "dbstep_strict_volume_gate",
                ]
            }
        )
    )
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
