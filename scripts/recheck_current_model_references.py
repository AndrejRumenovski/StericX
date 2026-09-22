#!/usr/bin/env python3
r"""Run unchanged independent model references against fresh CLI/domain captures.

  .venv/bin/python3 scripts/recheck_current_model_references.py \
      --observations /tmp/live-observations --cli-run /tmp/cli/runs/candidate \
      --cli-build-snapshot /tmp/build-snapshot --output /tmp/model-references \
      --components math reactions --kraken-recheck /tmp/independent-kraken

The Kraken analysis phase may be complete while its expensive references are
still running: its successful command receipt, fresh SUT stream and every
analysis artifact are verified. No historical raw model outputs are substituted
for captures. Historical serialized models used as actual CLI inputs are checked
against newly fitted models, allowing only their process creation timestamp.
Subprocess/SUT execution is forbidden inside all independent reference workers.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from datetime import UTC, datetime
from decimal import Decimal, InvalidOperation
from itertools import zip_longest
from pathlib import Path

import recheck_current_scientific_references as common

MATH = (
    "models_audit",
    "models_screen",
    "models_domain",
    "models_extended",
    "models_discrepancies",
    "models_outer_folds",
    "models_evaluation_edges",
)
REACTIONS = ("models_reactions", "models_finalize_evidence")
FIXED_DATA = (
    "data/reactions.sigpack",
    "data/reactions_raw.csv",
    "docs/study_001/published_model.json",
    "docs/study_002/official_kraken_comparison.csv",
    "docs/study_002/native_model.json",
    "docs/study_003/official_kraken_lmo_comparison.csv",
    "docs/study_003/quantum_model.json",
    "docs/study_011/rankings.csv",
)
SOURCES = (
    "ni_hda_reaction_data.csv",
    "ni_hda_kraken_features.csv",
    "crosscoupling_preprint_si.txt",
)
METADATA_INPUTS = {
    "objective_minimize.json",
    "objective_maximize_magnitude.json",
    "objective_unspecified.json",
    "marginal_ci_counterexample_model.json",
}
LIMITS = [
    "This reruns numerical references; it does not establish experimental validity.",
    "Known failures remain. Exact historical equality is not a scientific PASS.",
    "Study001/002/003 fits and Study011 rankings are independently reanalyzed "
    "historical evidence, not new native fits or new chemical outcomes.",
    "Formula-identity helper observations are not rerun; no frozen study SUT "
    "is executed by this helper.",
    "All33 model CLI captures are staged. Independent screening equations use "
    "the final v2 campaign; earlier six v1 captures are retained without a "
    "separate repeated reference calculation. Two deck exports are covered "
    "by the CLI oracle, not a new independent model equation here.",
    "Reaction references retain primary SI/source-version and ensemble limits. "
    "Their current-descriptor input must come from new candidate Kraken analysis.",
]


def load_cli(
    run: Path, built: dict, sealed, snapshot: Path
) -> tuple[dict, dict, object]:
    root = run.parent.parent
    helper = common.REPO / "scripts/check_scientific_cli_equivalence.py"
    spec = importlib.util.spec_from_file_location("candidate_cli_verifier", helper)
    if spec is None or spec.loader is None:
        raise ValueError("CLI verifier unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    plan = common.read(root / "manifest.json")
    identity = common.read(run / "identity.json")
    complete = common.read(run / "complete.json")
    if identity["manifest_sha256"] != common.sha(root / "manifest.json"):
        raise ValueError("CLI plan changed after observation")
    if identity["harness_sha256"] != common.sha(helper):
        raise ValueError("CLI observation helper changed")
    module.verify(plan)
    if common.sha(Path(identity["binary"])) != identity["binary_sha256"]:
        raise ValueError("CLI binary changed after observation")
    observed = {c["id"]: module.fingerprint(run / c["id"]) for c in plan["cases"]}
    observed["python"] = common.sha(run / "python/python.json")
    if observed != complete["results"]:
        raise ValueError("CLI captures changed after observation")
    receipt = common.bind_cli_build(snapshot, identity, built, sealed)
    return (
        plan,
        {
            "identity": identity,
            "build_receipt": receipt,
            "receipts": [
                common.record(run / n) for n in ("identity.json", "complete.json")
            ],
            "plan": common.record(root / "manifest.json"),
            "helper": common.record(helper),
        },
        module,
    )


def without_creation(value):
    value = json.loads(json.dumps(value))
    if isinstance(value, dict) and isinstance(value.get("created"), dict):
        value["created"].pop("created_utc", None)
    return value


def stage_raw(work: Path, sealed, run: Path, plan: dict, observed: dict) -> dict:
    raw = work / "models/raw"
    raw.mkdir(parents=True)
    model_cases = [c for c in plan["cases"] if c["id"].startswith("models_")]
    if len(model_cases) != 33:
        raise ValueError("Expected all33 audited model CLI cases")
    replacements = {}
    for case in model_cases:
        name = case["id"].removeprefix("models_")
        src = run / case["id"]
        dst = raw / name
        dst.mkdir()
        for file in ("stdout", "stderr"):
            common.candidate_copy(common.record(src / file), dst / file)
        for artifact in sorted((src / "artifacts").iterdir()):
            if not artifact.is_file():
                raise ValueError(f"Unexpected nested model artifact: {artifact}")
            common.candidate_copy(common.record(artifact), dst / artifact.name)
        status = common.read(src / "result.json")["returncode"]
        common.write_new(
            dst / "command.json",
            {
                "argv": case["argv"],
                "exit_code": status,
                "source_capture": str(src),
            },
        )
        replacements.update(
            {new: old for old, new in case["output_path_replacements"].items()}
        )
    lane = observed["lanes"].get("model_domain")
    if lane is None or lane["returncode"] != 0 or lane["coverage"]["rows"] != 33:
        raise ValueError("All33 fresh model-domain observer records are required")
    common.candidate_copy(lane["stdout"], raw / "domain_observer.jsonl")
    common.candidate_copy(lane["stderr"], raw / "domain_observer.stderr")
    if common.sha(Path(lane["requests"]["path"])) != common.sha(
        sealed.check("models/inputs/domain_requests.jsonl")
    ):
        raise ValueError("Model-domain observer inputs differ from audited requests")
    # Screens were actually run on these fixed serialized models. A changed fit
    # would require a new downstream capture, not a mismatched reference pairing.
    model_checks = []
    for name in ("synthetic_linear", "native_ni_hda"):
        reference = sealed.check(f"models/raw/{name}/portable.json")
        candidate = raw / name / "portable.json"
        same = without_creation(common.read(reference)) == without_creation(
            common.read(candidate)
        )
        model_checks.append(
            {
                "name": name,
                "same_scientific_content": same,
                "serialized_cli_fixture": common.record(reference),
                "fresh_fit_artifact": common.record(candidate),
            }
        )
        if not same:
            raise ValueError(
                f"New fit differs from downstream fixture: {name}; "
                "recapture downstream commands"
            )
    return {
        "model_cli_cases": [c["id"] for c in model_cases],
        "domain_records": 33,
        "path_aliases": replacements,
        "downstream_fixture_binding": model_checks,
    }


def stage_kraken_analysis(parent: Path, work: Path, observed: dict) -> dict:
    started = common.read(parent / "started.json")
    binding = started["candidate_observation_manifest"]
    source_observations = common.read(common.verify(binding))
    if source_observations["build_manifest"] != observed["build_manifest"]:
        raise ValueError("Kraken analysis belongs to a different candidate build")
    receipts = []
    for path in sorted((parent / "commands").glob("*.json")):
        row = common.read(path)
        if Path(row["argv"][1]).name == "kraken_analyze.py":
            if row["returncode"] != 0 or row["timed_out"]:
                raise ValueError("Fresh Kraken analysis did not finish successfully")
            common.verify(row["stdout"])
            common.verify(row["stderr"])
            receipts.append(common.record(path))
    if len(receipts) != 1:
        raise ValueError(
            "Exactly one completed new Kraken analysis receipt is required"
        )
    root = parent / "recalculation/kraken"
    sut = common.read(root / "sut/manifest_frozen.json")
    if sut["candidate_lane"] != observed["lanes"]["kraken"]:
        raise ValueError("Kraken phase is not based on this candidate's fresh stream")
    common.verify(sut["raw_output_files"]["stdout.jsonl"])
    analysis = common.read(root / "analysis/manifest.json")
    if analysis["frozen_SUT_manifest_sha256"] != common.sha(
        root / "sut/manifest_frozen.json"
    ):
        raise ValueError("Kraken analysis source binding changed")
    for name, digest in analysis["outputs"].items():
        if common.sha(root / "analysis" / name) != digest:
            raise ValueError(f"New Kraken analysis artifact changed: {name}")
    # The sealed auxiliary table is this exact projection, but kraken_analyze.py
    # itself does not emit it. Reconstruct it from NEW analysis, never old rows.
    import pandas as pd

    source = root / "analysis/comparisons.csv"
    frame = pd.read_csv(source)
    projected = frame.loc[
        frame.metric.eq("percent_buried_volume_min"),
        ["molecule_id", "sut", "reference", "conformers"],
    ]
    if projected.empty or projected.molecule_id.duplicated().any():
        raise ValueError("New Kraken percentage minima are empty or duplicated")
    table = work / "kraken/analysis/percent_buried_volume_min_by_ligand.csv"
    table.parent.mkdir(parents=True)
    with table.open("x") as stream:
        projected.to_csv(stream, index=False)
    return {
        "parent_started": common.record(parent / "started.json"),
        "analysis_receipt": receipts[0],
        "analysis_manifest": common.record(root / "analysis/manifest.json"),
        "fresh_descriptor_table": common.record(table),
        "projection_source": common.record(source),
        "projection": "metric=percent_buried_volume_min; "
        "columns molecule_id,sut,reference,conformers; pandas to_csv(index=False)",
        "scope": "New analysis complete; full Kraken references may still be active",
    }


def worker(work: Path, name: str, binary: Path) -> None:
    """Redirect observations/paths, not equations; fail any attempted SUT call."""
    sys.path.insert(0, str(work / "scripts"))
    import models_audit as core

    # platform.platform() can lazily call uname. Cache runtime identity before
    # forbidding subprocesses; no SUT is invoked by this metadata acquisition.
    core.platform.platform()
    original_freeze = core.freeze
    metadata_log = work.parent / "input_metadata_normalizations.jsonl"

    def guarded_freeze(path, data):
        path = Path(path)
        if not path.resolve().is_relative_to(work):
            raise ValueError(f"Reference attempted write outside its new tree: {path}")
        if path.exists() and path.read_bytes() != data and path.name in METADATA_INPUTS:
            a, b = json.loads(path.read_bytes()), json.loads(data)
            if without_creation(a) == without_creation(b):
                with metadata_log.open("a") as stream:
                    stream.write(
                        json.dumps(
                            {
                                "path": str(path.relative_to(work)),
                                "allowed_field": "/created/created_utc",
                                "reason": "creation time; scientific fields exact",
                            }
                        )
                        + "\n"
                    )
                return
        original_freeze(path, data)

    def no_sut(*args, **kwargs):
        raise RuntimeError(
            "Missing fresh capture: SUT/subprocess execution forbidden in worker"
        )

    core.freeze = guarded_freeze
    core.run = no_sut
    core.BIN = binary
    core.subprocess.Popen = no_sut
    core.subprocess.run = no_sut
    module = core if name == "models_audit" else __import__(name)
    if name == "models_reactions":
        module.ni_hda()
        module.crosscoupling()
    else:
        module.main()


def normalize(value, aliases: dict, cli_module, old: Path, work: Path, key=""):
    if isinstance(value, dict):
        return {
            k: normalize(v, aliases, cli_module, old, work, k) for k, v in value.items()
        }
    if isinstance(value, list):
        return [normalize(v, aliases, cli_module, old, work, key) for v in value]
    if isinstance(value, str):
        for actual, original in sorted(aliases.items(), key=lambda p: -len(p[0])):
            value = value.replace(actual, original)
        if key == "stdout":
            value = cli_module.scientific_stdout(value.encode()).decode()
        return common.normalize(value, old, work)
    return value


def compare(sealed, work: Path, output: Path, aliases: dict, cli_module) -> dict:
    summary = []
    with (output / "historical_value_differences.jsonl").open("x") as dst:
        for path in sorted((work / "models/results").iterdir()):
            if path.suffix not in {".csv", ".json"}:
                continue
            relative = str(path.relative_to(work))
            before_path = sealed.check(relative)
            changes = count = 0
            with before_path.open() as before, path.open() as after:
                rows = (
                    [(json.load(before), json.load(after))]
                    if path.suffix == ".json"
                    else zip_longest(csv.DictReader(before), csv.DictReader(after))
                )
                for count, (a, b) in enumerate(rows, 1):
                    a, b = [
                        normalize(v, aliases, cli_module, sealed.audit, work)
                        for v in (a, b)
                    ]
                    for delta in common.differences(a, b):
                        if path.suffix == ".csv" and "historical_present" not in delta:
                            x, y = delta["historical"], delta["candidate"]
                            try:
                                if (
                                    isinstance(x, str)
                                    and isinstance(y, str)
                                    and Decimal(x).is_finite()
                                    and Decimal(y).is_finite()
                                    and Decimal(x) == Decimal(y)
                                ):
                                    continue
                            except InvalidOperation:
                                pass
                        changes += 1
                        dst.write(
                            json.dumps(
                                {"file": relative, "record": count, **delta},
                                allow_nan=False,
                            )
                            + "\n"
                        )
            summary.append(
                {
                    "file": relative,
                    "records_compared": count,
                    "changed_values_or_structure": changes,
                    "historical": common.record(before_path),
                    "candidate": common.record(path),
                }
            )
    return {
        "files": summary,
        "changed_values_or_structure": sum(
            r["changed_values_or_structure"] for r in summary
        ),
        "normalizations": [
            "recorded CLI output path aliases",
            "audit root paths",
            "exactly equal decimal CSV formatting",
            "CLI oracle's explicit process stdout metric keys",
        ],
        "no_scientific_tolerance": True,
    }


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        worker(Path(sys.argv[2]).resolve(strict=True), sys.argv[3], Path(sys.argv[4]))
        return
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    for name in ("observations", "cli-run", "cli-build-snapshot", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument(
        "--components", choices=["math", "reactions"], nargs="+", default=["math"]
    )
    parser.add_argument("--kraken-recheck", type=Path)
    parser.add_argument("--audit", type=Path, default=common.AUDIT)
    parser.add_argument("--command-timeout", type=float, default=7200)
    args = parser.parse_args()
    if len(set(args.components)) != len(args.components) or args.command_timeout <= 0:
        parser.error("Components must be unique and timeout positive")
    if "reactions" in args.components and args.kraken_recheck is None:
        parser.error("Reaction references require --kraken-recheck with new analysis")
    audit, capture, cli, snapshot = [
        p.resolve(strict=True)
        for p in (args.audit, args.observations, args.cli_run, args.cli_build_snapshot)
    ]
    output = args.output.resolve()
    protected = [audit, capture, cli.parent.parent, snapshot]
    if args.kraken_recheck:
        protected.append(args.kraken_recheck.resolve(strict=True))
    if any(output.is_relative_to(p) or p.is_relative_to(output) for p in protected):
        raise ValueError("Output must be separate from audit and capture trees")
    if cli.is_relative_to(audit):
        raise ValueError("Fresh CLI captures required")
    output.mkdir(parents=True, exist_ok=False)
    work = output / "recalculation"
    work.mkdir()
    sealed = common.SealedInputs(audit)
    start = {
        "kind": "new_independent_model_reference_recheck",
        "created_utc": datetime.now(UTC).isoformat(),
        "helper": common.record(Path(__file__)),
        "shared_helper": common.record(Path(common.__file__)),
        "components": args.components,
        "candidate_observation_manifest": common.record(capture / "manifest.json"),
        "limitations": LIMITS,
    }
    common.write_new(output / "started.json", start)
    result = {**start, "complete": False, "scientific_pass_claimed": False}
    try:
        observed, built, _ = common.observations(capture, audit, [])
        result["reference_environment"] = common.reference_environment(sealed)
        expected_sklearn = common.read(sealed.check("manifest_initial.json"))["python"][
            "packages"
        ]["scikit-learn"]
        if common.importlib.metadata.version("scikit-learn") != expected_sklearn:
            raise ValueError("Scikit-learn differs from the pinned audit environment")
        plan, provenance, cli_module = load_cli(cli, built, sealed, snapshot)
        result["cli_provenance"] = provenance
        result["fresh_capture_coverage"] = stage_raw(work, sealed, cli, plan, observed)
        names = list(MATH if "math" in args.components else ()) + list(
            REACTIONS if "reactions" in args.components else ()
        )
        for name in dict.fromkeys(["models_audit", *names]):
            sealed.copy("scripts/" + name + ".py", work / "scripts" / (name + ".py"))
        for relative in sealed.files:
            if relative.startswith("models/inputs/") and not relative.endswith(
                (
                    "evaluate_edges_manifest_before_execution.json",
                    "current_native_percent_buried_volume_min.csv",
                )
            ):
                sealed.copy(relative, work / relative)
        for name in SOURCES:
            sealed.copy("models/sources/" + name, work / "models/sources" / name)
        for name in FIXED_DATA:
            sealed.copy("frozen/repository/" + name, work / "frozen/repository" / name)
        if "reactions" in args.components:
            result["new_kraken_analysis"] = stage_kraken_analysis(
                args.kraken_recheck.resolve(strict=True), work, observed
            )
        commands = [
            [
                str(Path(__file__).resolve()),
                "_worker",
                str(work),
                name,
                provenance["identity"]["binary"],
            ]
            for name in names
        ]
        common.write_new(
            output / "reference_plan.json",
            {
                "commands": commands,
                "inputs": list(sealed.used.values()),
                "forbidden": "All subprocess/SUT execution inside reference workers",
            },
        )
        result["commands"] = common.run_commands(commands, output, args.command_timeout)
        result["historical_comparison"] = compare(
            sealed,
            work,
            output,
            result["fresh_capture_coverage"]["path_aliases"],
            cli_module,
        )
        for item in sealed.used.values():
            common.verify(item)
        for key in ("helper", "shared_helper", "candidate_observation_manifest"):
            common.verify(start[key])
        common.verify(sealed.manifest)
        load_cli(cli, built, sealed, snapshot)
        result["complete"] = True
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        result["original_sealed_manifest"] = sealed.manifest
        result["used_sealed_files"] = list(sealed.used.values())
        common.write_new(output / "result.json", result)
        files = [common.record(p) for p in sorted(output.rglob("*")) if p.is_file()]
        common.write_new(
            output / "manifest.json",
            {"kind": "independent_model_reference_files", "files": files},
        )
    print(
        json.dumps(
            {
                "output": str(output),
                "complete": result["complete"],
                "changed_values_or_structure": result.get(
                    "historical_comparison", {}
                ).get("changed_values_or_structure"),
                "scientific_pass_claimed": False,
                "error": result.get("error"),
            },
            indent=2,
        )
    )
    if not result["complete"]:
        raise SystemExit(1)
    if result["historical_comparison"]["changed_values_or_structure"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
