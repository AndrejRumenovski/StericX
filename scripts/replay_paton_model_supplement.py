#!/usr/bin/env python3
"""Run unchanged model supplement equations using explicit receipt-bound paths.

The historical controller discovers inputs by directory-name prefixes. This
adapter removes that naming assumption; numerical helpers and gates are reused.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / "docs/scientific_remediation/models"
REFERENCE = HERE / "supplement_recheck.py"
SPEC = importlib.util.spec_from_file_location("paton_supplement_reference", REFERENCE)
supplement = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(supplement)
base = supplement.base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["reference-run", "build", "cli-run", "cli-plan", "demo", "output"]:
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    prior_path = args.reference_run.resolve() / "manifest.json"
    prior = base.replay.load_manifest(prior_path)
    if not prior["complete"] or not prior["numerical_contract_checks_passed"]:
        raise ValueError("Requires a complete successful model reference receipt")
    base.verify_records(prior)
    bound = {Path(item["path"]).resolve(): item for item in prior["bindings"]}
    build_path = args.build.resolve() / "manifest.json"
    cli_manifest = args.cli_run.resolve() / "manifest.json"
    cli_root = cli_manifest.parent
    plan_path = args.cli_plan.resolve() / "manifest.json"
    demo_path = args.demo.resolve() / "receipt.json"
    for path in [build_path, cli_manifest, plan_path, demo_path, HERE / "recheck.py"]:
        if path not in bound or base.oracle.file_record(path) != bound[path]:
            raise ValueError(f"Input is not bound by the original receipt: {path}")
    base.replay.verify_build(build_path)
    output = base.replay.new_output(args.output)
    helper_record = base.oracle.file_record(Path(__file__))
    reference_helper = base.oracle.file_record(REFERENCE)
    shutil.copyfile(Path(__file__), output / "supplement_recheck.py")
    sealed = base.refs.SealedInputs(base.replay.AUDIT)
    core = base.load_reference(output / "reference", sealed)
    check, details = base.Checks(), []
    for case in base.read(plan_path)["cases"]:
        root = cli_root / case["id"]
        if not case["id"].startswith("models_") or case["argv"][0] != "fit":
            continue
        if base.read(root / "result.json")["returncode"]:
            continue
        portable = root / "artifacts/portable.json"
        details.append(
            supplement.geometry(
                check,
                core,
                case["id"],
                case["argv"],
                root / "artifacts/report.json",
                portable if portable.exists() else None,
            )
        )
    for command in base.read(demo_path)["commands"]:
        argv = command["executed_argv"][1:]
        if argv[0] == "fit":
            report = Path(base.option(argv, "--output"))
            details.append(
                supplement.geometry(
                    check,
                    core,
                    "corrected/" + report.parent.name,
                    argv,
                    report,
                    Path(base.option(argv, "--portable-model")),
                )
            )
        elif argv[0] == "screen" and base.option(argv, "--format") == "json":
            name = (
                "corrected/"
                + Path(command["redirected_stdout"]["path"]).parent.name
                + "/"
                + Path(command["redirected_stdout"]["path"]).stem
            )
            details.append(
                supplement.screen(
                    check, core, name, argv, Path(command["stdout"]["path"])
                )
            )
    base.save(output / "comparisons.json", check.items)
    base.save(output / "independent_results.json", details)
    base.verify_records(prior)
    base.refs.verify(helper_record)
    base.refs.verify(reference_helper)
    base.replay.verify_build(build_path)
    failed = [item for item in check.items if not item["passed"]]
    base.replay.manifest(
        output / "manifest.json",
        dict(
            kind="corrected_model_reference_coverage_supplement",
            complete=True,
            numerical_contract_checks_passed=not failed,
            scientific_global_pass=False,
            reference_receipt=base.oracle.file_record(prior_path),
            build=bound[build_path],
            helper=helper_record,
            source_reference=bound[HERE / "recheck.py"],
            unchanged_supplement_helper=reference_helper,
            comparisons=len(check.items),
            failed=failed,
            fitted_geometry_cases=sum("means" in item for item in details),
            corrected_screen_cases=sum("eligible_ids" in item for item in details),
            tolerance_policy=(
                "Original receipt f64 equation bounds: atol1e-12,rtol2e-12; "
                "exact schema,identities,counts,portable-vs-fit geometry."
            ),
            limitations=(
                "Same numerical scope as original receipt; "
                "no new experimental validity claim."
            ),
            files=[
                base.oracle.file_record(p)
                for p in sorted(output.rglob("*"))
                if p.is_file()
            ],
        ),
    )
    print(
        {"comparisons": len(check.items), "failed": len(failed), "output": str(output)}
    )
    if failed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
