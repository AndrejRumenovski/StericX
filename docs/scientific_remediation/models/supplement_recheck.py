#!/usr/bin/env python3
"""Close fitted-geometry and full corrected-screen coverage, without recapturing SUT.

Consumes an immutable successful model reference receipt and the actual CLI/demo
inputs it binds. Reference arrays come from packed training records, not stored
model geometry. Greedy panel selection uses a full independent pairwise distance
matrix and the public report's declared rank-normalized objective.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("model_recheck", HERE / "recheck.py")
base = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(base)


def geometry(check, core, name, argv, report_path, portable_path):
    data = np.fromfile(base.option(argv, "--data"), dtype="<f4").reshape(-1, 16)
    metadata = base.rows(base.option(argv, "--metadata"))
    selected_rows = [
        i for i, row in enumerate(metadata) if row["Dataset_Split"].lower() == "train"
    ]
    report = base.read(report_path)
    selected = report["selected_feature_indices"]
    x = core.expand(data[selected_rows])[:, selected]
    means, scales = x.mean(0), x.std(0)
    scales[scales < 1e-12] = 1
    expected = (x - means) / scales
    g = report["training_geometry"]
    check.equal(name + "/feature_indices", g["feature_indices"], selected)
    check.equal(name + "/observations", g["observations"], len(selected_rows))
    check.equal(name + "/parameters", g["parameters"], len(selected) + 1)
    check.equal(
        name + "/training_labels",
        g["training_labels"],
        [metadata[i]["Reaction_ID"] for i in selected_rows],
    )
    check.equal(
        name + "/training_points_count",
        len(g["standardized_training_points"]),
        len(expected),
    )
    check.equal(name + "/means_count", len(g["means"]), len(means))
    check.equal(name + "/scales_count", len(g["scales"]), len(scales))
    for i in range(len(selected)):
        check.close(f"{name}/mean/{i}", g["means"][i], means[i])
        check.close(f"{name}/scale/{i}", g["scales"][i], scales[i])
    for i, row in enumerate(expected):
        check.equal(
            f"{name}/point_width/{i}",
            len(g["standardized_training_points"][i]),
            len(row),
        )
        for j, value in enumerate(row):
            check.close(
                f"{name}/point/{i}/{j}", g["standardized_training_points"][i][j], value
            )
    if portable_path:
        model = base.read(portable_path)
        check.equal(name + "/portable_geometry_exact", model["training_geometry"], g)
        check.equal(
            name + "/portable_selected_exact",
            model["selected_feature_indices"],
            selected,
        )
    return dict(
        name=name,
        training_rows=selected_rows,
        means=means.tolist(),
        scales=scales.tolist(),
        standardized_points=expected.tolist(),
    )


def screen(check, core, name, argv, stdout):
    model_path = Path(argv[1])
    model, report = base.read(model_path), base.read(stdout)
    library_path = Path(base.option(argv, "--library"))
    library = base.rows(library_path)
    records = np.fromfile(model_path.parent / "reactions.sigpack", dtype="<f4").reshape(
        -1, 16
    )
    check.equal(name + "/library_record_count", len(library), len(records))
    ids = [row["Reaction_ID"] for row in library]
    check.equal(name + "/unique_source_ids", len(set(ids)), len(ids))
    x = core.expand(records)
    scores = x @ np.asarray(model["weights"], dtype=np.float32).astype(float)
    tested_path = base.option(argv, "--exclude-tested")
    tested = set()
    if tested_path:
        tested = {row["Reaction_ID"] for row in base.rows(tested_path)}
        check.equal(
            name + "/tested_file_sha",
            report["exclusion"]["sha256"],
            base.oracle.sha(Path(tested_path)),
        )
        check.equal(
            name + "/tested_matched",
            report["exclusion"]["matched"],
            len(set(ids) & tested),
        )
        check.equal(
            name + "/tested_unresolved",
            report["exclusion"]["unresolved"],
            sorted(tested - set(ids)),
        )
    eligible = [i for i, ident in enumerate(ids) if ident not in tested]
    # These exact documented examples declare maximize, have all required
    # finite descriptors, and do not request the optional domain-only filter.
    check.equal(
        name + "/objective", model["inference"]["response"]["optimization"], "maximize"
    )
    check.equal(name + "/no_domain_filter_requested", "--in-domain-only" in argv, False)
    ranked = sorted(eligible, key=lambda i: (-scores[i], ids[i], i))
    top = int(base.option(argv, "--top", len(ranked)))
    take = min(top, len(ranked))
    selection = []
    if "--diverse" in argv:
        # Independently reconstruct standardization from the packed training
        # matrix and metadata, rather than using the model's stored geometry.
        g = model["training_geometry"]
        training_labels = set(g["training_labels"])
        training = [i for i, ident in enumerate(ids) if ident in training_labels]
        cols = model["selected_feature_indices"]
        means, scales = x[training][:, cols].mean(0), x[training][:, cols].std(0)
        scales[scales < 1e-12] = 1
        points = (x[ranked][:, cols] - means) / scales
        distances = cdist(points, points, metric="euclidean")
        scale = float(np.sqrt(np.sum(np.ptp(points, axis=0) ** 2)))
        desirability = np.linspace(1, 0, len(ranked)) if len(ranked) > 1 else np.ones(1)
        weight = float(base.option(argv, "--diversity-weight"))
        chosen = []
        while len(chosen) < take:
            trials = []
            for position, source_index in enumerate(ranked):
                if position in chosen:
                    continue
                nearest = float(distances[position, chosen].min()) if chosen else None
                normalized = (
                    min(nearest / scale, 1) if nearest is not None and scale else 0.0
                )
                objective = (
                    float(desirability[position])
                    if not chosen
                    else float(
                        (1 - weight) * desirability[position] + weight * normalized
                    )
                )
                trials.append(
                    (
                        -objective,
                        ids[source_index],
                        source_index,
                        position,
                        nearest,
                        normalized,
                    )
                )
            best = min(trials)
            chosen.append(best[3])
            selection.append(
                dict(
                    ligand=best[1],
                    step=len(chosen),
                    ranking_position=best[3] + 1,
                    desirability=float(desirability[best[3]]),
                    separation=best[4],
                    separation_normalized=best[5],
                    objective=-best[0],
                )
            )
        expected_ids = [ids[ranked[position]] for position in chosen]
        check.close(name + "/diversity_scale", report["diversity_scale"], scale)
    else:
        expected_ids = [ids[i] for i in ranked[:take]]
    check.equal(
        name + "/full_returned_id_order",
        [hit["ligand"] for hit in report["hits"]],
        expected_ids,
    )
    check.equal(name + "/library_size", report["library_size"], len(library))
    check.equal(name + "/screened_eligible", report["screened"], len(eligible))
    check.equal(name + "/returned", report["returned"], take)
    check.equal(name + "/actual_hit_count", len(report["hits"]), take)
    check.equal(name + "/skipped", report["skipped"], 0)
    check.equal(name + "/no_missing_descriptor_exclusions", report["excluded"], [])
    check.equal(name + "/no_domain_filtered", report["domain_filtered"], [])
    for i, expected in enumerate(selection):
        actual = report["hits"][i]["selection"]
        for key, value in expected.items():
            if key == "ligand":
                continue
            if value is None or key in ("step", "ranking_position"):
                check.equal(f"{name}/selection/{i}/{key}", actual[key], value)
            else:
                check.close(f"{name}/selection/{i}/{key}", actual[key], value)
    return dict(
        name=name,
        eligible_ids=[ids[i] for i in eligible],
        independent_order=expected_ids,
        independent_selection=selection,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prior_path = args.reference_run.resolve() / "manifest.json"
    prior = base.replay.load_manifest(prior_path)
    if not prior["complete"] or not prior["numerical_contract_checks_passed"]:
        raise ValueError(
            "Requires successful complete corrected model reference receipt"
        )
    base.verify_records(prior)
    bound = {Path(item["path"]).resolve(): item for item in prior["bindings"]}
    if bound[HERE / "recheck.py"]["sha256"] != base.oracle.sha(HERE / "recheck.py"):
        raise ValueError("Executing reference helper differs from prior receipt")
    build_path = next(
        p
        for p in bound
        if p.parent.name.startswith("validated_build_") and p.name == "manifest.json"
    )
    base.replay.verify_build(build_path)
    cli_manifest = next(
        p
        for p in bound
        if p.name == "manifest.json" and p.parent.name.startswith("corrected_v")
    )
    cli_root = cli_manifest.parent
    plan_path = next(
        p for p in bound if p.parent.name == "cli_v1" and p.name == "manifest.json"
    )
    demo_path = next(p for p in bound if p.name == "receipt.json")
    output = base.replay.new_output(args.output)
    helper_record = base.oracle.file_record(Path(__file__))
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
            geometry(
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
                geometry(
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
                screen(check, core, name, argv, Path(command["stdout"]["path"]))
            )
    base.save(output / "comparisons.json", check.items)
    base.save(output / "independent_results.json", details)
    base.verify_records(prior)
    base.refs.verify(helper_record)
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
