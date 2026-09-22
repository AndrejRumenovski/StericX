#!/usr/bin/env python3
"""Independent equations on corrected, source-bound model observations.

This is a new staging/comparison adapter, not a replacement for the sealed audit.
The sealed NumPy/SciPy/sklearn reference functions are copied unchanged. Historical
serialized models are only inputs to their actual replayed screen commands; they
never replace newly fitted models or newly captured SUT observations.
"""

from __future__ import annotations

import argparse
import copy
import csv
import importlib.metadata
import importlib.util
import json
import math
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts"))
import capture_remediation_cli as capture  # noqa: E402
import replay_scientific_remediation as replay  # noqa: E402

oracle = replay.oracle
refs = replay.refs
AUDIT_SHA = "c31b655641a352402b14dc2d4535261a9de9c9004979c5e571d78e45cac82f03"


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def read(path):
    return json.loads(Path(path).read_text())


def records_in(value):
    if isinstance(value, dict):
        if {"path", "sha256"} <= value.keys():
            yield value
        else:
            for item in value.values():
                yield from records_in(item)
    elif isinstance(value, list):
        for item in value:
            yield from records_in(item)


def verify_records(value):
    for item in records_in(value):
        refs.verify(item)


def option(argv, name, default=None):
    return argv[argv.index(name) + 1] if name in argv else default


def rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


class Checks:
    """Save every residual; nonfinite/absent numeric observations never pass."""

    def __init__(self):
        self.items = []

    def equal(self, name, actual, expected):
        self.items.append(
            dict(name=name, actual=actual, expected=expected, passed=actual == expected)
        )

    def close(self, name, actual, expected, *, atol=1e-12, rtol=2e-12):
        finite = (
            isinstance(actual, (int, float))
            and not isinstance(actual, bool)
            and math.isfinite(actual)
            and math.isfinite(float(expected))
        )
        error = abs(float(actual) - float(expected)) if finite else None
        self.items.append(
            dict(
                name=name,
                actual=actual,
                expected=float(expected),
                absolute_error=error,
                atol=atol,
                rtol=rtol,
                passed=bool(finite and error <= atol + rtol * abs(expected)),
            )
        )


def verify_cli(run, build_path, built):
    receipt = replay.load_manifest(run / "manifest.json")
    if (
        receipt.get("kind") != "corrected_cli_capture"
        or not receipt["complete_capture"]
    ):
        raise ValueError("Complete corrected CLI capture required")
    if receipt["build_manifest"] != oracle.file_record(build_path):
        raise ValueError("CLI and reference build identities differ")
    if receipt["binary"] != built["binaries"]["native"]:
        raise ValueError("CLI native binary differs from frozen build")
    verify_records(receipt)
    if receipt["helper"]["sha256"] != oracle.sha(Path(capture.__file__)):
        raise ValueError("Executing CLI verifier differs from capture helper")
    plan_path = refs.verify(receipt["plan"])
    plan = capture.verify_plan(plan_path.parent)
    expected_files = {
        str(p.resolve())
        for p in run.rglob("*")
        if p.is_file()
        and p not in (run / "manifest.json", run / "manifest.json.sha256")
    }
    if {item["path"] for item in receipt["files"]} != expected_files:
        raise ValueError("CLI raw file inventory is incomplete or has extra files")
    if set(receipt["cases"]) != {case["id"] for case in plan["cases"]}:
        raise ValueError("CLI case coverage changed")
    for case in plan["cases"]:
        if receipt["cases"][case["id"]]["fingerprint"] != capture.cli.fingerprint(
            run / case["id"]
        ):
            raise ValueError("CLI raw fingerprint differs: " + case["id"])
    model_cases = [case for case in plan["cases"] if case["id"].startswith("models_")]
    if len(model_cases) != 33:
        raise ValueError("All 33 original model CLI requests are required")
    return model_cases


def verify_demo(path, build_path, built, bundle):
    data = read(path)
    verify_records(data)
    if data["input_bindings"] != data["input_bindings_after"]:
        raise ValueError("Documented demo inputs changed during execution")
    by_path = {Path(item["path"]).resolve(): item for item in data["input_bindings"]}
    if by_path[build_path]["sha256"] != oracle.sha(build_path):
        raise ValueError("Demo build mismatch")
    if data["binary"]["sha256"] != built["binaries"]["native"]["sha256"]:
        raise ValueError("Demo native executable mismatch")
    bundle_path = bundle / "manifest.json"
    if by_path[bundle_path]["sha256"] != oracle.sha(bundle_path):
        raise ValueError("Demo corrected input bundle mismatch")
    verify_records(read(bundle_path))
    if not data["all_commands_succeeded"] or len(data["commands"]) != 17:
        raise ValueError("All 17 documented corrected commands required")
    inventory = {Path(item["path"]).resolve() for item in data["outputs"]}
    for command in data["commands"]:
        if command["returncode"] != 0:
            raise ValueError("Failed documented command")
        if oracle.sha(Path(command["executed_argv"][0])) != data["binary"]["sha256"]:
            raise ValueError("Documented command ran different executable")
        for flag in (
            "--data",
            "--metadata",
            "--output",
            "--predictions",
            "--portable-model",
        ):
            value = option(command["executed_argv"], flag)
            if value and Path(value).resolve() not in inventory | set(by_path):
                raise ValueError("Unbound documented command input/output: " + value)
    return data


def load_reference(work, sealed):
    # Execute the unchanged original module only in a new staging tree. Its
    # import-time directory creation and domain main() cannot touch the audit.
    (work / "models").mkdir(parents=True)
    for name in ("models_audit", "models_domain", "models_screen"):
        sealed.copy("scripts/" + name + ".py", work / "scripts" / (name + ".py"))
    sys.path.insert(0, str(work / "scripts"))
    spec = importlib.util.spec_from_file_location(
        "models_audit", work / "scripts/models_audit.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["models_audit"] = module
    spec.loader.exec_module(module)
    return module


def no_sut(*args, **kwargs):
    raise RuntimeError(
        "Independent reference worker attempted SUT/subprocess execution"
    )


def independent_fit(
    core, check, name, data_path, metadata_path, report_path, portable_path=None
):
    np = core.np
    raw = np.fromfile(data_path, dtype="<f4").reshape(-1, 16)
    metadata = rows(metadata_path)
    if len(raw) != len(metadata):
        raise ValueError("Packed records/metadata cardinality mismatch")
    train = [
        i for i, row in enumerate(metadata) if row["Dataset_Split"].lower() == "train"
    ]
    records = raw[train]
    x, y = core.expand(records), records[:, 6].astype(float)
    report = read(report_path)
    selected = report["selected_feature_indices"]
    check.equal(name + "/training_count", report["training_count"], len(train))
    check.equal(name + "/selected_features", selected, core.select(x, y))
    beta = core.linear(x[:, selected], y)
    fitted = core.predict(x[:, selected], beta)
    for i, column in enumerate([0, *selected]):
        check.close(
            f"{name}/coefficient/{column}",
            report["weights"][column],
            beta[i],
            atol=float(2 * abs(np.spacing(np.float32(beta[i])))),
            rtol=0,
        )
    groups = [metadata[i]["Ligand_Group"] for i in train]
    predictions = {
        "training": fitted,
        "fixed_feature_loo": core.cv(x, y, selected),
        "fixed_feature_group_loo": core.cv(x, y, selected, groups=groups),
    }
    for key, prediction in predictions.items():
        for metric, expected in core.metrics(y, prediction).items():
            if metric == "n":
                check.equal(f"{name}/{key}/count", report[key]["count"], expected)
            else:
                check.close(
                    f"{name}/{key}/{metric}",
                    report[key][metric],
                    expected,
                    atol=8 * float(np.finfo(np.float32).eps),
                    rtol=8 * float(np.finfo(np.float32).eps),
                )
    for kind, grid in (
        ("ridge", [1e-4, 1e-3, 1e-2, 0.1, 1, 10]),
        ("lasso", [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5]),
    ):
        prediction = core.cv(x, y, selected, kind, grid)
        for metric, expected in core.metrics(y, prediction).items():
            if metric != "n":
                check.close(
                    f"{name}/{kind}_fixed_feature_nested_alpha/{metric}",
                    report[kind + "_baseline"]["nested_loo"][metric],
                    expected,
                    atol=8 * float(np.finfo(np.float32).eps),
                    rtol=8 * float(np.finfo(np.float32).eps),
                )
    z = (x[:, selected] - x[:, selected].mean(0)) / x[:, selected].std(0)
    design = np.column_stack([np.ones(len(y)), z])
    inv = np.linalg.inv(design.T @ design)
    geometry = report["training_geometry"]
    for i, values in enumerate(inv):
        for j, expected in enumerate(values):
            check.close(
                f"{name}/xtx_inverse/{i}/{j}", geometry["xtx_inverse"][i][j], expected
            )
    rse = float(np.sqrt(np.sum((y - fitted) ** 2) / (len(y) - len(selected) - 1)))
    check.close(
        name + "/residual_standard_error",
        geometry["residual_standard_error"],
        rse,
        atol=8 * float(np.finfo(np.float32).eps),
    )
    if portable_path:
        portable = read(portable_path)
        check.equal(name + "/portable_weights", portable["weights"], report["weights"])
        bootstrap = portable.get("uncertainty")
        if bootstrap and bootstrap.get("replicates"):
            intervals = np.quantile(
                np.asarray(bootstrap["replicates"]), [0.025, 0.975], axis=0
            ).T
            for i, pair in enumerate(intervals):
                for side, expected in zip(("lower_95", "upper_95"), pair, strict=True):
                    check.close(
                        f"{name}/bootstrap/{i}/{side}",
                        report["coefficient_intervals"][i][side],
                        expected,
                    )
    return {
        "name": name,
        "training_rows": train,
        "selected": selected,
        "independent_coefficients": beta.tolist(),
        "independent_fitted": fitted.tolist(),
        "full_pipeline_loo_reference_only": core.metrics(
            y, core.cv(x, y, selected, reselect=True)
        ),
        "full_pipeline_reference_scope": (
            "Reference-only train-mean fallback on empty selection; not a "
            "completed native workflow."
        ),
    }


def has_number(value):
    if isinstance(value, dict):
        return any(has_number(item) for item in value.values())
    if isinstance(value, list):
        return any(has_number(item) for item in value)
    return isinstance(value, (float, int)) and not isinstance(value, bool)


def exact_historical_numbers(check, name, actual, historical, path=""):
    """No numerical allowance for unchanged valid historical model arithmetic."""
    if isinstance(historical, dict):
        for key, value in historical.items():
            # These two calculations intentionally changed (SI kinetics and
            # ordinary covariance reconstructed from stored training points).
            if key in {"predicted_ee_percent", "mahalanobis_distance"}:
                continue
            if key not in actual:
                if has_number(value):
                    check.equal(
                        name + "/historical_numeric_field" + path + "/" + key,
                        False,
                        True,
                    )
                continue  # textual/schema changes are outside the numeric check
            exact_historical_numbers(check, name, actual[key], value, path + "/" + key)
    elif isinstance(historical, list) and has_number(historical):
        check.equal(name + "/historical_length" + path, len(actual), len(historical))
        for i, (left, right) in enumerate(zip(actual, historical, strict=False)):
            exact_historical_numbers(check, name, left, right, path + "/" + str(i))
    elif isinstance(historical, (int, float)) and not isinstance(historical, bool):
        check.equal(name + "/historical_exact" + path, actual, historical)


def input_features(core, library, packed=None):
    np = core.np
    source = rows(library)
    if packed:
        raw = np.fromfile(packed, dtype="<f4").reshape(-1, 16)
        if len(raw) != len(source):
            raise ValueError("Corrected screen library/packed rows mismatch")
        return {
            row["Reaction_ID"]: core.expand(raw[i : i + 1])[0]
            for i, row in enumerate(source)
        }
    answer = {}
    for row_number, row in enumerate(source, start=2):
        values = [
            row.get(key, "")
            for key in (
                "sterimol_l",
                "sterimol_b1",
                "sterimol_b5",
                "nbo_charge",
                "ir_frequency",
            )
        ]
        # A missing required descriptor is separately checked against exclusions.
        raw = np.zeros((1, 16), dtype=np.float32)
        for i, value in enumerate(values):
            raw[0, i] = float(value) if value else 0.0
        answer[row.get("Reaction_ID", f"row_{row_number}")] = core.expand(raw)[0]
    return answer


def independent_screen(core, check, name, model_path, library, output, packed=None):
    np, stats = core.np, core.stats
    model, report = read(model_path), read(output)
    check.equal(
        name + "/model_input_hash",
        report["provenance"]["model_sha256"],
        oracle.sha(model_path),
    )
    check.equal(
        name + "/library_input_hash",
        report["provenance"]["library_sha256"],
        oracle.sha(library),
    )
    features = input_features(core, library, packed)
    geometry = model["training_geometry"]
    selected = model["selected_feature_indices"]
    train = np.asarray(geometry["standardized_training_points"])
    design = np.column_stack([np.ones(len(train)), train])
    inv = np.linalg.inv(design.T @ design)
    covariance = np.atleast_2d(np.cov(train, rowvar=False))
    covariance_inv = np.linalg.inv(covariance)
    weights = np.asarray(model["weights"], dtype=np.float32).astype(float)
    refs_by_id = {}
    for hit in report["hits"]:
        ident = hit["ligand"]
        row = features[ident]
        prediction = float(row @ weights)
        z = (row[selected] - geometry["means"]) / geometry["scales"]
        h = float(np.r_[1, z] @ inv @ np.r_[1, z])
        t = stats.t.isf(0.025, geometry["observations"] - geometry["parameters"])
        half = float(t * geometry["residual_standard_error"] * np.sqrt(1 + h))
        distance = float(np.sqrt(((train - z) ** 2).sum(1)).min())
        mahalanobis = float(np.sqrt(z @ covariance_inv @ z))
        expected = dict(
            predicted_ddg_kcal_mol=prediction,
            leverage=h,
            prediction_interval_low=prediction - half,
            prediction_interval_high=prediction + half,
            nearest_training_distance=distance,
            mahalanobis_distance=mahalanobis,
        )
        for key, value in expected.items():
            check.close(f"{name}/{ident}/{key}", hit[key], value)
        # Source descriptor checks are independent of the screen's reported
        # descriptor values. Corrected packed records are same-build SUT values;
        # this is the training-versus-screen input contract, not a geometry oracle.
        names = {
            "sterimol_l": 1,
            "sterimol_b1": 2,
            "sterimol_b5": 3,
            "nbo_charge": 4,
            "ir_frequency": 7,
        }
        for descriptor in hit["descriptors"]:
            check.close(
                f"{name}/{ident}/descriptor/{descriptor['name']}",
                descriptor["value"],
                row[names[descriptor["name"]]],
                atol=0,
                rtol=0,
            )
        bootstrap = model.get("uncertainty")
        if bootstrap and bootstrap.get("replicates"):
            samples = (
                np.asarray(bootstrap["replicates"]) @ row[bootstrap["column_indices"]]
            )
            interval = np.quantile(samples, [0.025, 0.975])
            for key, value in zip(("lower", "upper"), interval, strict=True):
                check.close(
                    f"{name}/{ident}/bootstrap/{key}", hit["uncertainty"][key], value
                )
            expected["bootstrap_interval"] = interval.tolist()
        if name.endswith("marginal_ci_counterexample"):
            envelope = (hit["coefficient_band_low"], hit["coefficient_band_high"])
            coverage = float(
                np.mean((samples >= envelope[0]) & (samples <= envelope[1]))
            )
            check.close(
                name + "/retained_envelope_empirical_coverage",
                coverage,
                0.904,
                atol=1e-12,
                rtol=0,
            )
            check.equal(
                name + "/no_universal_coverage_claim",
                "without a universal coverage guarantee" in report["uncertainty_note"],
                True,
            )
        refs_by_id[ident] = expected
    check.equal(name + "/nonempty_hits", bool(report["hits"]), True)
    check.equal(name + "/returned_count", report["returned"], len(report["hits"]))
    if not report["diversity_applied"]:

        def score(ident):
            value = refs_by_id[ident]["predicted_ddg_kcal_mol"]
            if report["ranking_order"] == "ascending":
                return value, ident
            if report["ranking_order"] == "magnitude_descending":
                return -abs(value), ident
            return -value, ident

        check.equal(
            name + "/reported_order",
            [h["ligand"] for h in report["hits"]],
            sorted(refs_by_id, key=score),
        )
    return {
        "name": name,
        "input_model": oracle.file_record(model_path),
        "library": oracle.file_record(library),
        "independent": refs_by_id,
        "diversity_order_scope": (
            "All returned scalar predictions/intervals checked; greedy "
            "diversity order is covered by retained CLI and focused tests, "
            "not independently reimplemented here."
        ),
    }


def domain_references(core, work, sealed, lane, check):
    source = sealed.check("models/inputs/domain_requests.jsonl")
    if oracle.sha(source) != lane["requests"]["sha256"] or lane["rows"] != 33:
        raise ValueError("Original 33 model observer requests changed")
    shutil.copyfile(source, work / "models/inputs/domain_requests.jsonl")
    shutil.copyfile(
        refs.verify(lane["stdout"]), work / "models/raw/domain_observer.jsonl"
    )
    spec = importlib.util.spec_from_file_location(
        "models_domain", work / "scripts/models_domain.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main()
    request_by_id = {row["id"]: row for row in oracle.requests(source)}
    adjusted_rows = []
    for row in read(work / "models/results/domain_comparison.json"):
        ident, actual, expected = row["id"], row["stericx"], row["independent"]
        if "leverage" in expected:
            np = core.np
            request = request_by_id[ident]
            geometry = request["geometry"]
            features = np.asarray(request["features"], dtype=np.float32).astype(float)
            z = (features[geometry["feature_indices"]] - geometry["means"]) / geometry[
                "scales"
            ]
            points = np.asarray(geometry["standardized_training_points"])
            design = np.column_stack([np.ones(len(points)), points])
            h = float(np.r_[1, z] @ np.linalg.inv(design.T @ design) @ np.r_[1, z])
            mahalanobis = float(
                np.sqrt(z @ np.linalg.inv(np.cov(points, rowvar=False)) @ z)
            )
            distance = float(np.sqrt(((points - z) ** 2).sum(1)).min())
            t = float(
                core.stats.t.isf(
                    0.025, geometry["observations"] - geometry["parameters"]
                )
            )
            s = geometry["residual_standard_error"]
            prediction = request["prediction"]
            adjusted = dict(
                leverage=h,
                mahalanobis=mahalanobis,
                nearest_distance=distance,
                pi=[
                    prediction - t * s * math.sqrt(1 + h),
                    prediction + t * s * math.sqrt(1 + h),
                ],
                ci=[
                    prediction - t * s * math.sqrt(h),
                    prediction + t * s * math.sqrt(h),
                ],
            )
            adjusted_rows.append(
                dict(
                    id=ident,
                    original_reference=expected,
                    f32_input_reference=adjusted,
                    features_received_by_f32_api=features.tolist(),
                    reason=(
                        "Same independent covariance/hat/Student-t equations, "
                        "applying the public feature input f32 conversion "
                        "before arithmetic."
                    ),
                )
            )
            expected = adjusted
        if "t_quantile" in expected:
            check.close(
                "domain/" + ident,
                actual["t_quantile"],
                expected["t_quantile"],
                atol=2e-10,
                rtol=2e-12,
            )
        elif "matrix_rank" in expected:
            check.equal(
                "domain/" + ident + "/ordinary_covariance_unavailable",
                "Err" in actual["mahalanobis"],
                True,
            )
            check.equal(
                "domain/" + ident + "/covariance_reason",
                "singular" in actual["mahalanobis"].get("Err", ""),
                True,
            )
        else:
            for key in ("leverage", "mahalanobis"):
                value = (
                    actual[key].get("Ok")
                    if isinstance(actual[key], dict)
                    else actual[key]
                )
                check.close("domain/" + ident + "/" + key, value, expected[key])
            for label, key in (
                ("pi", "prediction_interval"),
                ("ci", "confidence_interval"),
            ):
                for i in range(2):
                    check.close(
                        f"domain/{ident}/{label}/{i}",
                        actual[key][i],
                        expected[label][i],
                    )
            check.close(
                "domain/" + ident + "/nearest",
                actual["assessment"]["nearest_training_distance"],
                expected["nearest_distance"],
            )
    save(work / "models/results/domain_f32_input_comparison.json", adjusted_rows)


def capture_observer_extra(root, binary, original_requests):
    near_path = (
        REPO / "docs/scientific_remediation/models/near_unit_alpha/requests.jsonl"
    )
    near_record = oracle.file_record(near_path)
    if (
        near_record["sha256"]
        != "d9be53f53bb0b3a2181a5584f45dac694b318bf47f89a6340760d078bffae9bf"
    ):
        raise ValueError(
            "Near-unit-alpha requests differ from preserved before witness"
        )
    requests = list(oracle.requests(near_path))
    if len(requests) != 10:
        raise ValueError("Expected all ten witnessed near-unit-alpha cases")
    for name, df, alpha in (
        ("df_zero", 0, 0.05),
        ("df_negative", -1, 0.05),
        ("alpha_zero", 1, 0),
        ("alpha_one", 1, 1),
        ("alpha_negative", 1, -0.1),
        ("alpha_above_one", 1, 1.1),
    ):
        requests.append(dict(id=name, op="model", df=df, alpha=alpha))
    base = copy.deepcopy(next(oracle.requests(original_requests)))
    for name in ("feature_index", "inverse_shape", "zero_scale", "point_shape"):
        case = copy.deepcopy(base)
        case["id"] = "malformed_" + name
        g = case["geometry"]
        if name == "feature_index":
            g["feature_indices"][0] = 8
        elif name == "inverse_shape":
            g["xtx_inverse"][0] = []
        elif name == "zero_scale":
            g["scales"][0] = 0
        else:
            g["standardized_training_points"][0] = []
        requests.append(case)
    root.mkdir()
    request_path = root / "requests.jsonl"
    request_path.write_text("".join(json.dumps(row) + "\n" for row in requests))
    with (
        request_path.open("rb") as stdin,
        (root / "stdout").open("xb") as stdout,
        (root / "stderr").open("xb") as stderr,
    ):
        run = subprocess.run(
            [str(binary)],
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            cwd=REPO,
            env=dict(os.environ, RAYON_NUM_THREADS="1", LC_ALL="C"),
            timeout=120,
        )
    save(
        root / "command.json",
        dict(
            binary=oracle.file_record(binary),
            returncode=run.returncode,
            original_near_unit_requests=oracle.file_record(near_path),
        ),
    )
    if run.returncode:
        raise ValueError("Additional model observer failed; raw logs retained")
    answers = list(oracle.requests(root / "stdout"))
    if len(answers) != len(requests):
        raise ValueError("Incomplete additional model observer outputs")
    for request, answer in zip(requests, answers, strict=True):
        if answer.get("id") != request["id"] and set(answer) != {"error"}:
            raise ValueError("Reordered or unidentified additional model output")
        # A caught public-API panic is retained as a failed comparison below.
        # It is never admitted as a successful invalid-input result.
    refs.verify(near_record)
    return requests, answers


def extra_references(requests, answers, check):
    # Analytic Cauchy and central Student-t density, stable as alpha approaches 1.
    # For these other df and masses <=1e-8, omitted relative cubic terms <2e-16.
    for request, actual in zip(requests, answers, strict=True):
        name = "additional/" + request["id"]
        check.equal(name + "/identified_response", actual.get("id"), request["id"])
        check.equal(name + "/no_observer_panic", actual.get("error"), None)
        if "error" in actual:
            continue
        if request["id"].startswith("malformed_"):
            check.equal(
                name + "/assessment_unknown", actual["assessment"]["verdict"], "unknown"
            )
            check.equal(
                name + "/explicit_unavailability",
                bool(actual["assessment"].get("unavailable")),
                True,
            )
            check.equal(
                name + "/nearest_unavailable",
                actual["assessment"]["nearest_training_distance"],
                None,
            )
            check.equal(
                name + "/prediction_interval_unavailable",
                actual["prediction_interval"],
                None,
            )
            check.equal(
                name + "/confidence_interval_unavailable",
                actual["confidence_interval"],
                None,
            )
            check.equal(name + "/leverage_unavailable", actual["leverage"], None)
            check.equal(
                name + "/covariance_error", "Err" in actual["mahalanobis"], True
            )
            continue
        df, alpha = request["df"], request["alpha"]
        if not (df > 0 and 0 < alpha < 1):
            check.equal(
                name + "/invalid_quantile_unavailable", actual["t_quantile"], None
            )
            continue
        mass = 1 - alpha
        if df == 1:
            expected = math.tan(math.pi * mass / 2)
        else:
            density = {
                2: 1 / (2 * math.sqrt(2)),
                3: 2 / (math.pi * math.sqrt(3)),
                4: 3 / 8,
                10: 315 / (256 * math.sqrt(10)),
            }[df]
            expected = mass / (2 * density)
        check.close(
            name + "/central_quantile",
            actual["t_quantile"],
            expected,
            atol=0,
            rtol=2e-12,
        )


def capture_contracts(root, binary, model_path, library):
    root.mkdir()
    model, source = read(model_path), rows(library)
    row = next(row for row in source if ";" in row["Conformer_XYZ_Paths"])
    result = []
    count = len(row["Conformer_XYZ_Paths"].split(";"))
    negative_weights = ";".join(["-1", *(["1"] * (count - 1))])
    nonfinite_weights = ";".join(["NaN", *(["1"] * (count - 1))])
    variants = [
        ("unknown", None, None),
        ("record_values", "supplied_record_values", None),
        ("single_geometry", "single_geometry", None),
        ("missing_weights", "supplied_weight_mean", ""),
        ("negative_weights", "supplied_weight_mean", negative_weights),
        (
            "zero_weights",
            "supplied_weight_mean",
            ";".join("0" for _ in row["Conformer_XYZ_Paths"].split(";")),
        ),
        ("mismatched_weights", "supplied_weight_mean", "1"),
        ("nonfinite_weights", "supplied_weight_mean", nonfinite_weights),
    ]
    for name, aggregation, weights in variants:
        out = root / name
        out.mkdir()
        changed = copy.deepcopy(model)
        if aggregation is None:
            changed["schema_version"] = 2
            changed.pop("descriptor_aggregation", None)
        else:
            changed["descriptor_aggregation"] = aggregation
        save(out / "model.json", changed)
        data = dict(row)
        if weights is not None:
            data["Conformer_Boltzmann_Weights"] = weights
        with (out / "library.csv").open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data))
            writer.writeheader()
            writer.writerow(data)
        argv = [
            str(binary),
            "screen",
            str(out / "model.json"),
            str(out / "library.csv"),
            "--format",
            "json",
        ]
        with (
            (out / "stdout").open("xb") as stdout,
            (out / "stderr").open("xb") as stderr,
        ):
            proc = subprocess.run(
                argv, cwd=REPO, stdout=stdout, stderr=stderr, timeout=120
            )
        item = dict(
            name=name,
            argv=argv,
            returncode=proc.returncode,
            stderr=(out / "stderr").read_text(),
            source_model=oracle.file_record(model_path),
            source_library=oracle.file_record(library),
        )
        save(out / "command.json", item)
        result.append(item)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("build", "observations", "cli_run", "demo", "input_bundle", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    args = parser.parse_args()
    args = argparse.Namespace(
        **{name: value.resolve() for name, value in vars(args).items()}
    )
    output = replay.new_output(args.output)
    try:
        run(args, output)
    except Exception as exc:
        save(output / "failure.json", dict(type=type(exc).__name__, detail=str(exc)))
        raise


def run(args, output):
    build_path = args.build / "manifest.json"
    observed_path = args.observations / "manifest.json"
    built, _ = replay.verify_build(build_path)
    if built["audit_manifest"]["sha256"] != AUDIT_SHA:
        raise ValueError("Unexpected historical audit seal")
    observed, observation_build = replay.verify_observations(observed_path)
    if observation_build != built:
        raise ValueError("Observer and CLI model build differ")
    cases = verify_cli(args.cli_run, build_path, built)
    demo_path = args.demo / "receipt.json"
    demo = verify_demo(demo_path, build_path, built, args.input_bundle)
    bindings = [
        oracle.file_record(p)
        for p in (
            build_path,
            observed_path,
            args.cli_run / "manifest.json",
            args.cli_run.parent.parent / "manifest.json",
            demo_path,
            args.input_bundle / "manifest.json",
            Path(__file__),
            Path(capture.__file__),
            Path(capture.cli.__file__),
            Path(replay.__file__),
            Path(oracle.__file__),
            Path(refs.__file__),
            REPO / "docs/scientific_remediation/models/near_unit_alpha/requests.jsonl",
            REPO
            / "docs/scientific_remediation/models/near_unit_alpha/before_manifest.json",
        )
    ]
    for path in (
        Path(__file__),
        Path(capture.__file__),
        Path(capture.cli.__file__),
        Path(replay.__file__),
        Path(oracle.__file__),
        Path(refs.__file__),
    ):
        destination = output / "helpers" / path.name
        destination.parent.mkdir(exist_ok=True)
        shutil.copyfile(path, destination)
    sealed = refs.SealedInputs(replay.AUDIT)
    work = output / "reference"
    core = load_reference(work, sealed)
    tolerance_evidence = [
        oracle.file_record(sealed.check(path))
        for path in (
            "models/REPORT.md",
            "models/results/domain_comparison.json",
            "models/results/synthetic_linear_math.json",
            "models/results/native_ni_hda_math.json",
        )
    ]
    tolerance_policy = {
        "unchanged_historical_arithmetic": (
            "Exact numeric equality and hit identity/order; excludes "
            "intentionally corrected ee and covariance distance."
        ),
        "independent_f64": {
            "atol": 1e-12,
            "rtol": 2e-12,
            "reason": (
                "Small well-conditioned matrices: allowance for finite f64 "
                "solve/order differences; original report residuals reach "
                "1.43e-13. Raw residuals retained."
            ),
        },
        "coefficient_f32": (
            "Absolute bound two ULPs at independent reference coefficient, "
            "accounting for stored f32 plus shortest-decimal output; zero "
            "relative allowance."
        ),
        "fit_metrics_f32": (
            "Eight f32 eps absolute and relative (9.5367431640625e-7); "
            "independent SVD folds versus stored f32 fold coefficients and "
            "small solve floor. Original valid model arithmetic "
            "additionally exact, so this cannot admit a historical "
            "numerical regression."
        ),
        "student_t_original": {
            "atol": 2e-10,
            "rtol": 2e-12,
            "reason": (
                "Original independent report records max main-quantile error "
                "1.73e-10; analytic extreme/near-unit Cauchy tests separately "
                "require zero absolute,2e-12 relative."
            ),
        },
        "near_unit_alpha": {
            "atol": 0,
            "rtol": 2e-12,
            "reason": (
                "Analytic Cauchy / Student-t central-density equations; same "
                "focused remediation test bound."
            ),
        },
        "provenance": tolerance_evidence,
    }
    save(output / "tolerance_policy.json", tolerance_policy)
    staged_reference_records = [
        oracle.file_record(p) for p in sorted((work / "scripts").glob("*.py"))
    ]
    case_copies = output / "cli"
    for case in cases:
        shutil.copytree(args.cli_run / case["id"], case_copies / case["id"])
    shutil.copytree(args.demo, output / "documented_demo")
    check = Checks()
    source_rows = rows(args.input_bundle / "reactions.csv")
    for name in ("screening", "model_format"):
        packed = core.np.fromfile(
            args.demo / name / "reactions.sigpack", dtype="<f4"
        ).reshape(-1, 16)
        check.equal(
            "corrected/" + name + "/record_count", len(packed), len(source_rows)
        )
        for i, source_row in enumerate(source_rows):
            for column, index in (
                ("NBO_Charge", 3),
                ("IR_Frequency", 4),
                ("Temp_K", 5),
                ("Exp_ddG_kcal_mol", 6),
            ):
                check.close(
                    f"corrected/{name}/record/{i}/{column}",
                    float(packed[i, index]),
                    float(core.np.float32(source_row[column])),
                    atol=0,
                    rtol=0,
                )
        model = read(args.demo / name / "model.json")
        check.equal(
            "corrected/" + name + "/aggregation",
            model["descriptor_aggregation"],
            "supplied_weight_mean",
        )
        check.equal(
            "corrected/" + name + "/response_temperature",
            model["inference"]["response"]["temperature_k"],
            353.15,
        )
    extra_requests, extra_answers = capture_observer_extra(
        output / "additional_observer",
        refs.verify(built["binaries"]["observer"]),
        refs.verify(observed["lanes"]["model_domain"]["requests"]),
    )
    contracts = capture_contracts(
        output / "input_contracts",
        refs.verify(built["binaries"]["native"]),
        args.demo / "model_format/model.json",
        args.input_bundle / "reactions.csv",
    )
    save(
        output / "started.json",
        dict(
            bindings=bindings,
            staged_reference_programs=staged_reference_records,
            source_build=built["binaries"],
            model_cases=[c["id"] for c in cases],
        ),
    )
    # All SUT execution is complete. Disallow accidental fallback execution by
    # any imported reference function, even if a required staged output is absent.
    original_run, original_popen = subprocess.run, subprocess.Popen
    subprocess.run = subprocess.Popen = core.run = no_sut
    results = []
    try:
        domain_references(core, work, sealed, observed["lanes"]["model_domain"], check)
        extra_references(extra_requests, extra_answers, check)
        for item in contracts:
            check.equal(
                "contract/" + item["name"] + "/rejected", item["returncode"] != 0, True
            )
            check.equal(
                "contract/" + item["name"] + "/reason",
                {
                    "unknown": "descriptor_aggregation=unknown",
                    "record_values": "descriptor_aggregation=supplied_record_values",
                    "single_geometry": "multiple conformers",
                    "missing_weights": "requires explicit Conformer_Boltzmann_Weights",
                    "negative_weights": "finite and non-negative",
                    "nonfinite_weights": "finite and non-negative",
                    "zero_weights": "positive sum",
                    "mismatched_weights": "equal non-zero lengths",
                }[item["name"]]
                in item["stderr"],
                True,
            )
        failed_outer = {2, 4, 5, 6, 7, 8, 9}
        for case in cases:
            ident, argv = case["id"], case["argv"]
            root = case_copies / ident
            status = read(root / "result.json")["returncode"]
            stderr = (root / "stderr").read_text()
            if argv[0] == "evaluate":
                check.equal(ident + "/reject_invalid", status != 0, True)
                reason = (
                    "must both be finite"
                    if not ident.endswith("finite_mismatch")
                    else "does not match"
                )
                check.equal(ident + "/error_reason", reason in stderr, True)
            elif ident == "models_objective_unspecified":
                check.equal(ident + "/reject_unspecified", status != 0, True)
                check.equal(
                    ident + "/reason", "does not record which direction" in stderr, True
                )
            elif (
                ident.startswith("models_native_outer_")
                and int(ident.rsplit("_", 1)[1]) in failed_outer
            ):
                check.equal(ident + "/retained_failure", status != 0, True)
                check.equal(
                    ident + "/reason",
                    "no non-constant descriptor improved" in stderr,
                    True,
                )
            else:
                check.equal(ident + "/success", status, 0)
                if status:
                    continue
                historical_raw = "models/raw/" + ident.removeprefix("models_")
                if argv[0] == "fit":
                    historical_report = sealed.check(historical_raw + "/report.json")
                    exact_historical_numbers(
                        check,
                        ident,
                        read(root / "artifacts/report.json"),
                        read(historical_report),
                    )
                    portable = root / "artifacts/portable.json"
                    results.append(
                        independent_fit(
                            core,
                            check,
                            ident,
                            Path(option(argv, "--data")),
                            Path(option(argv, "--metadata")),
                            root / "artifacts/report.json",
                            portable if portable.exists() else None,
                        )
                    )
                elif argv[0] == "screen":
                    historical_screen = read(sealed.check(historical_raw + "/stdout"))
                    current_screen = read(root / "stdout")
                    check.equal(
                        ident + "/historical_hit_identity_order",
                        [h["ligand"] for h in current_screen["hits"]],
                        [h["ligand"] for h in historical_screen["hits"]],
                    )
                    exact_historical_numbers(
                        check, ident, current_screen, historical_screen
                    )
                    results.append(
                        independent_screen(
                            core,
                            check,
                            ident,
                            Path(argv[1]),
                            Path(argv[2]),
                            root / "stdout",
                        )
                    )
        for command in demo["commands"]:
            argv = command["executed_argv"][1:]
            if argv[0] == "fit":
                results.append(
                    independent_fit(
                        core,
                        check,
                        "corrected/" + Path(option(argv, "--output")).parent.name,
                        Path(option(argv, "--data")),
                        Path(option(argv, "--metadata")),
                        Path(option(argv, "--output")),
                        Path(option(argv, "--portable-model")),
                    )
                )
            elif argv[0] == "screen" and option(argv, "--format") == "json":
                name = (
                    "corrected/"
                    + Path(argv[1]).parent.name
                    + "/"
                    + Path(command["redirected_stdout"]["path"]).stem
                )
                results.append(
                    independent_screen(
                        core,
                        check,
                        name,
                        Path(argv[1]),
                        Path(option(argv, "--library")),
                        Path(command["stdout"]["path"]),
                        Path(argv[1]).parent / "reactions.sigpack",
                    )
                )
    finally:
        subprocess.run, subprocess.Popen = original_run, original_popen
    save(output / "independent_results.json", results)
    save(output / "comparisons.json", check.items)
    for binding in bindings + staged_reference_records + list(sealed.used.values()):
        refs.verify(binding)
    replay.verify_build(build_path)
    verify_cli(args.cli_run, build_path, built)
    verify_demo(demo_path, build_path, built, args.input_bundle)
    # Observation identity and all lane bytes were checked at admission; rehash
    # every stream at completion without a second expensive full schema parse.
    verify_records(observed)
    failed = [row for row in check.items if not row["passed"]]
    packages = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "pandas", "scikit-learn")
    }
    replay.manifest(
        output / "manifest.json",
        dict(
            kind="corrected_model_independent_reference_recheck",
            complete=True,
            numerical_contract_checks_passed=not failed,
            scientific_global_pass=False,
            original_cli_cases=33,
            original_model_observer_rows=33,
            additional_model_observer_rows=len(extra_requests),
            additional_input_contract_cases=len(contracts),
            comparisons=len(check.items),
            tolerance_policy=tolerance_policy,
            failed_comparisons=failed,
            bindings=bindings,
            sealed_reference_inputs=list(sealed.used.values()),
            staged_reference_programs=staged_reference_records,
            environment=dict(
                python=sys.version,
                executable=oracle.file_record(Path(sys.executable).resolve()),
                packages=packages,
                platform=platform.platform(),
            ),
            limitations=[
                (
                    "Numerical contracts do not establish experimental validity or "
                    "calibration."
                ),
                (
                    "Original seven failed full-pipeline folds remain failed; no "
                    "native score is invented."
                ),
                (
                    "Bootstrap replicate generation and permutation RNG are not "
                    "independently reimplemented; propagation and quantiles are "
                    "checked from actual replicates."
                ),
                (
                    "All screen scalar predictions and intervals use the model "
                    "actually passed to that command. Legacy input models are not "
                    "substituted for corrected fitted models."
                ),
                (
                    "Corrected screen descriptors are checked against same-build "
                    "packed training inputs. Independent geometry accuracy is a "
                    "separate reference gate."
                ),
                (
                    "Greedy diversity order remains covered by focused tests; it is "
                    "not independently reimplemented in this adapter."
                ),
            ],
            files=[
                oracle.file_record(path)
                for path in sorted(output.rglob("*"))
                if path.is_file()
            ],
        ),
    )
    print(
        json.dumps(
            dict(
                output=str(output),
                comparisons=len(check.items),
                failed=len(failed),
            )
        )
    )
    if failed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
