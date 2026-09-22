#!/usr/bin/env python3
"""Observe LIVE Rust source without modifying the sealed scientific audit.

This is an exact regression oracle, NOT a claim that the audited science passes.
Known failures, classifications, and independent residual evidence remain locked.
No observed maximum error is promoted to a tolerance for other inputs.

Commands (every destination must be new):
  prepare --output .stericx/scientific-oracle
  build --oracle .stericx/scientific-oracle --output /tmp/baseline-observer
  observe --oracle .stericx/scientific-oracle --build /tmp/baseline-observer \
      --output /tmp/baseline-observations --require-frozen-fidelity
  compare --baseline /tmp/baseline-observations \
      --candidate /tmp/candidate-observations --output /tmp/exact-comparison.json

Build links the selected source tree and appends the unchanged visibility adapter
to a byte-identical copy of its current BV source. The copied implementation is
the SUT, never an independent reference. Observation can select --lanes for a
bounded check; comparison reports that scope and never calls it full coverage.
CLI commands and Python preparation require separate oracles.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import struct
import subprocess
import sys
import time
import tomllib
from collections import Counter
from datetime import UTC, datetime
from itertools import zip_longest
from pathlib import Path

REPO = next(
    parent
    for parent in Path(__file__).resolve().parents
    if (parent / "Cargo.toml").is_file()
    and (parent / "docs/scientific_accuracy_audit").is_dir()
)
DEFAULT_AUDIT = REPO / "docs/scientific_accuracy_audit"
PURPOSE = (
    "Exact implementation equivalence with preserved documented findings; "
    "not independent scientific validation or a scientific pass."
)
LANES = {
    "geometry": (
        "geometry/requests.jsonl",
        "geometry/sut_outputs.jsonl",
        "geometry/sut_stderr.txt",
        844,
    ),
    "geometry_focused": (
        "geometry/focused/requests.jsonl",
        "geometry/focused/sut_outputs.jsonl",
        "geometry/focused/stderr.txt",
        15,
    ),
    "geometry_alignment": (
        "geometry/alignment/requests.jsonl",
        "geometry/alignment/sut_outputs.jsonl",
        "geometry/alignment/stderr.txt",
        4,
    ),
    "geometry_alignment_rotation": (
        "geometry/alignment/rotation/requests.jsonl",
        "geometry/alignment/rotation/sut_outputs.jsonl",
        "geometry/alignment/rotation/stderr.txt",
        4,
    ),
    "geometry_numerical": (
        "geometry/numerical/requests.jsonl",
        "geometry/numerical/sut_outputs.jsonl",
        "geometry/numerical/stderr.txt",
        11,
    ),
    "geometry_historical": (
        "geometry/historical/requests.jsonl",
        "geometry/historical/sut_outputs.jsonl",
        "geometry/historical/stderr.txt",
        20,
    ),
    "kraken": (
        "kraken/prepared_all/requests.jsonl",
        "kraken/sut/stdout.jsonl",
        "kraken/sut/stderr.txt",
        31721,
    ),
    "kinetics": (
        "kinetics/inputs/kinetics.json",
        "kinetics/frozen_outputs/kinetics.jsonl",
        "kinetics/frozen_outputs/kinetics.stderr",
        195,
    ),
    "aggregation": (
        "kinetics/inputs/aggregation.json",
        "kinetics/frozen_outputs/aggregation.jsonl",
        "kinetics/frozen_outputs/aggregation.stderr",
        11,
    ),
    "model_domain": (
        "models/inputs/domain_requests.jsonl",
        "models/raw/domain_observer.jsonl",
        "models/raw/domain_observer.stderr",
        33,
    ),
}
ALL_LANES = (*LANES, "kraken_all_bins")
SCHEMA_VERSION = 3
# Individually reviewed against the sealed public errors and unchanged adapter:
# donor_neighbor_indices fails before dump construction for exactly these rows.
# KRAKEN:1299:54318 has a public asymmetry error but MUST retain a successful dump.
KRAKEN_DUMP_ERRORS = {
    key: {
        "error": (
            "donor must be trivalent for the quadrant frame, "
            "found 4 bonded substituents"
        )
    }
    for key in (
        "KRAKEN:1281:54241",
        "KRAKEN:1281:54235",
        "KRAKEN:1907:63291",
        "KRAKEN:1907:63292",
        "KRAKEN:1907:63293",
        "KRAKEN:1907:63294",
        "KRAKEN:1907:63295",
        "KRAKEN:1907:63296",
    )
}
BV_FIELDS = (
    "buried_volume percent_buried_volume qvbur_min qvbur_max max_delta_qvbur "
    "ovbur_min ovbur_max near_vbur far_vbur"
).split()
KINETIC_FIELDS = (
    "input_ddg input_temperature rate major_percent minor_percent "
    "r_percent s_percent ee_percent"
).split()
AGGREGATE_FIELDS = (
    "vbur_boltz vbur_min vbur_max vbur_delta qvbur_min_boltz qvbur_max_boltz "
    "max_delta_qvbur_boltz max_delta_qvbur_min max_delta_qvbur_max "
    "max_delta_qvbur_delta "
    "max_delta_qvbur_vburminconf near_vbur_boltz far_vbur_boltz"
).split()


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    ).encode()


def write_new(path: Path, value: object) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def file_record(path: Path) -> dict:
    path = path.resolve(strict=True)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)}


def verify(record: dict) -> Path:
    path = Path(record["path"])
    if (
        not path.is_file()
        or path.stat().st_size != record["bytes"]
        or sha(path) != record["sha256"]
    ):
        raise ValueError(f"Locked file changed or missing: {path}")
    return path


def new_directory(path: Path, audit: Path) -> Path:
    path = path.resolve()
    if path.is_relative_to(audit.resolve()):
        raise ValueError("New evidence must be outside the historical scientific audit")
    path.mkdir(parents=True, exist_ok=False)
    return path


def tool_version(argv: list[str], cwd: Path | None = None) -> dict:
    result = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, check=False)
    return {
        "argv": argv,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def runtime() -> dict:
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
    }


def load_oracle(root: Path) -> tuple[dict, dict]:
    root = root.resolve(strict=True)
    manifest = read(root / "manifest.json")
    if (
        manifest.get("kind") != "live_scientific_oracle_inputs"
        or manifest.get("schema_version") != SCHEMA_VERSION
    ):
        raise ValueError("Not a prepared live scientific oracle")
    current_helper = file_record(Path(__file__))
    verify(manifest["helper"])
    if any(
        manifest["helper"][key] != current_helper[key] for key in ("sha256", "bytes")
    ):
        raise ValueError("Prepared oracle must use this exact current helper")
    for item in manifest["locked_files"]:
        verify(item)
    evidence = read(verify(manifest["evidence_lock"]))
    for item in evidence["files"]:
        verify(item)
    audit = Path(manifest["audit"]).resolve(strict=True)
    sealed = {
        item["path"]: item
        for item in read(verify(evidence["sealed_manifest"]))["files"]
    }
    if evidence["sealed_manifest"] != file_record(audit / "manifest_final.json"):
        raise ValueError("Independent evidence uses another sealed inventory")
    required_evidence = [
        file_record(audit / name) for name in evidence_paths(audit, sealed)
    ]
    if evidence["files"] != required_evidence:
        raise ValueError("Independent evidence inventory is incomplete or substituted")
    for item in required_evidence:
        relative = str(Path(item["path"]).relative_to(audit))
        if (item["sha256"], item["bytes"]) != (
            sealed[relative]["sha256"],
            sealed[relative]["bytes"],
        ):
            raise ValueError(f"Independent evidence is not sealed: {relative}")
    if set(manifest["lanes"]) != set(ALL_LANES):
        raise ValueError("Prepared oracle lane contract is incomplete")
    for name, lane in manifest["lanes"].items():
        input_name, output_name, stderr_name, count = LANES[
            "kraken" if name == "kraken_all_bins" else name
        ]
        for key, relative in (
            ("source", input_name),
            (
                "public_expected_stdout"
                if name == "kraken_all_bins"
                else "expected_stdout",
                output_name,
            ),
            (
                "public_expected_stderr"
                if name == "kraken_all_bins"
                else "expected_stderr",
                stderr_name,
            ),
        ):
            expected = file_record(audit / relative)
            if lane[key] != expected or (expected["sha256"], expected["bytes"]) != (
                sealed[relative]["sha256"],
                sealed[relative]["bytes"],
            ):
                raise ValueError(f"Oracle lane not tied to sealed inputs: {name}/{key}")
        if lane["requests"] != file_record(root / "requests" / f"{name}.jsonl"):
            raise ValueError(f"Unverified request identity: {name}")
        original = (
            json.loads(Path(lane["source"]["path"]).read_text())
            if input_name.endswith(".json")
            else requests(Path(lane["source"]["path"]))
        )
        actual = requests(Path(lane["requests"]["path"]))
        for before, after in zip_longest(original, actual):
            if name == "kraken_all_bins" and before is not None:
                before.update(dump=True, include_points=False)
            if before is None or after is None or before != after:
                raise ValueError(f"Prepared request derivation changed: {name}")
        if lane["expected_rows"] != count or lane["expected_returncode"] != 0:
            raise ValueError(f"Prepared process/row contract changed: {name}")
        expected_errors = (
            KRAKEN_DUMP_ERRORS
            if name == "kraken_all_bins"
            else frozen_dump_errors(Path(lane["expected_stdout"]["path"]))
        )
        if lane["allowed_dump_errors"] != expected_errors:
            raise ValueError(f"Reviewed private error contract changed: {name}")
        if name == "kraken_all_bins" and lane["required_private_coverage"] != {
            "successful_dumps": 31713,
            "observed_frames": 95139,
            "individual_bins": 1141668,
        }:
            raise ValueError("Full-corpus private-bin contract changed")
        if name == "kraken_all_bins" and lane.get("require_finite_private") is not True:
            raise ValueError("Full-corpus private numeric contract changed")
    claims = read(audit / "claims.json")["claims"]
    if evidence["classification_counts"] != dict(
        Counter(row["status"] for row in claims)
    ):
        raise ValueError("Independent classifications changed")
    if evidence["claims"] != [
        {key: row[key] for key in ("id", "status", "result", "limitations")}
        for row in claims
    ]:
        raise ValueError("Per-case independent findings changed")
    return manifest, evidence


def requests(path: Path):
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip():
                raise ValueError(f"Empty JSONL record: {path}:{number}")
            yield json.loads(line)


def require_fields(
    value: object, fields: list[str] | tuple[str, ...], context: str
) -> None:
    if not isinstance(value, dict) or not set(fields) <= value.keys():
        raise ValueError(f"Missing required observation fields: {context}")


def numeric_array(
    value: object, count: int, context: str, *, finite: bool = False
) -> None:
    if (
        not isinstance(value, list)
        or len(value) != count
        or any(
            x is not None and (isinstance(x, bool) or not isinstance(x, (float, int)))
            for x in value
        )
    ):
        raise ValueError(f"Invalid numeric observation array: {context}")
    if finite and any(x is None or not math.isfinite(x) for x in value):
        raise ValueError(f"Missing finite private numeric observations: {context}")


def float_observation(value: object, fields: list[str], context: str) -> None:
    if isinstance(value, dict) and "error" in value:
        if not isinstance(value["error"], str) or not value["error"]:
            raise ValueError(f"Invalid error observation: {context}")
        return
    require_fields(value, [*fields, "_bits", "_nonfinite"], context)
    if not isinstance(value["_bits"], dict) or set(value["_bits"]) != set(fields):
        raise ValueError(f"Missing exact bit observations: {context}")
    if not isinstance(value["_nonfinite"], dict) or not set(value["_nonfinite"]) <= set(
        fields
    ):
        raise ValueError(f"Invalid nonfinite observations: {context}")
    for field in fields:
        bits = value["_bits"][field]
        if not isinstance(bits, str) or len(bits) != 8:
            raise ValueError(f"Invalid f32 bits: {context}/{field}")
        try:
            decoded = struct.unpack("!f", bytes.fromhex(bits))[0]
        except (ValueError, struct.error) as exc:
            raise ValueError(f"Invalid f32 bits: {context}/{field}") from exc
        if math.isfinite(decoded):
            if (
                type(value[field]) not in (float, int)
                or value[field] != decoded
                or field in value["_nonfinite"]
            ):
                raise ValueError(f"Number/bit observation mismatch: {context}/{field}")
            if decoded == 0.0 and math.copysign(1, value[field]) != math.copysign(
                1, decoded
            ):
                raise ValueError(f"Signed zero/bit mismatch: {context}/{field}")
        elif value[field] is not None or value["_nonfinite"].get(field) != (
            "NaN" if math.isnan(decoded) else "inf" if decoded > 0 else "-inf"
        ):
            raise ValueError(f"Nonfinite/bit observation mismatch: {context}/{field}")


def validate_result(request: dict, result: dict) -> None:
    key = request["id"]
    operation = request.get("op", "geometry")
    if operation == "geometry":
        require_fields(
            result,
            [
                "atoms",
                "bonded_neighbors",
                "center",
                "sterimol_bond",
                "sterimol_dummy_raw",
                "sterimol_coordination",
                "pyramidalization",
            ],
            key,
        )
        if not isinstance(result["atoms"], list) or len(result["atoms"]) != len(
            request["atoms"]
        ):
            raise ValueError(f"Missing atom observations: {key}")
        for atom in result["atoms"]:
            require_fields(
                atom, ["element", "position", "vdw_radius", "covalent_radius"], key
            )
            numeric_array(atom["position"], 3, key)
        if not isinstance(result["bonded_neighbors"], list):
            raise ValueError(f"Missing donor neighbor observations: {key}")
        for name in ("sterimol_bond", "sterimol_dummy_raw", "sterimol_coordination"):
            float_observation(result[name], ["l", "b1", "b5"], f"{key}/{name}")
        float_observation(result["pyramidalization"], ["pyr_p", "pyr_alpha"], key)
        if not request.get("sterimol_only"):
            require_fields(result, ["buried_volume"], key)
            float_observation(result["buried_volume"], BV_FIELDS, key)
    elif operation == "kinetics":
        float_observation(result, KINETIC_FIELDS, key)
    elif operation == "aggregate":
        float_observation(result, AGGREGATE_FIELDS, key)
        if "error" not in result:
            require_fields(result, ["conformer_count"], key)
    elif operation == "model":
        fields = (
            ["t_quantile"]
            if "df" in request
            else ["calibration"]
            if "points" in request
            else [
                "design",
                "leverage",
                "mahalanobis",
                "prediction_interval",
                "confidence_interval",
                "t_multiplier",
                "assessment",
            ]
        )
        require_fields(result, fields, key)
    elif operation == "raw_occupancy":
        require_fields(
            result,
            [
                "buried_volume",
                "near_vbur",
                "far_vbur",
                "quadrants",
                "octants",
                "point_bins",
            ],
            key,
        )
        numeric_array(result["quadrants"], 4, key)
        numeric_array(result["octants"], 8, key)
        if len(result["point_bins"]) != len(request["points"]):
            raise ValueError(f"Missing raw point-bin observations: {key}")
    else:
        raise ValueError(f"Unsupported oracle operation: {operation}")


def validate_dump(
    request: dict, result: dict, allowed_errors: dict, *, finite: bool = False
) -> bool:
    key = request["id"]
    dump = result.get("dump")
    if not isinstance(dump, dict):
        raise ValueError(f"Missing requested private observations: {key}")
    if "error" in dump:
        if "error" not in result["buried_volume"] or allowed_errors.get(key) != dump:
            raise ValueError(
                "Unreviewed private dump error or successful public BV without bins: "
                f"{key}"
            )
        return False
    if key in allowed_errors:
        raise ValueError(f"Reviewed private error unexpectedly changed: {key}")
    require_fields(
        dump, ["center", "neighbors", "grid_count", "points", "orientations"], key
    )
    numeric_array(dump["center"], 3, key, finite=finite)
    neighbors = dump["neighbors"]
    if (
        not isinstance(neighbors, list)
        or len(neighbors) != 3
        or any(
            type(index) is not int or not 0 <= index < len(request["atoms"])
            for index in neighbors
        )
        or len(set(neighbors)) != 3
    ):
        raise ValueError(f"Missing private frame neighbors: {key}")
    if type(dump["grid_count"]) is not int or dump["grid_count"] < (1 if finite else 0):
        raise ValueError(f"Invalid private grid count: {key}")
    if request.get("include_points"):
        if (
            not isinstance(dump["points"], list)
            or len(dump["points"]) != dump["grid_count"]
        ):
            raise ValueError(f"Missing requested private grid points: {key}")
        for point in dump["points"]:
            numeric_array(point, 3, key, finite=finite)
    elif dump["points"] is not None:
        raise ValueError(f"Unexpected private point projection: {key}")
    orientations = dump["orientations"]
    if not isinstance(orientations, list) or len(orientations) != 3:
        raise ValueError(f"Expected all three frames: {key}")
    atom_count = sum(
        request.get("config", {}).get("include_hydrogens", False)
        or atom["element"].upper() != "H"
        for atom in request["atoms"]
    )
    for plane, frame in zip(neighbors, orientations, strict=True):
        require_fields(
            frame,
            [
                "plane",
                "basis",
                "aligned_atoms",
                "buried_volume",
                "near_vbur",
                "far_vbur",
                "quadrants",
                "octants",
            ],
            key,
        )
        if (
            type(frame["plane"]) is not int
            or frame["plane"] != plane
            or not isinstance(frame["basis"], list)
            or len(frame["basis"]) != 3
        ):
            raise ValueError(f"Incomplete or reordered private frames: {key}")
        for axis in frame["basis"]:
            numeric_array(axis, 3, key, finite=finite)
        numeric_array(frame["quadrants"], 4, key, finite=finite)
        numeric_array(frame["octants"], 8, key, finite=finite)
        numeric_array(
            [frame[name] for name in ("buried_volume", "near_vbur", "far_vbur")],
            3,
            key,
            finite=finite,
        )
        if (
            not isinstance(frame["aligned_atoms"], list)
            or len(frame["aligned_atoms"]) != atom_count
        ):
            raise ValueError(f"Missing aligned atom observations: {key}")
        for atom in frame["aligned_atoms"]:
            require_fields(atom, ["position", "radius_squared"], key)
            numeric_array(atom["position"], 3, key, finite=finite)
            numeric_array([atom["radius_squared"]], 1, key, finite=finite)
    return True


def coverage(
    request_path: Path,
    output_path: Path,
    allowed_dump_errors: dict | None = None,
    *,
    require_finite_private: bool = False,
) -> dict:
    """Require one ordered response per input; retain negative findings explicitly."""
    seen = set()
    identifiers = hashlib.sha256()
    errors = []
    nonfinite = []
    operation_states = Counter()
    frames = bins = 0
    successful_dumps = 0
    for number, (request, result) in enumerate(
        zip_longest(requests(request_path), requests(output_path)), 1
    ):
        if request is None or result is None:
            raise ValueError(f"Request/response count mismatch at row {number}")
        key = request.get("id")
        if not isinstance(key, str) or key in seen or result.get("id") != key:
            raise ValueError(
                f"Missing, duplicate, or out-of-order ID at row {number}: {key!r}"
            )
        seen.add(key)
        identifiers.update(canonical(key))
        validate_result(request, result)
        if "error" in result:
            errors.append({"id": key, "operation": "request", "error": result["error"]})
            operation_states["request/error"] += 1
        for operation, value in result.items():
            if isinstance(value, dict) and "error" in value:
                errors.append(
                    {"id": key, "operation": operation, "error": value["error"]}
                )
                operation_states[f"{operation}/error"] += 1
            elif isinstance(value, dict):
                operation_states[f"{operation}/observed"] += 1
            if isinstance(value, dict) and value.get("_nonfinite"):
                nonfinite.append(
                    {"id": key, "operation": operation, "values": value["_nonfinite"]}
                )
        if result.get("_nonfinite"):
            nonfinite.append(
                {"id": key, "operation": "request", "values": result["_nonfinite"]}
            )
        if request.get("dump") and not request.get("sterimol_only"):
            if validate_dump(
                request,
                result,
                allowed_dump_errors or {},
                finite=require_finite_private,
            ):
                successful_dumps += 1
                frames += 3
                bins += 36
    return {
        "rows": len(seen),
        "ordered_id_sha256": identifiers.hexdigest(),
        "operation_states": dict(sorted(operation_states.items())),
        "errors": errors,
        "nonfinite": nonfinite,
        "observed_frames": frames,
        "individual_bins": bins,
        "successful_dumps": successful_dumps,
    }


def public_digest(path: Path) -> str:
    """For the explicitly extended bin lane only, omit its added 'dump' field."""
    digest = hashlib.sha256()
    for row in requests(path):
        row.pop("dump", None)
        digest.update(canonical(row))
    return digest.hexdigest()


def evidence_paths(audit: Path, sealed: dict) -> list[str]:
    exact = {
        "CLAIMS.md",
        "claims.json",
        "SCIENTIFIC_ACCURACY_AUDIT.md",
        "REPRODUCE.md",
        "manifest_initial.json",
        "observer/manifest.json",
        "observer/Cargo.toml",
        "observer/Cargo.lock",
        "observer/src/main.rs",
        "observer/src/model_observation.rs",
        "geometry/observer_append.rs",
        "kraken/prepared_all/manifest.json",
        "kraken/prepared_all/inventory.jsonl",
        "kraken/primary/ligands_manifest.json",
    }
    roots = (
        "scripts/",
        "kraken/scripts/",
        "kinetics/results/",
        "models/results/",
        "kraken/analysis/",
        "kraken/interpretation/",
        "kraken/delta_interpretation/",
    )
    for name in sealed:
        path = Path(name)
        if name.startswith(roots) and path.suffix in {".py", ".json", ".jsonl", ".csv"}:
            exact.add(name)
        if name.startswith("geometry/") and (
            path.suffix == ".md"
            or path.name
            in {
                "reference_results.json",
                "reference_results.jsonl",
                "reference_selfcheck.json",
                "descriptor_comparisons.csv",
                "bin_comparisons.csv",
                "invariance.csv",
                "convergence.csv",
                "failures.csv",
                "metrics.json",
                "scientific_scale.json",
            }
        ):
            exact.add(name)
        if name.startswith("kraken/") and path.name in {
            "reference_raw.jsonl",
            "manifest_frozen.json",
        }:
            exact.add(name)
        if name.startswith(("kinetics/", "models/", "kraken/")) and path.name in {
            "REPORT.md",
            "KINETICS_CONFORMERS.md",
            "DISCREPANCIES.md",
        }:
            exact.add(name)
    missing = exact - sealed.keys()
    if missing:
        raise ValueError(
            f"Required evidence absent from sealed inventory: {sorted(missing)}"
        )
    return sorted(exact)


def frozen_dump_errors(path: Path) -> dict:
    return {
        row["id"]: row["dump"]
        for row in requests(path)
        if isinstance(row.get("dump"), dict) and "error" in row["dump"]
    }


def prepare(args) -> None:
    audit = args.audit.resolve(strict=True)
    root = new_directory(args.output, audit)
    (root / "helper").mkdir()
    helper = root / "helper/check_current_scientific_equivalence.py"
    shutil.copyfile(Path(__file__), helper)
    sealed_path = audit / "manifest_final.json"
    sealed = {item["path"]: item for item in read(sealed_path)["files"]}
    locked = {}

    def lock(relative):
        item = file_record(audit / relative)
        expected = sealed.get(relative)
        if expected is None or (item["sha256"], item["bytes"]) != (
            expected["sha256"],
            expected["bytes"],
        ):
            raise ValueError(
                f"Historical evidence differs from sealed inventory: {relative}"
            )
        locked[item["path"]] = item
        return item

    lanes = {}
    (root / "requests").mkdir()
    for name, (input_name, output_name, stderr_name, count) in LANES.items():
        source = lock(input_name)
        target = root / "requests" / f"{name}.jsonl"
        if input_name.endswith(".json"):
            with target.open("x") as stream:
                for row in json.loads(Path(source["path"]).read_text()):
                    stream.write(json.dumps(row) + "\n")
        else:
            shutil.copyfile(source["path"], target)
        expected = lock(output_name)
        allowed_errors = frozen_dump_errors(Path(expected["path"]))
        observed = coverage(target, Path(expected["path"]), allowed_errors)
        if observed["rows"] != count:
            raise ValueError(
                f"Unexpected frozen coverage: {name}: {observed['rows']} != {count}"
            )
        lanes[name] = {
            "requests": file_record(target),
            "source": source,
            "expected_stdout": expected,
            "expected_stderr": lock(stderr_name),
            "expected_returncode": 0,
            "expected_coverage": observed,
            "allowed_dump_errors": allowed_errors,
            "expected_rows": count,
            "comparison": "complete output bytes",
        }
        locked[str(target)] = file_record(target)
    target = root / "requests/kraken_all_bins.jsonl"
    with target.open("x") as stream:
        for row in requests(Path(lanes["kraken"]["requests"]["path"])):
            row.update(dump=True, include_points=False)
            stream.write(json.dumps(row) + "\n")
    lanes["kraken_all_bins"] = {
        "requests": file_record(target),
        "source": lanes["kraken"]["source"],
        "expected_rows": 31721,
        "public_expected_stdout": lanes["kraken"]["expected_stdout"],
        "public_expected_stderr": lanes["kraken"]["expected_stderr"],
        "expected_returncode": 0,
        "allowed_dump_errors": KRAKEN_DUMP_ERRORS,
        "require_finite_private": True,
        "required_private_coverage": {
            "successful_dumps": 31713,
            "observed_frames": 95139,
            "individual_bins": 1141668,
        },
        "comparison": (
            "new baseline/candidate complete bytes; frozen public projection fidelity"
        ),
    }
    public_errors = {
        row["id"]: row["buried_volume"]
        for row in requests(Path(lanes["kraken"]["expected_stdout"]["path"]))
        if "error" in row["buried_volume"]
    }
    if any(
        public_errors.get(key) != error for key, error in KRAKEN_DUMP_ERRORS.items()
    ):
        raise ValueError("Reviewed dump errors do not match sealed public evidence")
    locked[str(target)] = file_record(target)
    evidence_files = [lock(name) for name in evidence_paths(audit, sealed)]
    claims = read(audit / "claims.json")["claims"]
    evidence = {
        "purpose": PURPOSE,
        "sealed_manifest": file_record(sealed_path),
        "classification_counts": dict(Counter(row["status"] for row in claims)),
        "claims": [
            {key: row[key] for key in ("id", "status", "result", "limitations")}
            for row in claims
        ],
        "files": evidence_files,
        "limits_policy": (
            "Per-case evidence and exact outputs locked; "
            "no global tolerance inferred from observed maxima."
        ),
        "reference_policy": (
            "Immutable independent results retained; "
            "this helper does not rerun or modify reference algorithms."
        ),
    }
    write_new(root / "independent_evidence_lock.json", evidence)
    locked[str(sealed_path)] = file_record(sealed_path)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "kind": "live_scientific_oracle_inputs",
        "purpose": PURPOSE,
        "created_utc": datetime.now(UTC).isoformat(),
        "audit": str(audit),
        "runtime": runtime(),
        "helper": file_record(helper),
        "lanes": lanes,
        "locked_files": list(locked.values()),
        "evidence_lock": file_record(root / "independent_evidence_lock.json"),
        "excluded_scope": [
            "native CLI commands",
            "Python preparation",
            "new independent reference calculations",
            "experimental validity",
        ],
    }
    write_new(root / "manifest.json", manifest)
    print(
        json.dumps(
            {
                "prepared": str(root),
                "lanes": list(lanes),
                "classification_counts": evidence["classification_counts"],
                "purpose": PURPOSE,
            },
            indent=2,
        )
    )


def source_records(source: Path) -> list[dict]:
    files = [
        source / "Cargo.toml",
        source / "Cargo.lock",
        *sorted((source / "src").rglob("*.rs")),
    ]
    if (source / "build.rs").exists():
        files.append(source / "build.rs")
    return [file_record(path) for path in files]


def check_dependency_locks(source: Path, adapter: Path) -> None:
    """Do not silently compile live source against older frozen dependencies."""

    def dependencies(path):
        packages = tomllib.loads(path.read_text())["package"]
        return {
            (row["name"], row["version"], row.get("source"), row.get("checksum"))
            for row in packages
            if row["name"] not in {"steric_x", "stericx-audit-observer"}
        }

    missing = dependencies(source) - dependencies(adapter)
    if missing:
        raise ValueError(
            "Observer lock differs from live library dependencies; "
            f"an explicit adapter-lock review is required: {sorted(missing)}"
        )


def build(args) -> None:
    oracle = args.oracle.resolve(strict=True)
    prepared, _ = load_oracle(oracle)
    audit = Path(prepared["audit"])
    original_source = args.source.resolve(strict=True)
    output = new_directory(args.output, audit)
    original_sources = source_records(original_source)
    # A later candidate may edit the working tree. Keep the actual compiled
    # baseline source immutable and verifiable independently of that worktree.
    source = output / "source"
    for item in original_sources:
        relative = Path(item["path"]).relative_to(original_source)
        destination = source / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(verify(item), destination)
    captured_sources = source_records(source)
    crate = output / "observer"
    (crate / "src").mkdir(parents=True)
    shutil.copyfile(audit / "observer/Cargo.lock", crate / "Cargo.lock")
    check_dependency_locks(source / "Cargo.lock", crate / "Cargo.lock")
    cargo = (audit / "observer/Cargo.toml").read_text()
    old = 'steric_x = { path = "../frozen/repository" }'
    if cargo.count(old) != 1:
        raise ValueError("Unrecognized frozen observer Cargo dependency")
    cargo = cargo.replace(old, "steric_x = { path = " + json.dumps(str(source)) + " }")
    (crate / "Cargo.toml").write_text(cargo)
    main = (audit / "observer/src/main.rs").read_text()
    old = '#[path = "../../geometry/frozen_buried_volume.rs"]'
    if main.count(old) != 1:
        raise ValueError("Unrecognized frozen private-observation module")
    (crate / "src/main.rs").write_text(
        main.replace(old, '#[path = "observed_buried_volume.rs"]')
    )
    shutil.copyfile(
        audit / "observer/src/model_observation.rs", crate / "src/model_observation.rs"
    )
    prefix = (source / "src/geometry/buried_volume.rs").read_bytes()
    append = (audit / "geometry/observer_append.rs").read_bytes()
    observed = crate / "src/observed_buried_volume.rs"
    observed.write_bytes(prefix + append)
    if (
        observed.read_bytes()[: len(prefix)] != prefix
        or observed.read_bytes()[len(prefix) :] != append
    ):
        raise ValueError("Live source prefix or unchanged visibility adapter mismatch")
    versions = {
        "rustc": tool_version(["rustc", "-vV"]),
        "cargo": tool_version(["cargo", "-V"]),
    }
    if any(item["returncode"] != 0 for item in versions.values()):
        raise ValueError("Rust toolchain is unavailable")
    target = output / "target"
    command = [
        "cargo",
        "build",
        "--locked",
        "--release",
        "--manifest-path",
        str(crate / "Cargo.toml"),
        "--target-dir",
        str(target),
    ]
    if not args.online:
        command.append("--offline")
    started = {
        "schema_version": SCHEMA_VERSION,
        "kind": "live_scientific_observer_build",
        "purpose": PURPOSE,
        "source": str(source),
        "original_source": str(original_source),
        "original_sources": original_sources,
        "oracle_manifest": file_record(oracle / "manifest.json"),
        "helper": prepared["helper"],
        "sources": captured_sources,
        "adapter_files": [
            file_record(path) for path in sorted(crate.rglob("*")) if path.is_file()
        ],
        "bv_source_sha256": hashlib.sha256(prefix).hexdigest(),
        "bv_prefix_bytes": len(prefix),
        "unchanged_append_sha256": hashlib.sha256(append).hexdigest(),
        "original_prefix_verified": True,
        "toolchain": versions,
        "runtime": runtime(),
        "command": command,
        "compiler_environment": {
            name: os.environ.get(name)
            for name in (
                "RUSTFLAGS",
                "CARGO_ENCODED_RUSTFLAGS",
                "RUSTC",
                "RUSTC_WRAPPER",
                "CARGO_BUILD_TARGET",
            )
        },
        "git": tool_version(["git", "rev-parse", "HEAD"], original_source),
    }
    write_new(output / "build_started.json", started)
    with (
        (output / "build.stdout").open("xb") as stdout,
        (output / "build.stderr").open("xb") as stderr,
    ):
        result = subprocess.run(
            command, cwd=source, stdout=stdout, stderr=stderr, check=False
        )
    for item in original_sources + captured_sources + started["adapter_files"]:
        verify(item)
    finished = {
        **started,
        "returncode": result.returncode,
        "stdout": file_record(output / "build.stdout"),
        "stderr": file_record(output / "build.stderr"),
    }
    if result.returncode == 0:
        finished["binary"] = file_record(target / "release/stericx-audit-observer")
    write_new(output / "manifest.json", finished)
    if result.returncode:
        raise ValueError(f"Build failed; logs retained in {output}")
    print(
        json.dumps(
            {"build": str(output), "binary": finished["binary"], "purpose": PURPOSE},
            indent=2,
        )
    )


def load_build(build_path: Path, oracle: Path, prepared: dict) -> dict:
    compiled = read(build_path / "manifest.json")
    if (
        compiled.get("kind") != "live_scientific_observer_build"
        or compiled.get("schema_version") != SCHEMA_VERSION
        or compiled.get("returncode") != 0
    ):
        raise ValueError("A successful live-source observer build is required")
    if compiled["oracle_manifest"] != file_record(oracle / "manifest.json"):
        raise ValueError(
            "Build and observations must use the identical prepared oracle"
        )
    if compiled["helper"] != prepared["helper"]:
        raise ValueError("Build and oracle helper identities differ")
    verify(compiled["helper"])
    if compiled["binary"] != file_record(
        build_path / "target/release/stericx-audit-observer"
    ):
        raise ValueError("Build binary is not the recorded build output")
    for item in compiled["sources"] + compiled["adapter_files"]:
        verify(item)
    source = build_path / "source"
    if Path(compiled["source"]) != source or compiled["sources"] != source_records(
        source
    ):
        raise ValueError("Compiled source snapshot inventory changed")
    original = {
        str(Path(item["path"]).relative_to(compiled["original_source"])): (
            item["sha256"],
            item["bytes"],
        )
        for item in compiled["original_sources"]
    }
    captured = {
        str(Path(item["path"]).relative_to(source)): (item["sha256"], item["bytes"])
        for item in compiled["sources"]
    }
    if original != captured:
        raise ValueError("Build snapshot differs from recorded original source")
    crate = build_path / "observer"
    if compiled["adapter_files"] != [
        file_record(path) for path in sorted(crate.rglob("*")) if path.is_file()
    ]:
        raise ValueError("Observation adapter inventory changed")
    audit = Path(prepared["audit"])
    expected_cargo = (
        (audit / "observer/Cargo.toml")
        .read_text()
        .replace(
            'steric_x = { path = "../frozen/repository" }',
            "steric_x = { path = " + json.dumps(str(source)) + " }",
        )
    )
    expected_main = (
        (audit / "observer/src/main.rs")
        .read_text()
        .replace(
            '#[path = "../../geometry/frozen_buried_volume.rs"]',
            '#[path = "observed_buried_volume.rs"]',
        )
    )
    if (crate / "Cargo.toml").read_text() != expected_cargo or (
        crate / "src/main.rs"
    ).read_text() != expected_main:
        raise ValueError(
            "Observer no longer links the captured source and visibility adapter"
        )
    for relative in ("Cargo.lock", "src/model_observation.rs"):
        if (crate / relative).read_bytes() != (
            audit / "observer" / relative
        ).read_bytes():
            raise ValueError(f"Frozen observer component changed: {relative}")
    prefix = (source / "src/geometry/buried_volume.rs").read_bytes()
    append = (audit / "geometry/observer_append.rs").read_bytes()
    if (
        (crate / "src/observed_buried_volume.rs").read_bytes() != prefix + append
        or compiled["bv_source_sha256"] != hashlib.sha256(prefix).hexdigest()
        or compiled["bv_prefix_bytes"] != len(prefix)
        or compiled["unchanged_append_sha256"] != hashlib.sha256(append).hexdigest()
        or compiled["original_prefix_verified"] is not True
    ):
        raise ValueError("Current BV prefix or unchanged append identity mismatch")
    for key in ("stdout", "stderr"):
        if compiled[key] != file_record(build_path / f"build.{key}"):
            raise ValueError("Build log identity mismatch")
    check_dependency_locks(source / "Cargo.lock", crate / "Cargo.lock")
    return compiled


def lane_coverage_and_fidelity(
    lane: dict, stdout: Path, stderr: Path
) -> tuple[dict, bool]:
    observed = coverage(
        verify(lane["requests"]),
        stdout,
        lane["allowed_dump_errors"],
        require_finite_private=lane.get("require_finite_private", False),
    )
    if observed["rows"] != lane["expected_rows"]:
        raise ValueError("Incomplete observation lane")
    for key, expected in lane.get("required_private_coverage", {}).items():
        if observed[key] != expected:
            raise ValueError(f"Incomplete required private coverage: {key}")
    if "expected_stdout" in lane:
        expected_coverage = coverage(
            verify(lane["requests"]),
            verify(lane["expected_stdout"]),
            lane["allowed_dump_errors"],
        )
        if lane["expected_coverage"] != expected_coverage:
            raise ValueError(
                "Prepared expected coverage differs from sealed output bytes"
            )
        faithful = (
            sha(stdout) == lane["expected_stdout"]["sha256"]
            and sha(stderr) == lane["expected_stderr"]["sha256"]
            and observed == expected_coverage
        )
    else:
        faithful = (
            public_digest(stdout)
            == public_digest(verify(lane["public_expected_stdout"]))
            and sha(stderr) == lane["public_expected_stderr"]["sha256"]
        )
    return observed, faithful


def observe(args) -> None:
    oracle = args.oracle.resolve(strict=True)
    prepared, evidence = load_oracle(oracle)
    build_path = args.build.resolve(strict=True)
    compiled = load_build(build_path, oracle, prepared)
    binary = verify(compiled["binary"])
    output = new_directory(args.output, Path(prepared["audit"]))
    environment = dict(os.environ)
    environment.update(
        LC_ALL="C",
        TZ="UTC",
        RAYON_NUM_THREADS=str(args.threads),
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
    )
    environment.pop("STERICX_PROFILE_PATH", None)
    names = args.lanes or list(ALL_LANES)
    if len(set(names)) != len(names):
        raise ValueError("Duplicate observation lanes")
    started = {
        "schema_version": SCHEMA_VERSION,
        "kind": "live_scientific_observations",
        "purpose": PURPOSE,
        "oracle_manifest": file_record(oracle / "manifest.json"),
        "evidence_lock": prepared["evidence_lock"],
        "build_manifest": file_record(build_path / "manifest.json"),
        "binary": compiled["binary"],
        "helper": prepared["helper"],
        "environment": {
            key: environment[key]
            for key in (
                "LC_ALL",
                "TZ",
                "RAYON_NUM_THREADS",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
            )
        },
        "classification_counts": evidence["classification_counts"],
        "requested_lanes": names,
        "full_observer_scope": set(names) == set(ALL_LANES),
        "excluded_scope": prepared["excluded_scope"],
    }
    write_new(output / "observation_started.json", started)
    results = {}
    failure = None
    try:
        for name in names:
            lane = prepared["lanes"][name]
            directory = output / name
            directory.mkdir()
            request_path = verify(lane["requests"])
            before = time.monotonic()
            timed_out = False
            with (
                request_path.open("rb") as stdin,
                (directory / "stdout.jsonl").open("xb") as stdout,
                (directory / "stderr.txt").open("xb") as stderr,
            ):
                try:
                    process = subprocess.run(
                        [str(binary)],
                        stdin=stdin,
                        stdout=stdout,
                        stderr=stderr,
                        env=environment,
                        timeout=args.timeout,
                        check=False,
                    )
                    returncode = process.returncode
                except subprocess.TimeoutExpired:
                    timed_out = True
                    returncode = None
            row = {
                "requests": lane["requests"],
                "stdout": file_record(directory / "stdout.jsonl"),
                "stderr": file_record(directory / "stderr.txt"),
                "returncode": returncode,
                "timed_out": timed_out,
                "seconds": time.monotonic() - before,
            }
            results[name] = row
            if timed_out or returncode != lane["expected_returncode"]:
                raise ValueError(
                    f"Observer failed in {name}: {returncode}, timeout={timed_out}"
                )
            row["coverage"], row["frozen_fidelity"] = lane_coverage_and_fidelity(
                lane, directory / "stdout.jsonl", directory / "stderr.txt"
            )
            if "expected_stdout" in lane:
                row["fidelity_scope"] = "complete bytes, process status and coverage"
            else:
                row["fidelity_scope"] = (
                    "original public fields; new private bins require "
                    "a baseline/candidate comparison"
                )
            write_new(directory / "result.json", row)
            print(
                f"{name}: {row['coverage']['rows']} rows; "
                f"frozen fidelity={row['frozen_fidelity']}",
                flush=True,
            )
            if args.require_frozen_fidelity and not row["frozen_fidelity"]:
                raise ValueError(f"Frozen observation fidelity mismatch: {name}")
        load_build(build_path, oracle, prepared)
    except Exception as exc:
        failure = str(exc)
    finished = {
        **started,
        "complete": failure is None,
        "error": failure,
        "lanes": results,
        "frozen_fidelity": failure is None
        and all(row["frozen_fidelity"] for row in results.values()),
    }
    write_new(output / "manifest.json", finished)
    if failure:
        raise ValueError(f"{failure}; partial observations retained at {output}")
    print(
        json.dumps(
            {
                "observations": str(output),
                "full_observer_scope": started["full_observer_scope"],
                "frozen_fidelity": finished["frozen_fidelity"],
                "purpose": PURPOSE,
            },
            indent=2,
        )
    )


def validate_observations(root: Path) -> tuple[dict, dict]:
    """Recompute admission from bytes; manifest flags are assertions, not evidence."""
    manifest = read(root / "manifest.json")
    if (
        manifest.get("kind") != "live_scientific_observations"
        or manifest.get("schema_version") != SCHEMA_VERSION
        or manifest.get("complete") is not True
        or manifest.get("error") is not None
    ):
        raise ValueError(
            "Complete current-schema live-source observations are required"
        )
    oracle = verify(manifest["oracle_manifest"]).parent
    prepared, evidence = load_oracle(oracle)
    build_path = verify(manifest["build_manifest"]).parent
    compiled = load_build(build_path, oracle, prepared)
    for key, expected in (
        ("binary", compiled["binary"]),
        ("helper", prepared["helper"]),
        ("evidence_lock", prepared["evidence_lock"]),
        ("classification_counts", evidence["classification_counts"]),
        ("excluded_scope", prepared["excluded_scope"]),
    ):
        if manifest[key] != expected:
            raise ValueError(
                f"Observation identity is not linked to its oracle/build: {key}"
            )
    names = manifest["requested_lanes"]
    if (
        not isinstance(names, list)
        or not names
        or len(names) != len(set(names))
        or not set(names) <= set(ALL_LANES)
        or set(names) != set(manifest["lanes"])
    ):
        raise ValueError(
            "Requested and observed lane sets must be nonempty and identical"
        )
    full_scope = set(names) == set(ALL_LANES)
    if manifest["full_observer_scope"] is not full_scope:
        raise ValueError("Forged full-observer coverage flag")
    for name in names:
        row = manifest["lanes"][name]
        lane = prepared["lanes"][name]
        if row["requests"] != lane["requests"]:
            raise ValueError(
                f"Observation request identity is not tied to prepared lane: {name}"
            )
        for key, filename in (("stdout", "stdout.jsonl"), ("stderr", "stderr.txt")):
            if row[key] != file_record(root / name / filename):
                raise ValueError(f"Observation stream identity changed: {name}/{key}")
        if (
            row["returncode"] != lane["expected_returncode"]
            or row["timed_out"] is not False
        ):
            raise ValueError(f"Unsuccessful observer process: {name}")
        observed, faithful = lane_coverage_and_fidelity(
            lane, verify(row["stdout"]), verify(row["stderr"])
        )
        if row["coverage"] != observed or row["frozen_fidelity"] is not faithful:
            raise ValueError(
                f"Recorded coverage/fidelity differs from current stream bytes: {name}"
            )
        if read(root / name / "result.json") != row:
            raise ValueError(f"Lane receipt and observation manifest differ: {name}")
    faithful = all(row["frozen_fidelity"] for row in manifest["lanes"].values())
    if manifest["frozen_fidelity"] is not faithful:
        raise ValueError("Forged aggregate frozen-fidelity flag")
    return manifest, prepared


def comparison(args) -> None:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    base_path = args.baseline.resolve(strict=True) / "manifest.json"
    candidate_path = args.candidate.resolve(strict=True) / "manifest.json"
    validated = {}
    for path in (base_path, candidate_path):
        if path not in validated:
            validated[path] = validate_observations(path.parent)
        _, prepared = validated[path]
        if output.is_relative_to(Path(prepared["audit"]).resolve()):
            raise ValueError("Comparison output must be outside the historical audit")
    baseline, candidate = validated[base_path][0], validated[candidate_path][0]
    if (
        baseline["oracle_manifest"] != candidate["oracle_manifest"]
        or baseline["evidence_lock"] != candidate["evidence_lock"]
    ):
        raise ValueError(
            "Both observations must share the exact frozen inputs "
            "and independent evidence lock"
        )
    if set(baseline["lanes"]) != set(candidate["lanes"]):
        raise ValueError(
            "Observation lane sets differ; missing coverage cannot be ignored"
        )
    rows = {}
    for name, expected in baseline["lanes"].items():
        actual = candidate["lanes"][name]
        checks = {
            "stdout_exact": expected["stdout"]["sha256"] == actual["stdout"]["sha256"],
            "stderr_exact": expected["stderr"]["sha256"] == actual["stderr"]["sha256"],
            "returncode_exact": expected["returncode"] == actual["returncode"],
            "requests_exact": expected["requests"] == actual["requests"],
            "coverage_exact": expected["coverage"] == actual["coverage"],
        }
        rows[name] = {
            **checks,
            "exact": all(checks.values()),
            "rows": expected["coverage"]["rows"],
        }
    equivalent = all(row["exact"] for row in rows.values())
    anchored = baseline["frozen_fidelity"] and candidate["frozen_fidelity"]
    result = {
        "schema_version": SCHEMA_VERSION,
        "purpose": PURPOSE,
        "baseline": file_record(base_path),
        "candidate": file_record(candidate_path),
        "helper": baseline["helper"],
        "equivalent": equivalent,
        "admitted": equivalent and anchored,
        "lanes": rows,
        "full_observer_scope": baseline["full_observer_scope"]
        and candidate["full_observer_scope"],
        "frozen_observations_preserved": anchored,
        "independent_evidence_lock": baseline["evidence_lock"],
        "covered_per_case_residuals_preserved": equivalent and anchored,
        "all_observer_residual_evidence_preserved": (
            equivalent
            and anchored
            and baseline["full_observer_scope"]
            and candidate["full_observer_scope"]
        ),
        "reference_kernels_reexecuted": False,
        "classification_counts": baseline["classification_counts"],
        "scientific_pass_claimed": False,
        "excluded_scope": baseline["excluded_scope"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    write_new(output, result)
    print(json.dumps(result, indent=2))
    if not equivalent or not anchored:
        raise ValueError(
            "Exact anchored equivalence required; "
            "comparison and raw observations retained"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    actions = parser.add_subparsers(dest="action", required=True)
    command = actions.add_parser(
        "prepare",
        help="Lock sealed fixtures, findings and independent residual evidence",
    )
    command.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    command.add_argument("--output", type=Path, required=True)
    command.set_defaults(function=prepare)
    command = actions.add_parser(
        "build",
        help="Build an adapter from the selected LIVE source; no frozen kernels",
    )
    command.add_argument("--oracle", type=Path, required=True)
    command.add_argument("--source", type=Path, default=REPO)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument(
        "--online",
        action="store_true",
        help="Permit Cargo dependency downloads (default is offline)",
    )
    command.set_defaults(function=build)
    command = actions.add_parser(
        "observe",
        help="Run selected frozen requests through a built live-source adapter",
    )
    command.add_argument("--oracle", type=Path, required=True)
    command.add_argument("--build", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--lanes", nargs="+", choices=ALL_LANES)
    command.add_argument("--threads", type=int, default=1)
    command.add_argument(
        "--timeout",
        type=float,
        default=600,
        help="Per-lane process time limit in seconds",
    )
    command.add_argument("--require-frozen-fidelity", action="store_true")
    command.set_defaults(function=observe)
    command = actions.add_parser(
        "compare", help="Compare two observations exactly; retain scientific failures"
    )
    command.add_argument("--baseline", type=Path, required=True)
    command.add_argument("--candidate", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.set_defaults(function=comparison)
    args = parser.parse_args()
    if getattr(args, "threads", 1) < 1 or getattr(args, "timeout", 1) <= 0:
        parser.error("threads and timeout must be positive")
    try:
        args.function(args)
    except (OSError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
