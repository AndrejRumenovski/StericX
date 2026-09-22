#!/usr/bin/env python3
"""Replay retained thermodynamics witnesses without changing the original audit.

Run from the repository with .venv/bin/python; native executables must be supplied.
Outputs are deliberately non-overwriting. This is a numerical/domain replay,
not a performance measurement or a global scientific validation claim.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import shutil
import struct
import subprocess
import sys
from decimal import Decimal as D
from decimal import localcontext
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
AUDIT = REPO / "docs/scientific_accuracy_audit/kinetics"
KB, NA, H = D("1.380649e-23"), D("6.02214076e23"), D("6.62607015e-34")
R = KB * NA / D(4184)
HARTREE_KCAL = D("627.5094740628974597648183556")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, (float, np.floating)) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: clean(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(item) for item in value]
    return value


def write(path, value):
    path.write_text(json.dumps(clean(value), indent=2, sort_keys=True) + "\n")


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def probabilities(energies, temperature):
    with localcontext() as context:
        context.prec = 80
        values = [D(float(value)) for value in energies]
        low = min(values)
        weights = [
            (-(value - low) / (R * D(float(temperature)))).exp() for value in values
        ]
        return [float(value / sum(weights)) for value in weights]


def capture_python(output):
    prep = module("remediated_prepare_data", output / "source/scripts/prepare_data.py")
    quantum = module("remediated_quantum", output / "source/scripts/stericx_quantum.py")
    records = {"crest": [], "mmff": [], "ee": []}
    for temperature in (298.15, 500):
        backend = object.__new__(quantum.QuantumBackend)
        backend.config = quantum.QuantumConfig(
            cache_dir=output / "unused", temperature_k=temperature
        )
        frames = [
            quantum.XyzFrame(
                ("C",), np.zeros((1, 3)), "analytic", -100 + value / float(HARTREE_KCAL)
            )
            for value in (0, 1)
        ]
        for name in ("default", "missing", "negative"):
            row = {"summary": name, "temperature": temperature}
            try:
                result = backend._conformer_thermodynamics(
                    frames, AUDIT / f"inputs/crest_{name}.log"
                )
                row["result"] = result
                row["reference_weights"] = probabilities([0, 1], temperature)
                row["max_abs_error"] = max(
                    abs(actual["boltzmann_weight"] - expected)
                    for actual, expected in zip(
                        result, row["reference_weights"], strict=True
                    )
                )
                assert name != "negative" and row["max_abs_error"] < 2e-12
            except quantum.QuantumBackendError as exc:
                row["error"] = str(exc)
                assert name == "negative"
            records["crest"].append(row)
    for case in json.loads((AUDIT / "inputs/ensembles.json").read_text()):
        energies = [float(value) for value in case["energies"]]
        statuses = case.get("statuses", [0] * len(energies))
        row = {"id": case["id"]}
        try:
            with (
                patch.object(
                    prep.AllChem,
                    "EmbedMultipleConfs",
                    return_value=tuple(range(len(energies))),
                ),
                patch.object(
                    prep.AllChem, "MMFFHasAllMoleculeParams", return_value=True
                ),
                patch.object(
                    prep.AllChem,
                    "MMFFOptimizeMoleculeConfs",
                    return_value=list(zip(statuses, energies, strict=True)),
                ),
            ):
                result = prep.embed_and_optimize(
                    "CP(C)C",
                    20260919,
                    len(energies),
                    0.1,
                    case.get("energy_window", 2000),
                    case["temperature"],
                )
        except ValueError as exc:
            assert case["id"] == "failed_all", case
            row["error"] = str(exc)
            records["mmff"].append(row)
            continue
        reference = probabilities(
            [energies[index] for index in result.conformer_ids], case["temperature"]
        )
        row.update(
            ids=result.conformer_ids,
            weights=result.boltzmann_weights,
            statuses=result.mmff_statuses,
            reference_weights=reference,
            max_abs_error=max(
                abs(actual - expected)
                for actual, expected in zip(
                    result.boltzmann_weights, reference, strict=True
                )
            ),
        )
        assert row["max_abs_error"] < 2e-14
        records["mmff"].append(row)
    cases = json.loads((AUDIT / "inputs/ee.json").read_text())
    for temperature in cases["temperatures"]:
        for entry in cases["values"]:
            value = float(entry)
            row = {"ee": value, "temperature": temperature}
            try:
                actual = float(prep.ee_to_ddg(pd.Series([value]), temperature).iloc[0])
                row["ddg"] = actual
                if math.isnan(value):
                    assert math.isnan(actual)
                else:
                    with localcontext() as context:
                        context.prec = 80
                        fraction = abs(D(value)) / 100
                        reference = float(
                            R * D(temperature) * ((1 + fraction) / (1 - fraction)).ln()
                        )
                    row.update(reference=reference, abs_error=abs(actual - reference))
                    assert abs(actual - reference) < 2e-11
            except ValueError as exc:
                row["error"] = str(exc)
                assert abs(value) >= 100 or math.isinf(value)
            records["ee"].append(row)
    source = pd.read_csv(REPO / "data/official/ni_hda_kraken.csv")
    normalized = prep.normalize_public_sigman(source)
    records["ni_hda"] = {
        "row_count": len(normalized),
        "temperatures": sorted(normalized.Temp_K.unique().tolist()),
        "published_targets_preserved": bool(
            np.allclose(
                normalized.Experimental_ddG,
                pd.to_numeric(source.ddG_abs, errors="coerce"),
                equal_nan=True,
            )
        ),
    }
    assert records["ni_hda"]["temperatures"] == [353.15]
    assert records["ni_hda"]["published_targets_preserved"]
    write(output / "python.json", records)
    return records


def capture_native(output, binary, observer):
    raw = output / "native"
    raw.mkdir()
    runs = []
    for path in sorted((AUDIT / "inputs").glob("parse_*.csv")):
        dest = raw / (path.stem + ".sigpack")
        command = [
            str(binary),
            "parse",
            "--csv",
            str(path),
            "--xyz-dir",
            str(AUDIT / "inputs"),
            "--output",
            str(dest),
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        row = {
            "id": path.stem,
            "command": command,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
        if result.returncode == 0:
            values = struct.unpack("=16f", dest.read_bytes())
            row.update(packed_sha256=digest(dest), packed_f32=values)
        runs.append(row)
    for difference in (0, 1, 20, -1):
        command = [
            str(binary),
            "simulate",
            "--ddg",
            str(difference),
            "--temp",
            "298.15",
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        assert result.returncode == 0 and "rate_constant" not in result.stdout
        runs.append(
            {
                "id": f"simulate_{difference}",
                "command": command,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
            }
        )
    write(raw / "cli_runs.json", runs)
    expected_success = {
        "parse_provided",
        "parse_weights_missing",
        "parse_tiny_weights",
        "parse_relative_offset",
        "parse_single",
    }
    for row in runs:
        if not row["id"].startswith("parse_"):
            continue
        assert (row["returncode"] == 0) == (row["id"] in expected_success), row
        if row["returncode"] == 0:
            weight = (
                0.5
                if row["id"] == "parse_weights_missing"
                else 0
                if row["id"] == "parse_single"
                else 0.75
            )
            expected = [3.7 + 2 * weight, 1.7, 3.47 + 2 * weight]
            assert (
                max(
                    abs(a - b)
                    for a, b in zip(row["packed_f32"][:3], expected, strict=True)
                )
                < 2e-6
            )
            assert row["packed_f32"][14] == (0 if row["id"] == "parse_single" else 1)
    for group in ("kinetics", "aggregation"):
        cases = json.loads((AUDIT / f"inputs/{group}.json").read_text())
        result = subprocess.run(
            [str(observer)],
            input="".join(json.dumps(case) + "\n" for case in cases),
            text=True,
            capture_output=True,
            check=True,
        )
        (raw / f"{group}.jsonl").write_text(result.stdout)
        (raw / f"{group}.stderr").write_text(result.stderr)
        rows = [json.loads(line) for line in result.stdout.splitlines()]
        assert [row["id"] for row in rows] == [case["id"] for case in cases]
    return runs


def analyze(output):
    rows = []
    for row in map(
        json.loads, (output / "native/kinetics.jsonl").read_text().splitlines()
    ):
        energy, temperature = row["input_ddg"], row["input_temperature"]
        if energy is None or temperature is None or temperature <= 0:
            assert all(
                row[field] is None
                for field in ("rate", "ee_percent", "r_percent", "s_percent")
            )
            continue
        if abs(energy / (float(R) * temperature)) > 1e6:
            continue
        with localcontext() as context:
            context.prec = 80
            log_ratio = D(energy) / (R * D(temperature))
            rate = KB * D(temperature) / H * (-log_ratio).exp()
            terms = [D(1), (-abs(log_ratio)).exp()]
            major, minor = [100 * term / sum(terms) for term in terms]
            r, s = (major, minor) if energy >= 0 else (minor, major)
            reference = {
                "rate": rate,
                "ee_percent": major - minor,
                "r_percent": r,
                "s_percent": s,
            }
        for field, expected in reference.items():
            actual = row[field]
            record = {
                "id": row["id"],
                "field": field,
                "actual": actual,
                "reference_decimal": str(expected),
            }
            if actual is not None:
                error = abs(actual - float(expected))
                record.update(
                    abs_error=error,
                    relative_error=(
                        error / float(expected) if float(expected) != 0 else 0
                    ),
                )
                # Absolute f32 rounding bound, including subnormal rates.
                assert error <= max(abs(float(expected)) * 1.3e-7, 7.1e-46), record
            else:
                assert field == "rate" and expected > D(
                    float(np.finfo(np.float32).max)
                ), record
            rows.append(record)
    write(output / "kinetics_comparison.json", rows)
    agg_inputs = {
        case["id"]: case
        for case in json.loads((AUDIT / "inputs/aggregation.json").read_text())
    }
    mapping = {
        "vbur_boltz": "buried_volume",
        "qvbur_min_boltz": "qvbur_min",
        "qvbur_max_boltz": "qvbur_max",
        "max_delta_qvbur_boltz": "max_delta_qvbur",
        "near_vbur_boltz": "near_vbur",
        "far_vbur_boltz": "far_vbur",
    }
    aggregated = []
    for row in map(
        json.loads, (output / "native/aggregation.jsonl").read_text().splitlines()
    ):
        case = agg_inputs[row["id"]]
        weights = [D(float(np.float32(value))) for value in case["weights"]]
        invalid = (
            not weights
            or any(not value.is_finite() or value < 0 for value in weights)
            or sum(weights) <= 0
            or any(
                not math.isfinite(float(value))
                for conformer in case["conformers"]
                for value in conformer.values()
            )
        )
        assert ("error" in row) == invalid, row
        if invalid:
            aggregated.append({"id": row["id"], "error": row["error"]})
            continue
        for field, source in mapping.items():
            reference = sum(
                D(float(np.float32(conformer.get(source, 0)))) * weight
                for conformer, weight in zip(case["conformers"], weights, strict=True)
            ) / sum(weights)
            error = abs(row[field] - float(reference))
            assert error <= max(float(abs(reference)) * 6e-8, 7.1e-46)
            aggregated.append(
                {
                    "id": row["id"],
                    "field": field,
                    "actual": row[field],
                    "reference_decimal": str(reference),
                    "abs_error": error,
                }
            )
    write(output / "aggregation_comparison.json", aggregated)
    return {
        "kinetic_scalar_comparisons": len(rows),
        "aggregation_scalar_or_rejection_comparisons": len(aggregated),
        "rate_max_relative_error_normal": max(
            row["relative_error"]
            for row in rows
            if row["field"] == "rate"
            and row["actual"] is not None
            and row["actual"] >= float(np.finfo(np.float32).tiny)
        ),
        "ee_max_abs_error": max(
            row.get("abs_error", 0) for row in rows if row["field"] == "ee_percent"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--observer", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for relative in [
        "scripts/prepare_data.py",
        "scripts/stericx_quantum.py",
        "src/reaction.rs",
        "src/kinetics/eyring.rs",
        "src/kinetics/mod.rs",
        "src/commands/simulate.rs",
        "src/commands/buried_volume.rs",
        "src/geometry/buried_volume.rs",
    ]:
        target = output / "source" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO / relative, target)
    shutil.copy2(__file__, output / "recheck.py")
    binary, observer = output / "stericx", output / "observer"
    shutil.copy2(args.binary, binary)
    shutil.copy2(args.observer, observer)
    historical = json.loads(
        (Path(__file__).parent / "historical_evidence.json").read_text()
    )
    for entry in historical["files"]:
        assert digest(REPO / entry["path"]) == entry["sha256"], entry
    python = capture_python(output)
    native = capture_native(output, binary, observer)
    summary = analyze(output)
    summary.update(
        python_cases={
            key: len(value) for key, value in python.items() if isinstance(value, list)
        },
        cli_cases=len(native),
        ni_hda=python["ni_hda"],
        status="PASS scoped numerical/domain replay",
    )
    write(output / "summary.json", summary)
    files = [
        {
            "path": str(path.relative_to(output)),
            "sha256": digest(path),
            "bytes": path.stat().st_size,
        }
        for path in sorted(output.rglob("*"))
        if path.is_file()
    ]
    write(
        output / "manifest.json",
        {
            "files": files,
            "historical_evidence_sha256": digest(
                Path(__file__).parent / "historical_evidence.json"
            ),
            "scope": (
                "Preserved audit inputs, copied Python source/executables; "
                "full build provenance is supplied by the coordinated final "
                "remediation freeze, not inferred from these source copies."
            ),
            "python_version": sys.version,
        },
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
