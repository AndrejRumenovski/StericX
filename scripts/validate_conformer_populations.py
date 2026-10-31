#!/usr/bin/env python3
"""Independent 80-digit population reference; see the predeclared protocol."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import platform
import struct
import subprocess
import sys
from datetime import UTC, datetime
from decimal import Decimal as D
from decimal import getcontext, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = (
    ROOT / "docs/paton_recommendations/conformers/POPULATION_VALIDATION_PROTOCOL.md"
)
ABS_LIMIT = 32 * sys.float_info.epsilon
REL_LIMIT = 2e-12


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact_float(bits):
    return D.from_float(struct.unpack(">d", struct.pack(">Q", bits))[0])


def constants():
    na = D("6.02214076e23")
    return D("1.380649e-23") * na, {
        "kcal_mol": D(4184),
        "kj_mol": D(1000),
        "hartree": D("4.359744722206e-18") * na,
        "electron_volt": D("1.602176634e-19") * na,
    }


def request(
    name, values, unit="kcal_mol", temperature=298.15, kind="electronic", gs=None
):
    return {
        "id": name,
        "temperature_k": temperature,
        "assumptions": "analytical synthetic states; g counts only unenumerated multiplicity",  # noqa: E501
        "energies": [
            {
                "value": float(value),
                "unit": unit,
                "kind": kind,
                "model": "synthetic mathematical validation",
                "provenance": f"{name}/state_{i}",
                "temperature_k": temperature if kind == "gibbs_free_energy" else None,
                "degeneracy": gs[i] if gs else 1,
            }
            for i, value in enumerate(values)
        ],
    }


def cases():
    valid, invalid = [], []
    with localcontext() as ctx:
        ctx.prec = 80
        r, conversions = constants()
        for unit, conversion in conversions.items():
            for temperature in [1.0, 77.0, 298.15, 1000.0, 1e6]:
                gap = r * D.from_float(temperature) * D(2).ln() / conversion
                valid.append(
                    request(
                        f"two_state_{unit}_{temperature}", [0, gap], unit, temperature
                    )
                )
            values = [D(v) * D(4184) / conversion for v in ["0", "0.1", "1", "3"]]
            valid.append(request(f"multistate_{unit}", values, unit))
            valid.append(request(f"permuted_{unit}", values[::-1], unit))
            offset = D(1000) if unit != "hartree" else D(-1000)
            valid.append(request(f"offset_{unit}", [v + offset for v in values], unit))
        for kind in ["electronic", "force_field", "gibbs_free_energy"]:
            valid.append(request(f"kind_{kind}", [0, 0.5, 2.0], kind=kind))
        valid.extend(
            [
                request("multiplicities", [0, 0, 0], gs=[1, 2, 4]),
                request("explicit_equal_states", [0] * 7),
                request("large_multiplicity", [0, 0, 1], gs=[1, 4294967295, 3]),
                request(
                    "cold_subnormal",
                    [0, 1],
                    temperature=float.fromhex("0x0.0000000000001p-1022"),
                ),
                request("large_gap_underflow", [0, 1e300]),
                request(
                    "opposite_extreme_energies",
                    [-sys.float_info.max, sys.float_info.max],
                    "kj_mol",
                    sys.float_info.max,
                ),
                request(
                    "large_ensemble",
                    [(i % 37) / 10.0 for i in range(4096)],
                    gs=[1 + i % 4 for i in range(4096)],
                ),
            ]
        )
    for label, field, value in [
        ("zero_temperature", "temperature_k", 0),
        ("negative_temperature", "temperature_k", -1),
        ("empty_assumptions", "assumptions", ""),
        ("empty_states", "energies", []),
    ]:
        case = request(label, [0, 1])
        case[field] = value
        invalid.append(case)
    for field, value in [
        ("unit", "hartree"),
        ("kind", "force_field"),
        ("model", "other"),
        ("provenance", ""),
        ("degeneracy", 0),
    ]:
        case = request(f"invalid_{field}", [0, 1])
        case["energies"][1][field] = value
        invalid.append(case)
    mismatch = request("gibbs_temperature_mismatch", [0, 1], kind="gibbs_free_energy")
    mismatch["temperature_k"] = 300.0
    invalid.append(mismatch)
    missing = copy.deepcopy(mismatch)
    missing["id"] = "gibbs_missing_temperature"
    missing["energies"][0]["temperature_k"] = None
    invalid.append(missing)
    return valid, invalid


def reference(case, response):
    with localcontext() as ctx:
        ctx.prec = 80
        r, conversions = constants()
        values = [exact_float(bits) for bits in response["energy_bits"]]
        temperature = exact_float(response["temperature_bits"])
        minimum = min(values)
        terms = []
        for value, state in zip(values, case["energies"], strict=True):
            exponent = (
                -(value - minimum) * conversions[state["unit"]] / (r * temperature)
            )
            terms.append(
                D(state["degeneracy"])
                * (exponent.exp() if exponent >= -10000 else D(0))
            )
        total = sum(terms)
        return [term / total for term in terms]


def main():
    getcontext().prec = 80
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    observer = args.observer.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    sources = [
        PROTOCOL,
        Path(__file__),
        ROOT / "src/ensemble/population.rs",
        ROOT / "src/ensemble/mod.rs",
        ROOT / "src/kinetics/eyring.rs",
        ROOT / "examples/population_observer.rs",
        ROOT / "Cargo.lock",
    ]
    # This receipt is written before executing the observer.
    receipt = {
        "started_utc": datetime.now(UTC).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "command": sys.argv,
        "observer": str(observer),
        "observer_sha256": sha(observer),
        "head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "sources": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "absolute_limit": ABS_LIMIT,
        "relative_limit": REL_LIMIT,
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    valid, invalid = cases()
    requests = valid + invalid
    raw = "".join(json.dumps(case, allow_nan=False) + "\n" for case in requests)
    (args.output / "requests.jsonl").write_text(raw)
    process = subprocess.run(
        [str(observer)], input=raw, text=True, capture_output=True, check=False
    )
    (args.output / "stdout.jsonl").write_text(process.stdout)
    (args.output / "stderr.txt").write_text(process.stderr)
    responses = [json.loads(line) for line in process.stdout.splitlines()]
    failures, rows, max_abs, max_rel = [], [], 0.0, 0.0
    if process.returncode or len(responses) != len(requests):
        failures.append(
            f"observer exit/count: {process.returncode}, {len(responses)}/{len(requests)}"  # noqa: E501
        )
    for index, (case, response) in enumerate(zip(requests, responses, strict=True)):
        expected_success = index < len(valid)
        row = {
            "id": case["id"],
            "expected_success": expected_success,
            "response_status": response["status"],
        }
        if (
            response.get("id") != case["id"]
            or (response["status"] == "success") != expected_success
        ):
            failures.append(f"{case['id']}: unexpected outcome")
        elif expected_success:
            expected = reference(case, response)
            actual = response["population"]["weights"]
            errors = [
                abs(D.from_float(a) - e) for a, e in zip(actual, expected, strict=True)
            ]
            relative = [
                error / e
                for error, e in zip(errors, expected, strict=True)
                if e >= D("1e-300")
            ]
            absolute_error = float(max(errors))
            relative_error = float(max(relative, default=D(0)))
            normalization_error = abs(sum(D.from_float(a) for a in actual) - D(1))
            requested_bits = [
                struct.unpack(">Q", struct.pack(">d", float(s["value"])))[0]
                for s in case["energies"]
            ]
            parse_differences = [
                {"state": i, "requested_bits": before, "accepted_bits": after}
                for i, (before, after) in enumerate(
                    zip(requested_bits, response["energy_bits"], strict=True)
                )
                if before != after
            ]
            valid_output = len(actual) == len(expected) and all(
                a >= 0 and a < float("inf") for a in actual
            )
            passed = (
                valid_output
                and absolute_error <= ABS_LIMIT
                and relative_error <= REL_LIMIT
                and normalization_error <= D.from_float(ABS_LIMIT)
            )
            row.update(
                reference_weights=[str(e) for e in expected],
                max_absolute_error=absolute_error,
                max_relative_error=relative_error,
                normalization_error=str(normalization_error),
                json_input_rounding_differences=parse_differences,
                passed=passed,
            )
            if not passed:
                failures.append(f"{case['id']}: numerical gate failed")
            max_abs = max(max_abs, absolute_error)
            max_rel = max(max_rel, relative_error)
        rows.append(row)
    report = {
        "valid_cases": len(valid),
        "invalid_cases": len(invalid),
        "returncode": process.returncode,
        "max_absolute_error": max_abs,
        "max_relative_error": max_rel,
        "failures": failures,
        "cases": rows,
        "files": {
            name: sha(args.output / name)
            for name in ["requests.jsonl", "stdout.jsonl", "stderr.txt", "receipt.json"]
        },
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in report.items() if k not in ["cases", "files"]}, indent=2
        )
    )
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
