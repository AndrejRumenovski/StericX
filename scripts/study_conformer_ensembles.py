#!/usr/bin/env python3
"""Fixed-cohort ensemble behavior study with unmodified Morfeus references."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import math
import struct
import sys
from collections import defaultdict
from contextlib import redirect_stderr, redirect_stdout
from datetime import UTC, datetime
from pathlib import Path

import morfeus
import numpy as np
from morfeus import BuriedVolume, Pyramidalization, Sterimol
from validate_conformer_populations import reference as population_reference

ROOT = Path(__file__).resolve().parents[1]
FIELDS = [
    "sterimol_l",
    "sterimol_b1",
    "sterimol_b5",
    "buried_volume",
    "percent_buried_volume",
    "qvbur_min",
    "qvbur_max",
    "max_delta_qvbur",
    "ovbur_min",
    "ovbur_max",
    "near_vbur",
    "far_vbur",
    "pyr_p",
    "pyr_alpha",
]
MODES = {
    "single": None,
    "arithmetic": "mean",
    "boltzmann": "boltzmann_mean",
    "minimum": "minimum",
    "maximum": "maximum",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def f32(value):
    return np.asarray(value, dtype=np.float32).astype(float)


def native_values(row):
    out = dict.fromkeys(FIELDS)
    for family, fields in [
        ("sterimol", {"l": "sterimol_l", "b1": "sterimol_b1", "b5": "sterimol_b5"}),
        ("pyramidalization", {"pyr_p": "pyr_p", "pyr_alpha": "pyr_alpha"}),
        ("buried_volume", {field: field for field in FIELDS[3:12]}),
    ]:
        outcome = row[family]
        if outcome["status"] == "success":
            values = (
                outcome["value"]["summary"]
                if family == "buried_volume"
                else outcome["value"]
            )
            out.update({name: float(f32(values[key])) for key, name in fields.items()})
    return out


def references(row, options):
    result, failures, frames = dict.fromkeys(FIELDS), {}, []
    if (
        row["geometry"]["status"] != "success"
        or row["buried_volume"]["status"] != "success"
    ):
        return (
            result,
            {"shared_inputs": "source geometry or common center unavailable"},
            frames,
        )
    geometry = row["geometry"]["value"]
    coordinates, radii = (
        f32(geometry["coordinates_angstrom"]),
        f32(geometry["radii_angstrom"]),
    )
    elements = geometry["elements"]
    details = row["buried_volume"]["value"]
    center = f32(details["center"])
    donor = options["donor_index"]
    neighbors = [frame["plane_atom_index"] for frame in details["orientations"]]
    extended = np.vstack([center, coordinates])
    extended_elements = ["Ni", *elements]
    try:
        s = Sterimol(
            extended_elements,
            extended.copy(),
            dummy_index=1,
            attached_index=donor + 2,
            radii=np.r_[0.0, radii],
            n_rot_vectors=3600,
        )
        result.update(
            sterimol_l=float(s.L_value),
            sterimol_b1=float(s.B_1_value),
            sterimol_b5=float(s.B_5_value),
        )
    except Exception as error:
        failures["sterimol"] = repr(error)
    try:
        p = Pyramidalization(
            coordinates, donor + 1, neighbor_indices=[i + 1 for i in neighbors]
        )
        result.update(pyr_p=float(p.P), pyr_alpha=float(p.alpha))
    except Exception as error:
        failures["pyramidalization"] = repr(error)
    try:
        config = options["volume"]
        for neighbor in neighbors:
            b = BuriedVolume(
                extended_elements,
                extended.copy(),
                metal_index=1,
                excluded_atoms=[1],
                radii=np.r_[0.0, radii * float(f32(config["radii_scale"]))],
                include_hs=config["include_hydrogens"],
                radius=float(f32(config["sphere_radius"])),
                density=float(f32(config["density"])),
                z_axis_atoms=[donor + 2],
                xz_plane_atoms=[neighbor + 2],
            )
            b.octant_analysis()
            q = [float(b.quadrants["buried_volume"][i]) for i in [1, 2, 3, 4]]
            o = {str(i): float(v) for i, v in b.octants["buried_volume"].items()}
            frames.append(
                {
                    "plane_atom_index": neighbor,
                    "buried_volume": float(b.buried_volume),
                    "percent": float(b.fraction_buried_volume * 100),
                    "quadrants": q,
                    "octants": o,
                    "near": math.fsum(o[str(i)] for i in [4, 5, 6, 7]),
                    "far": math.fsum(o[str(i)] for i in [0, 1, 2, 3]),
                }
            )
        qs = [q for frame in frames for q in frame["quadrants"]]
        os = [v for frame in frames for v in frame["octants"].values()]
        result.update(
            buried_volume=math.fsum(f["buried_volume"] for f in frames) / 3,
            percent_buried_volume=math.fsum(f["percent"] for f in frames) / 3,
            near_vbur=math.fsum(f["near"] for f in frames) / 3,
            far_vbur=math.fsum(f["far"] for f in frames) / 3,
            qvbur_min=min(qs),
            qvbur_max=max(qs),
            ovbur_min=min(os),
            ovbur_max=max(os),
            max_delta_qvbur=max(
                abs(f["quadrants"][i] - f["quadrants"][i - 1])
                for f in frames
                for i in range(4)
            ),
        )
    except Exception as error:
        failures["buried_volume"] = repr(error)
    for field, value in result.items():
        if value is not None and not math.isfinite(value):
            failures[field] = "nonfinite reference result"
            result[field] = None
    return result, failures, frames


def tolerances(row, options):
    if (
        row["geometry"]["status"] != "success"
        or row["buried_volume"]["status"] != "success"
    ):
        return dict.fromkeys(FIELDS)
    g = row["geometry"]["value"]
    center = f32(row["buried_volume"]["value"]["center"])
    coordinates, radii = f32(g["coordinates_angstrom"]), f32(g["radii_angstrom"])
    rel = coordinates - center
    scale = max(1.0, float(np.max(np.linalg.norm(rel, axis=1) + radii)))
    axis = rel[options["donor_index"]]
    axis = axis / np.linalg.norm(axis)
    rho = float(np.max(np.linalg.norm(rel - np.outer(rel @ axis, axis), axis=1)))
    eps = 2.0**-23
    volume = 4 * math.pi * float(f32(options["volume"]["sphere_radius"])) ** 3 / 3
    limits = {f: 32 * eps * volume for f in FIELDS}
    limits.update(
        sterimol_l=64 * eps * scale,
        sterimol_b5=64 * eps * scale,
        sterimol_b1=64 * eps * scale + rho * (math.pi / 360 + math.pi / 3599),
        percent_buried_volume=32 * eps * 100,
        pyr_p=1e-5,
        pyr_alpha=1e-3,
    )
    return limits


def metrics(rows):
    pairs = [
        (r["reference"], r["native"])
        for r in rows
        if r["reference"] is not None and r["native"] is not None
    ]
    result = {"n": len(pairs), "missing_pairs": len(rows) - len(pairs)}
    if not pairs:
        return result
    x, y = np.asarray(pairs, dtype=float).T
    residual = y - x
    ss = float(np.sum((x - x.mean()) ** 2))
    slope = (
        float(np.sum((x - x.mean()) * (y - y.mean())) / ss)
        if ss > 0 and len(x) > 1
        else None
    )
    result.update(
        mae=float(np.mean(np.abs(residual))),
        rmse=float(np.sqrt(np.mean(residual**2))),
        maximum_absolute_error=float(np.max(np.abs(residual))),
        median_absolute_error=float(np.median(np.abs(residual))),
        r_squared=1 - float(np.sum(residual**2)) / ss if ss > 0 else None,
        slope=slope,
        intercept=float(y.mean() - slope * x.mean()) if slope is not None else None,
    )
    return result


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=ROOT / "docs/paton_recommendations/conformers/sensitivity_v1",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    input_manifest = json.loads((args.inputs / "manifest.json").read_text())
    protocol = ROOT / "docs/paton_recommendations/conformers/ENSEMBLE_STUDY_PROTOCOL.md"
    historical_path = ROOT / "data/official/ni_hda_kraken.csv"
    reference_root = Path(morfeus.__file__).parent
    receipt = {
        "created_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "python": sys.version,
        "packages": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "scipy", "morfeus-ml"]
        },
        "protocol_sha256": sha(protocol),
        "script_sha256": sha(__file__),
        "population_reference_sha256": sha(
            Path(__file__).with_name("validate_conformer_populations.py")
        ),
        "input_manifest_sha256": sha(args.inputs / "manifest.json"),
        "historical_source_sha256": sha(historical_path),
        "reference_sources": {
            name: sha(reference_root / name)
            for name in [
                "buried_volume.py",
                "sterimol.py",
                "pyramidalization.py",
                "geometry.py",
                "utils.py",
            ]
        },
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    official = defaultdict(list)
    for row in csv.DictReader(historical_path.open()):
        official[row[""]].append(row)
    all_states, all_reductions, historical, raw, failures, investigations = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for ligand in input_manifest["ligands"]:
        identity = ligand["ligand_id"]
        print(f"Reference {identity}", flush=True)
        path = args.inputs / identity / "report.json"
        if sha(path) != ligand["files"]["report.json"]:
            raise ValueError(f"changed input report {identity}")
        report = json.loads(path.read_text())
        options = report["descriptors"]["options"]
        rows = report["descriptors"]["conformers"]
        energy_origin = report["sensitivity"]["statistics"]["population"]["origin"]

        def bits(value):
            return struct.unpack(">Q", struct.pack(">d", value))[0]

        weights = [
            float(w)
            for w in population_reference(
                {"energies": energy_origin["energies"]},
                {
                    "energy_bits": [
                        bits(e["value"]) for e in energy_origin["energies"]
                    ],
                    "temperature_bits": bits(energy_origin["temperature_k"]),
                },
            )
        ]
        reference_rows, native_rows, limits = [], [], []
        for row in rows:
            stdout, stderr = io.StringIO(), io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                ref, errors, frames = references(row, options)
            native, limit = native_values(row), tolerances(row, options)
            reference_rows.append(ref)
            native_rows.append(native)
            limits.append(limit)
            raw.append(
                {
                    "ligand_id": identity,
                    "index": row["index"],
                    "reference": ref,
                    "errors": errors,
                    "volume_frames": frames,
                    "stdout": stdout.getvalue(),
                    "stderr": stderr.getvalue(),
                }
            )
            if errors:
                failures.append(
                    {"ligand_id": identity, "index": row["index"], "errors": errors}
                )
            for field in FIELDS:
                residual = (
                    native[field] - ref[field]
                    if native[field] is not None and ref[field] is not None
                    else None
                )
                investigate = residual is not None and abs(residual) > limit[field]
                item = {
                    "ligand_id": identity,
                    "index": row["index"],
                    "descriptor": field,
                    "native": native[field],
                    "reference": ref[field],
                    "residual": residual,
                    "triage_limit": limit[field],
                    "investigation_required": investigate,
                }
                all_states.append(item)
                if investigate:
                    investigations.append(item)
        for field in FIELDS:
            values = [r[field] for r in reference_rows]
            valid = all(v is not None for v in values)
            reductions = {
                "single": values[0],
                "arithmetic": math.fsum(values) / len(values) if valid else None,
                "boltzmann": math.fsum(
                    v * w for v, w in zip(values, weights, strict=True)
                )
                if valid
                else None,
                "minimum": min(values) if valid else None,
                "maximum": max(values) if valid else None,
            }
            for mode, key in MODES.items():
                native = (
                    native_rows[0][field]
                    if key is None
                    else report["sensitivity"]["statistics"]["descriptors"][field][key]
                )
                reference = reductions[mode]
                bound = max(
                    (limit[field] for limit in limits if limit[field] is not None),
                    default=None,
                )
                residual = (
                    native - reference
                    if native is not None and reference is not None
                    else None
                )
                all_reductions.append(
                    {
                        "ligand_id": identity,
                        "descriptor": field,
                        "mode": mode,
                        "native": native,
                        "reference": reference,
                        "residual": residual,
                        "triage_limit": bound,
                        "investigation_required": residual is not None
                        and abs(residual) > bound + 1e-12,
                    }
                )
                if field == "max_delta_qvbur":
                    target_rows = official[identity.removeprefix("SIG-NIHDA-")]
                    target = (
                        float(target_rows[0]["vbur_max_delta_qvbur_min"])
                        if len(target_rows) == 1
                        else None
                    )
                    historical.append(
                        {
                            "ligand_id": identity,
                            "mode": mode,
                            "native": native,
                            "reference": target,
                            "residual": native - target
                            if native is not None and target is not None
                            else None,
                            "same_reduction_definition": mode == "minimum",
                            "comparison_kind": "historical_workflow_estimator_not_identical_geometry",  # noqa: E501
                        }
                    )
    write_csv(args.output / "conformer_residuals.csv", all_states)
    write_csv(args.output / "ensemble_residuals.csv", all_reductions)
    write_csv(args.output / "historical_residuals.csv", historical)
    (args.output / "reference_raw.json").write_text(
        json.dumps(raw, indent=2, allow_nan=False) + "\n"
    )
    summary = {
        "ligands": len(input_manifest["ligands"]),
        "conformers": len(raw),
        "reference_failures": failures,
        "conformer_investigations": investigations,
        "ensemble_investigations": [
            r for r in all_reductions if r["investigation_required"]
        ],
        "independent_metrics": {
            field: {
                mode: metrics(
                    [
                        r
                        for r in all_reductions
                        if r["descriptor"] == field and r["mode"] == mode
                    ]
                )
                for mode in MODES
            }
            for field in FIELDS
        },
        "historical_metrics": {
            mode: metrics([r for r in historical if r["mode"] == mode])
            for mode in MODES
        },
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    (args.output / "manifest.json").write_text(
        json.dumps(
            {
                path.name: sha(path)
                for path in sorted(args.output.iterdir())
                if path.is_file()
            },
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            {
                "ligands": summary["ligands"],
                "conformers": summary["conformers"],
                "failed_reference_calls": len(failures),
                "conformer_investigations": len(investigations),
                "ensemble_investigations": len(summary["ensemble_investigations"]),
                "historical_metrics": summary["historical_metrics"],
            },
            indent=2,
        )
    )
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
