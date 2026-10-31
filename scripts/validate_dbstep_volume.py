#!/usr/bin/env python3
"""DBSTEP volume comparison: captured native lattice and independent lattices."""

from __future__ import annotations

import argparse
import copy
import io
import json
import math
import subprocess
import sys
from contextlib import redirect_stderr, redirect_stdout
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from dbstep import Dbstep, constants, sterics
from study_conformer_ensembles import f32, metrics, sha, write_csv

ROOT = Path(__file__).resolve().parents[1]


def load_cases():
    inputs = ROOT / "docs/paton_recommendations/conformers/sensitivity_v1"
    manifest = json.loads((inputs / "manifest.json").read_text())
    cases = []
    for ligand in manifest["ligands"]:
        path = inputs / ligand["ligand_id"] / "report.json"
        if sha(path) != ligand["files"]["report.json"]:
            raise ValueError(f"changed report {path}")
        report = json.loads(path.read_text())
        for row in report["descriptors"]["conformers"]:
            cases.append(
                {
                    "id": f"{ligand['ligand_id']}/{row['index']}",
                    "row": row,
                    "options": report["descriptors"]["options"],
                    "source": str(path.relative_to(ROOT)),
                    "source_sha256": sha(path),
                    "variant": "base",
                }
            )
    previous = ROOT / "docs/validation/dbstep/sterimol_v1"
    previous_manifest = json.loads((previous / "manifest.json").read_text())
    for name in [
        "forward",
        "backward",
        "symmetric",
        "rotated_backward",
        "permuted_backward",
    ]:
        relative = f"witnesses/{name}/native.json"
        path = previous / relative
        if sha(path) != previous_manifest[relative]:
            raise ValueError(f"changed witness {path}")
        report = json.loads(path.read_text())
        cases.append(
            {
                "id": name,
                "row": report["descriptors"]["conformers"][0],
                "options": report["descriptors"]["options"],
                "source": str(path.relative_to(ROOT)),
                "source_sha256": sha(path),
                "variant": "base",
            }
        )
    extras = [
        (f"SIG-NIHDA-{ligand}/0", "include_hydrogens", True)
        for ligand in [401, 498, 2067]
    ]
    extras += [("forward", "density", value) for value in [0.08, 0.00125]]
    extras += [("SIG-NIHDA-401/0", "sphere_radius", value) for value in [2.5, 4.5]]
    extras += [("SIG-NIHDA-401/0", "radii_scale", 1.0)]
    for identity, key, value in extras:
        case = copy.deepcopy(next(c for c in cases if c["id"] == identity))
        case["id"] += f"/{key}={value}"
        case["variant"] = key
        case["options"]["volume"][key] = value
        cases.append(case)
    return cases


def grid_key(config):
    return f"R={float(f32(config['sphere_radius']))};density={float(f32(config['density']))}"  # noqa: E501


def request(case, include_points):
    geometry = case["row"]["geometry"]["value"]
    details = case["row"]["buried_volume"]["value"]
    return {
        "id": case["id"],
        "op": "geometry_topology",
        "donor": case["options"]["donor_index"],
        "neighbors": [p["plane_atom_index"] for p in details["orientations"]],
        "center": f32(details["center"]).tolist(),
        "config": case["options"]["volume"],
        "dump": True,
        "include_points": include_points,
        "atoms": [
            {"element": e, "position": f32(p).tolist(), "radius": float(f32(r))}
            for e, p, r in zip(
                geometry["elements"],
                geometry["coordinates_angstrom"],
                geometry["radii_angstrom"],
                strict=True,
            )
        ],
    }


def nominal_spacing(config):
    # Only the nominal cell-volume label uses this model of the native spacing.
    # All actual reference point locations come from the captured Rust grid.
    r, density, pi = (
        np.float32(config["sphere_radius"]),
        np.float32(config["density"]),
        np.float32(math.pi),
    )
    radius_cubed = np.float32(np.float32(r * r) * r)
    volume = np.float32(
        np.float32(np.float32(np.float32(4) * pi) * radius_cubed) / np.float32(3)
    )
    side = max(
        2,
        round(
            float(
                np.cbrt(
                    np.float32(
                        np.float32(np.float32(volume / density) * np.float32(6)) / pi
                    )
                )
            )
        ),
    )
    return float(np.float32(np.float32(2) * r / np.float32(side - 1)))


def matched_case(case, captured, points):
    config = case["options"]["volume"]
    radius = float(f32(config["sphere_radius"]))
    exact_sphere_volume = 4 * math.pi * radius**3 / 3
    spacing = nominal_spacing(config)
    source_atoms = captured["atoms"]
    retained = [
        a
        for a in source_atoms
        if config["include_hydrogens"] or a["element"].upper() != "H"
    ]
    scaled = np.asarray(
        [
            np.float32(a["vdw_radius"]) * np.float32(config["radii_scale"])
            for a in retained
        ],
        dtype=np.float32,
    ).astype(float)
    rows, errors, frames = [], [], []
    if case["variant"] == "base":
        expected = case["row"]["buried_volume"]["value"]["summary"]
        for field, value in expected.items():
            if (
                np.float32(value).tobytes()
                != np.float32(captured["buried_volume"][field]).tobytes()
            ):
                errors.append(
                    {
                        "case_id": case["id"],
                        "error": "native public summary changed",
                        "field": field,
                    }
                )
    for frame in captured["dump"]["orientations"]:
        atoms = frame["aligned_atoms"]
        if len(atoms) != len(scaled):
            raise ValueError("native H-policy atom mapping differs")
        for atom, r in zip(atoms, scaled, strict=True):
            if (
                np.float32(r * r).tobytes()
                != np.float32(atom["radius_squared"]).tobytes()
            ):
                raise ValueError(
                    "scaled radius does not match captured squared-radius intermediate"
                )
        coordinates = f32([a["position"] for a in atoms])
        options = SimpleNamespace(
            grid=spacing, radius=radius, verbose=False, debug=False
        )
        occupied, tree, grid_occupied_volume = sterics.occupied(
            points, coordinates, scaled, np.zeros(3), options
        )
        percent, shell = sterics.buried_vol(
            occupied, tree, np.zeros(3), radius, 0.0, options
        )
        n_sphere = len(tree.query_ball_point(np.zeros(3), radius))
        n_occupied = int(np.count_nonzero(np.sum(occupied**2, axis=1) <= radius**2))
        if n_sphere != len(points):
            errors.append(
                {
                    "case_id": case["id"],
                    "plane": frame["plane"],
                    "error": "sphere denominator differs",
                    "native": len(points),
                    "dbstep": n_sphere,
                }
            )
        reference = percent / 100 * exact_sphere_volume
        native = float(f32(frame["buried_volume"]))
        residual = native - reference
        limit = 32 * 2.0**-23 * exact_sphere_volume
        row = {
            "case_id": case["id"],
            "variant": case["variant"],
            "plane_atom_index": frame["plane"],
            "native": native,
            "reference": reference,
            "residual": residual,
            "arithmetic_limit": limit,
            "native_sphere_points": len(points),
            "dbstep_sphere_points": n_sphere,
            "dbstep_occupied_points": n_occupied,
            "dbstep_percent": percent,
            "nominal_spacing": spacing,
            "dbstep_grid_buried_volume": n_occupied * spacing**3,
            "dbstep_sampled_sphere_volume": n_sphere * spacing**3,
            "analytic_sphere_volume": exact_sphere_volume,
            "passed": abs(residual) <= limit and n_sphere == len(points),
        }
        rows.append(row)
        frames.append(
            {
                "plane_atom_index": frame["plane"],
                "scaled_radii": scaled.tolist(),
                "retained_source_atom_indices": [
                    i
                    for i, a in enumerate(source_atoms)
                    if config["include_hydrogens"] or a["element"].upper() != "H"
                ],
                "grid_occupied_volume_returned": float(grid_occupied_volume),
                "percent_shell": float(shell),
            }
        )
    summaries = []
    for field, reference, limit in [
        (
            "buried_volume",
            math.fsum(r["reference"] for r in rows) / 3,
            32 * 2.0**-23 * exact_sphere_volume,
        ),
        (
            "percent_buried_volume",
            math.fsum(r["dbstep_percent"] for r in rows) / 3,
            32 * 2.0**-23 * 100,
        ),
    ]:
        native = float(f32(captured["buried_volume"][field]))
        summaries.append(
            {
                "case_id": case["id"],
                "variant": case["variant"],
                "descriptor": field,
                "native": native,
                "reference": reference,
                "residual": native - reference,
                "arithmetic_limit": limit,
                "passed": abs(native - reference) <= limit,
            }
        )
    return rows, summaries, errors, frames


def independent_lattice(case, captured, spacing):
    config = case["options"]["volume"]
    center = f32(captured["center"])
    coordinates = f32([a["position"] for a in captured["atoms"]])
    elements = [a["element"] for a in captured["atoms"]]
    radius, scale = (
        float(f32(config["sphere_radius"])),
        float(f32(config["radii_scale"])),
    )
    settings = dict(
        atom1=1,
        volume=True,
        sterimol=False,
        surface="vdw",
        radii="bondi",
        scalevdw=scale,
        noH=not config["include_hydrogens"],
        nometals=False,
        r=radius,
        grid=spacing,
        quiet=True,
    )
    result = Dbstep.from_coords(
        ["Bq", *elements], np.vstack([center, coordinates]), name=case["id"], **settings
    )
    if result.options.radius != radius:
        raise ValueError("DBSTEP did not apply the requested sphere radius")
    db_radii = (
        np.array([constants.bondi[element] for element in result.atoms])
        * result.options.SCALE_VDW
    )
    translated = (
        np.asarray(result.coords)
        - np.asarray(result.coords)[result.options.spec_atom_1 - 1]
    )
    bounds = sterics.max_dim(translated, db_radii, result.options)
    axes = [
        np.round(np.linspace(low, high, int(1 + round((high - low) / spacing))), 8)
        for low, high in [
            (bounds[0], bounds[1]),
            (bounds[2], bounds[3]),
            (bounds[4], bounds[5]),
        ]
    ]
    n_sphere = sterics.count_grid_points_in_sphere(*axes, np.zeros(3), radius)
    percent = float(result.bur_vol)
    occupied_estimate = percent / 100 * n_sphere
    occupied_count = round(occupied_estimate)
    if abs(occupied_estimate - occupied_count) > 1e-7:
        raise ValueError(
            "public fraction does not reconstruct an integer occupied count on the declared grid"  # noqa: E501
        )
    sphere_volume = 4 * math.pi * radius**3 / 3
    rows = []
    for field, value in [
        ("percent_buried_volume", percent),
        ("buried_volume", percent / 100 * sphere_volume),
    ]:
        native = float(f32(captured["buried_volume"][field]))
        rows.append(
            {
                "case_id": case["id"],
                "variant": case["variant"],
                "spacing_angstrom": spacing,
                "descriptor": field,
                "native": native,
                "reference": value,
                "residual": native - value,
                "same_nominal_surface": not config["include_hydrogens"],
                "earliest_convention_difference": "hydrogen_radius_and_lattice"
                if config["include_hydrogens"]
                else "floating_radius_representation_and_lattice",
            }
        )
    details = {
        "case_id": case["id"],
        "spacing": spacing,
        "settings": settings,
        "retained_elements": result.atoms.tolist(),
        "retained_coordinates": np.asarray(result.coords).tolist(),
        "scaled_radii": db_radii.tolist(),
        "grid_axes": [a.tolist() for a in axes],
        "sphere_points": n_sphere,
        "occupied_points_reconstructed_from_fraction": occupied_count,
        "percent_buried_volume": percent,
        "grid_buried_volume_derived": occupied_count * spacing**3,
        "sampled_sphere_volume": n_sphere * spacing**3,
        "whole_molecule_grid_volume_NOT_buried_volume": float(result.occ_vol),
    }
    return rows, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build", type=Path, default=ROOT / ".stericx/paton/dbstep_build_v2"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    build_path = args.build / "manifest.json"
    if sha(build_path) != build_path.with_suffix(".json.sha256").read_text().strip():
        raise ValueError("changed build manifest")
    build = json.loads(build_path.read_text())
    for item in (
        build["sources"] + build["adapter_files"] + list(build["binaries"].values())
    ):
        if sha(item["path"]) != item["sha256"]:
            raise ValueError(f"changed native source/executable {item['path']}")
    cases = load_cases()
    seen, requests = set(), []
    for case in cases:
        key = grid_key(case["options"]["volume"])
        requests.append(request(case, key not in seen))
        seen.add(key)
    raw_requests = "".join(json.dumps(r, allow_nan=False) + "\n" for r in requests)
    (output / "requests.jsonl").write_text(raw_requests)
    reference_root = ROOT / "docs/validation/dbstep"
    receipt = {
        "created_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "script_sha256": sha(__file__),
        "protocol_sha256": sha(reference_root / "BURIED_VOLUME_PROTOCOL.md"),
        "build_manifest_sha256": sha(build_path),
        "observer": build["binaries"]["observer"],
        "reference_environment_sha256": sha(reference_root / "environment.json"),
        "requests_sha256": sha(output / "requests.jsonl"),
        "cases": [
            {k: c[k] for k in ["id", "variant", "source", "source_sha256"]}
            for c in cases
        ],
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    process = subprocess.run(
        [build["binaries"]["observer"]["path"]],
        input=raw_requests,
        text=True,
        capture_output=True,
        check=False,
    )
    (output / "native_stdout.jsonl").write_text(process.stdout)
    (output / "native_stderr.txt").write_text(process.stderr)
    observed = [json.loads(line) for line in process.stdout.splitlines()]
    if process.returncode or len(observed) != len(cases):
        raise RuntimeError("incomplete native capture; retained raw streams")
    grids, frame_rows, summary_rows, api_rows, raw, errors = {}, [], [], [], [], []
    for case, captured in zip(cases, observed, strict=True):
        print(f"DBSTEP volume {case['id']}", flush=True)
        stdout, stderr = io.StringIO(), io.StringIO()
        try:
            if captured["id"] != case["id"]:
                raise ValueError("native response order/identity changed")
            key = grid_key(case["options"]["volume"])
            if captured["dump"]["points"] is not None:
                grids[key] = f32(captured["dump"]["points"])
            points = grids[key]
            if len(points) != captured["dump"]["grid_count"]:
                raise ValueError("grid count changed within fixed configuration")
            with redirect_stdout(stdout), redirect_stderr(stderr):
                frames, summaries, native_errors, reference_frames = matched_case(
                    case, captured, points
                )
                frame_rows.extend(frames)
                summary_rows.extend(summaries)
                errors.extend(native_errors)
                independent = []
                for spacing in [0.20, 0.10, 0.05]:
                    api, data = independent_lattice(case, captured, spacing)
                    api_rows.extend(api)
                    independent.append(data)
            raw.append(
                {
                    "case_id": case["id"],
                    "matched_frames": reference_frames,
                    "independent_lattices": independent,
                    "stdout": stdout.getvalue(),
                    "stderr": stderr.getvalue(),
                }
            )
        except (Exception, SystemExit) as error:
            errors.append(
                {
                    "case_id": case["id"],
                    "error": repr(error),
                    "stdout": stdout.getvalue(),
                    "stderr": stderr.getvalue(),
                }
            )
    write_csv(output / "matched_frame_residuals.csv", frame_rows)
    write_csv(output / "matched_summary_residuals.csv", summary_rows)
    write_csv(output / "independent_lattice_residuals.csv", api_rows)
    (output / "reference_raw.json").write_text(
        json.dumps(raw, indent=2, allow_nan=False) + "\n"
    )
    failed_frames = [r for r in frame_rows if not r["passed"]]
    failed_summaries = [r for r in summary_rows if not r["passed"]]
    summary = {
        "cases": len(cases),
        "matched_frames": len(frame_rows),
        "errors": errors,
        "failed_matched_frames": failed_frames,
        "failed_matched_summaries": failed_summaries,
        "matched_metrics": {
            field: metrics([r for r in summary_rows if r["descriptor"] == field])
            for field in ["buried_volume", "percent_buried_volume"]
        },
        "independent_lattice_metrics": {
            str(spacing): {
                field: metrics(
                    [
                        r
                        for r in api_rows
                        if r["spacing_angstrom"] == spacing
                        and r["descriptor"] == field
                        and r["same_nominal_surface"]
                    ]
                )
                for field in ["buried_volume", "percent_buried_volume"]
            }
            for spacing in [0.20, 0.10, 0.05]
        },
    }
    summary["finite_lattice_gate_pass"] = (
        not errors
        and not failed_frames
        and not failed_summaries
        and len(frame_rows) == 3 * len(cases)
    )
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    (output / "manifest.json").write_text(
        json.dumps(
            {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()}, indent=2
        )
        + "\n"
    )
    print(
        json.dumps(
            {
                "cases": len(cases),
                "frames": len(frame_rows),
                "errors": errors,
                "failed_matched_frames": len(failed_frames),
                "failed_matched_summaries": len(failed_summaries),
                "finite_lattice_gate_pass": summary["finite_lattice_gate_pass"],
            },
            indent=2,
        )
    )
    return not summary["finite_lattice_gate_pass"]


if __name__ == "__main__":
    raise SystemExit(main())
