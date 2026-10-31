#!/usr/bin/env python3
"""Unmodified DBSTEP workflow and supplied-radii Sterimol comparisons."""
from __future__ import annotations
import argparse
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime, timezone
import io
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.spatial.transform import Rotation
from dbstep import Dbstep, calculator, constants, sterics
from morfeus import Sterimol
from study_conformer_ensembles import f32, metrics, sha, write_csv

ROOT = Path(__file__).resolve().parents[1]


def prepare_witnesses(output):
    forward = np.array([[0, 0, 1], [1.4, 0, 1.45], [-0.7, 1.212, 1.45], [-0.7, -1.212, 1.45]])
    backward = np.vstack([forward, [0, 0, -5]])
    symmetric = np.array([[0, 0, 1], [1.4, 0, 1], [-0.7, 0.7 * math.sqrt(3), 1], [-0.7, -0.7 * math.sqrt(3), 1]])
    rotation = Rotation.from_euler("xyz", [0.37, -0.53, 0.41])
    translation = np.array([2.5, -1.3, 4.0])
    permutation = [4, 2, 0, 3, 1]
    definitions = [
        ("forward", forward, np.zeros(3), list(range(4))),
        ("backward", backward, np.zeros(3), list(range(5))),
        ("symmetric", symmetric, np.zeros(3), list(range(4))),
        ("rotated_backward", rotation.apply(backward) + translation, translation, list(range(5))),
        ("permuted_backward", backward[permutation], np.zeros(3), permutation),
    ]
    witnesses = []
    for name, coordinates, center, mapping in definitions:
        directory = output / "witnesses" / name
        directory.mkdir(parents=True)
        atoms = ["P" if original == 0 else "C" for original in mapping]
        xyz = f"{len(atoms)}\nsynthetic support witness: {name}\n" + "".join(
            f"{element} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}\n" for element, point in zip(atoms, coordinates))
        (directory / "input.xyz").write_text(xyz)
        options = {"axis": "coordination", "donor_index": mapping.index(0), "center": center.tolist()}
        (directory / "options.json").write_text(json.dumps(options, indent=2) + "\n")
        witnesses.append({"id": name, "directory": str(directory), "atom_mapping_to_original": mapping,
                          "input_sha256": sha(directory / "input.xyz"), "options_sha256": sha(directory / "options.json")})
    return witnesses


def compare_case(case):
    row, options = case["row"], case["options"]
    geometry = row["geometry"]["value"]
    elements = geometry["elements"]
    coords, radii = f32(geometry["coordinates_angstrom"]), f32(geometry["radii_angstrom"])
    center = f32(row["buried_volume"]["value"]["center"])
    donor = options["donor_index"]
    extended = np.vstack([center, coords])
    reference_atoms = ["Bq", *elements]
    kwargs = dict(atom1=1, atom2=donor + 2, sterimol=True, volume=False, surface="vdw", measure="classic",
                  radii="bondi", scalevdw=1.0, noH=False, nometals=False, quiet=True)
    default = Dbstep.from_coords(reference_atoms, extended.copy(), name=case["id"], **kwargs)
    relative = extended - center
    aligned = calculator.rotate_mol(relative.copy(), 1, relative[donor + 1].copy())
    supplied = np.r_[0.0, radii]
    l, b5, b1, _ = sterics.get_classic_sterimol(aligned, supplied, reference_atoms)
    reference = {"supplied_radii_kernel": {"l": float(l), "b1": float(b1), "b5": float(b5)},
                 "default_api": {"l": float(default.L), "b1": float(default.Bmin), "b5": float(default.Bmax)}}
    m = Sterimol(["Ni", *elements], extended.copy(), dummy_index=1, attached_index=donor + 2,
                 radii=supplied, n_rot_vectors=3600)
    morfeus = {"l": float(m.L_value_uncorrected), "b1": float(m.B_1_value), "b5": float(m.B_5_value)}
    original = {key: float(f32(value)) for key, value in row["sterimol"]["value"].items()}
    native = {**original, "l": original["l"] - float(f32(0.4))}
    rel = coords - center
    axis = rel[donor] / np.linalg.norm(rel[donor])
    rho = float(np.max(np.linalg.norm(rel - np.outer(rel @ axis, axis), axis=1)))
    scale = max(1.0, float(np.max(np.linalg.norm(rel, axis=1) + radii)))
    arithmetic = 64 * 2.0**-23 * scale + 2 * math.sqrt(3) * 1e-8
    default_radii = np.array([constants.bondi[element] for element in reference_atoms])
    radii_difference = float(np.max(np.abs(default_radii - supplied)))
    comparisons, envelopes = [], {}
    for lane, values in reference.items():
        lane_radii = supplied if lane == "supplied_radii_kernel" else default_radii
        lane_coords = aligned if lane == "supplied_radii_kernel" else np.asarray(default.aligned_coords)
        positive = float(np.max(lane_coords[:, 2] + lane_radii))
        negative = float(np.max(-lane_coords[:, 2] + lane_radii))
        envelopes[lane] = {"forward": positive, "backward": negative}
        for field in ["l", "b1", "b5"]:
            common = field != "l" or positive >= negative - arithmetic
            bound = arithmetic + (rho * math.pi / 180 if field == "b1" else 0)
            if lane == "default_api": bound += radii_difference
            residual = native[field] - values[field]
            classification = "within_predeclared_bound" if common and abs(residual) <= bound else (
                "reference_convention_difference" if not common else "unresolved")
            comparisons.append({"case_id": case["id"], "source_kind": case["source_kind"], "lane": lane,
                                "descriptor": field, "native_reported": original[field], "native": native[field],
                                "reference": values[field], "residual": residual, "morfeus_common": morfeus[field],
                                "same_axial_definition": common, "arithmetic_budget": arithmetic,
                                "radius_convention_bound": radii_difference if lane == "default_api" else 0,
                                "acceptance_bound": bound, "classification": classification})
    evidence = {"case_id": case["id"], "source_kind": case["source_kind"], "source": case["source"],
                "native_reported": original, "native_compared": native, "reference": reference, "morfeus_common": morfeus,
                "dbstep_settings": kwargs, "elements": elements, "coordinates_angstrom": coords.tolist(),
                "center_angstrom": center.tolist(), "native_donor_index": donor, "dbstep_atom1": 1, "dbstep_atom2": donor + 2,
                "atom_mapping_to_original": case["atom_mapping_to_original"],
                "aligned_kernel_coordinates": aligned.tolist(), "aligned_default_coordinates": np.asarray(default.aligned_coords).tolist(),
                "supplied_radii": supplied.tolist(), "default_radii": default_radii.tolist(), "axial_envelopes": envelopes}
    return comparisons, evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output, engine = args.output.resolve(), args.engine.resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference_root = ROOT / "docs/validation/dbstep"
    inputs = ROOT / "docs/paton_recommendations/conformers/sensitivity_v1"
    manifest = json.loads((inputs / "manifest.json").read_text())
    witnesses = prepare_witnesses(output)
    receipt = {"created_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
               "script_sha256": sha(__file__), "protocol_sha256": sha(reference_root / "STERIMOL_PROTOCOL.md"),
               "reference_environment_sha256": sha(reference_root / "environment.json"), "input_manifest_sha256": sha(inputs / "manifest.json"),
               "engine": str(engine), "engine_sha256": sha(engine), "witnesses": witnesses}
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    cases, errors = [], []
    for ligand in manifest["ligands"]:
        path = inputs / ligand["ligand_id"] / "report.json"
        if sha(path) != ligand["files"]["report.json"]: raise ValueError(f"changed source report: {path}")
        report = json.loads(path.read_text())
        for row in report["descriptors"]["conformers"]:
            cases.append({"id": f"{ligand['ligand_id']}/{row['index']}", "row": row, "options": report["descriptors"]["options"],
                          "source_kind": "retained_mmff", "source": str(path.relative_to(ROOT)),
                          "atom_mapping_to_original": list(range(len(row["geometry"]["value"]["elements"])))})
    for witness in witnesses:
        directory = Path(witness["directory"])
        command = [str(engine), str(directory / "input.xyz"), str(directory / "native"),
                   "--options", str(directory / "options.json"), "--ligand-id", witness["id"]]
        process = subprocess.run(command, text=True, capture_output=True, check=False)
        (directory / "execution.json").write_text(json.dumps({"command": command, "returncode": process.returncode,
            "stdout": process.stdout, "stderr": process.stderr}, indent=2) + "\n")
        if process.returncode:
            errors.append({"case_id": witness["id"], "stage": "native", "error": process.stderr}); continue
        report = json.loads((directory / "native.json").read_text())
        cases.append({"id": witness["id"], "row": report["descriptors"]["conformers"][0],
                      "options": report["descriptors"]["options"], "source_kind": "synthetic",
                      "source": str(directory / "input.xyz"), "atom_mapping_to_original": witness["atom_mapping_to_original"]})
    comparisons, evidence = [], []
    for case in cases:
        print(f"DBSTEP {case['id']}", flush=True)
        stdout, stderr = io.StringIO(), io.StringIO()
        try:
            with redirect_stdout(stdout), redirect_stderr(stderr):
                rows, raw = compare_case(case)
            comparisons.extend(rows)
            raw.update(stdout=stdout.getvalue(), stderr=stderr.getvalue())
            evidence.append(raw)
        except (Exception, SystemExit) as error:
            errors.append({"case_id": case["id"], "stage": "reference", "error": repr(error),
                           "stdout": stdout.getvalue(), "stderr": stderr.getvalue()})
    write_csv(output / "residuals.csv", comparisons)
    (output / "case_evidence.json").write_text(json.dumps(evidence, indent=2, allow_nan=False) + "\n")
    summary = {"requested_cases": 56 + len(witnesses), "successful_cases": len(evidence), "errors": errors,
               "unresolved": [r for r in comparisons if r["classification"] == "unresolved"],
               "convention_differences": [r for r in comparisons if r["classification"] == "reference_convention_difference"],
               "metrics": {lane: {field: {
                   "all_cases": metrics([r for r in comparisons if r["lane"] == lane and r["descriptor"] == field]),
                   "common_axial_definition": metrics([r for r in comparisons if r["lane"] == lane and r["descriptor"] == field and r["same_axial_definition"]])}
                   for field in ["l", "b1", "b5"]} for lane in ["default_api", "supplied_radii_kernel"]}}
    summary["conditional_numerical_gate_pass"] = not errors and not summary["unresolved"]
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    files = {str(path.relative_to(output)): sha(path) for path in sorted(output.rglob("*")) if path.is_file()}
    (output / "manifest.json").write_text(json.dumps(files, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ["requested_cases", "successful_cases", "errors", "conditional_numerical_gate_pass"]}, indent=2))
    return not summary["conditional_numerical_gate_pass"]


if __name__ == "__main__":
    raise SystemExit(main())
