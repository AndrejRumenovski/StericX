#!/usr/bin/env python3
"""Retain point/atom witnesses for matched-table f32/f64 occupancy differences.

This diagnostic reconstructs the documented f32 predicate from captured aligned
atoms and obtains the actual grid from the frozen observer. The independent
float64 frame/grid/union equations remain unchanged. It adds no acceptance bound.
"""

import argparse
import itertools
import json
import math
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from sealed_reference import ROOT, load_reference, refuse_audit_output, sha

reference = load_reference()


def rows(path):
    with path.open() as handle:
        for line in handle:
            yield json.loads(line)


def same_f32(left, right):
    return struct.pack("f", left) == struct.pack("f", right)


def main(args):
    refuse_audit_output(args.output)
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.output / Path(__file__).name)
    obs_manifest = args.observations / "manifest.json"
    manifest = json.loads(obs_manifest.read_text())
    if (
        sha(obs_manifest)
        != obs_manifest.with_suffix(".json.sha256").read_text().strip()
    ):
        raise ValueError("Changed observation seal")
    lane = manifest["lanes"]["kraken_topology_all_bins"]
    build_path = Path(manifest["build_manifest"]["path"])
    build = json.loads(build_path.read_text())
    binary_record = build["binaries"]["observer"]
    binary = Path(binary_record["path"])
    if sha(binary) != binary_record["sha256"]:
        raise ValueError("Changed frozen observer")
    receipt = json.loads((args.references / "complete.json").read_text())
    ref_path = args.references / "reference_results.jsonl"
    if sha(ref_path) != receipt["reference_results_sha256"]:
        raise ValueError("Changed references")
    for rec in [manifest["build_manifest"], lane["requests"], lane["stdout"]]:
        if sha(Path(rec["path"])) != rec["sha256"]:
            raise ValueError("Changed captured input")
    inputs = {
        str(p): sha(p)
        for p in [
            obs_manifest,
            build_path,
            binary,
            Path(lane["requests"]["path"]),
            Path(lane["stdout"]["path"]),
            ref_path,
            Path(__file__),
            ROOT / "docs/scientific_accuracy_audit/scripts/geometry_reference.py",
        ]
    }
    selected = []
    n_all = 0
    for req, obs, ref in itertools.zip_longest(
        rows(Path(lane["requests"]["path"])),
        rows(Path(lane["stdout"]["path"])),
        rows(ref_path),
    ):
        if (
            any(x is None for x in [req, obs, ref])
            or not req["id"] == obs["id"] == ref["id"]
        ):
            raise ValueError("Mismatched full-corpus IDs")
        n_all += 1
        expected_radii = reference.radii_for(req)
        if any(
            not same_f32(a["vdw_radius"], r)
            for a, r in zip(obs["atoms"], expected_radii, strict=True)
        ):
            continue
        expected_frames = {
            x["plane"]: x["independent"] for x in ref["buried_volume"]["orientations"]
        }
        changed = []
        sphere = 4 * math.pi * req["config"]["sphere_radius"] ** 3 / 3
        for frame in obs["dump"]["orientations"]:
            expected = expected_frames[frame["plane"]]
            inferred = [
                round(v / (sphere / 8) * total)
                for v, total in zip(
                    frame["octants"], expected["octant_grid_counts"], strict=True
                )
            ]
            if inferred != expected["octant_occupied_counts"]:
                changed.append(frame["plane"])
        if changed:
            selected.append((req, obs, ref, changed))
    if n_all != 31721 or not selected:
        raise ValueError("Incomplete or empty boundary selection")
    request = dict(selected[0][0], include_points=True)
    request_path = args.output / "grid_request.jsonl"
    request_path.write_text(json.dumps(request) + "\n")
    with (
        request_path.open("rb") as source,
        (args.output / "grid_observation.jsonl").open("xb") as out,
        (args.output / "grid_stderr.txt").open("xb") as err,
    ):
        proc = subprocess.run(
            [str(binary)], stdin=source, stdout=out, stderr=err, check=False
        )
    if proc.returncode:
        raise ValueError("Grid observer failed")
    grid_observation = json.loads((args.output / "grid_observation.jsonl").read_text())
    grid_record = dict(grid_observation)
    grid_record["dump"] = dict(grid_record["dump"], points=None)
    if grid_record != selected[0][1]:
        raise ValueError("Frozen observer changed an existing scientific observation")
    native_points = np.asarray(grid_observation["dump"]["points"], dtype=np.float32)
    cloud = reference.grid(3.5, 0.01)
    if native_points.shape != cloud.shape or len(cloud) != 15408:
        raise ValueError("Unexpected default lattice")
    # Same ordered integer lattice sites; differences are represented coordinates.
    if not np.array_equal(np.signbit(native_points), np.signbit(cloud)):
        raise ValueError("Different regional point membership")
    sphere32 = np.float32(
        np.float32(np.float32(4) * np.float32(np.pi))
        * np.float32(np.float32(3.5) ** 3)
        / np.float32(3)
    )
    witnesses = []
    frames_reconstructed = points_changed = 0
    with (
        (args.output / "selected_inputs.jsonl").open("x") as pair_file,
        (args.output / "point_atom_witnesses.jsonl").open("x") as output,
    ):
        for req, obs, ref, changed in selected:
            pair_file.write(
                json.dumps({"request": req, "observation": obs, "reference": ref})
                + "\n"
            )
            xyz = np.array(
                [atom["position"] for atom in obs["atoms"]], dtype=np.float64
            )
            center = np.array(obs["center"])
            radii = reference.radii_for(req) * req["config"]["radii_scale"]
            include = np.array([a["element"].upper() != "H" for a in req["atoms"]])
            indices = np.flatnonzero(include)
            z = reference.unit(center - xyz[req["donor"]])
            expected_frames = {
                x["plane"]: x["independent"]
                for x in ref["buried_volume"]["orientations"]
            }
            for frame in obs["dump"]["orientations"]:
                if frame["plane"] not in changed:
                    continue
                p = xyz[frame["plane"]] - center
                x = reference.unit(p - np.dot(p, z) * z)
                aligned = (xyz - center) @ np.array([x, np.cross(z, x), z]).T
                native_atoms = frame["aligned_atoms"]
                if len(native_atoms) != len(indices):
                    raise ValueError("Different atom-inclusion policy")
                native_diffs, reference_diffs = [], []
                for atom, atom_index in zip(native_atoms, indices, strict=True):
                    delta32 = native_points - np.array(
                        atom["position"], dtype=np.float32
                    )
                    d32 = (
                        delta32[:, 0] * delta32[:, 0] + delta32[:, 1] * delta32[:, 1]
                    ) + delta32[:, 2] * delta32[:, 2]
                    native_diffs.append(d32.astype(float) - atom["radius_squared"])
                    delta64 = cloud - aligned[atom_index]
                    reference_diffs.append(
                        np.sum(delta64 * delta64, axis=1) - radii[atom_index] ** 2
                    )
                native_diffs, reference_diffs = (
                    np.array(native_diffs),
                    np.array(reference_diffs),
                )
                native_mask = native_diffs.min(axis=0) <= 0
                independent_mask = reference_diffs.min(axis=0) <= 0
                # Independently require cKDTree union agreement, the sealed kernel.
                tree_mask = np.zeros(len(cloud), bool)
                tree = cKDTree(cloud)
                for atom_index in indices:
                    tree_mask[
                        tree.query_ball_point(aligned[atom_index], radii[atom_index])
                    ] = True
                if not np.array_equal(tree_mask, independent_mask):
                    raise ValueError(
                        "Analytic witness differs from unchanged cKDTree reference"
                    )
                expected = expected_frames[frame["plane"]]
                if independent_mask.sum() != expected["occupied_total"]:
                    raise ValueError(
                        "Independent reconstruction differs from sealed reference"
                    )
                observed_v = (
                    np.float32(np.float32(native_mask.sum()) / np.float32(len(cloud)))
                    * sphere32
                )
                if not same_f32(observed_v, frame["buried_volume"]):
                    raise ValueError(
                        "Native reconstruction differs from captured volume"
                    )
                counts = []
                for zi in [1, -1]:
                    for xi, yi in [(1, 1), (-1, 1), (-1, -1), (1, -1)]:
                        bin_mask = (
                            (native_points[:, 0] * xi > 0)
                            & (native_points[:, 1] * yi > 0)
                            & (native_points[:, 2] * zi > 0)
                        )
                        count = int((native_mask & bin_mask).sum())
                        counts.append(count)
                        volume = np.float32(
                            np.float32(np.float32(count) / np.float32(bin_mask.sum()))
                            * sphere32
                        ) / np.float32(8)
                        if not same_f32(volume, frame["octants"][len(counts) - 1]):
                            raise ValueError(
                                "Native reconstruction differs from captured bin"
                            )
                flipped = np.flatnonzero(native_mask != independent_mask)
                frames_reconstructed += 1
                points_changed += len(flipped)
                witnesses.append(
                    {
                        "id": req["id"],
                        "plane": frame["plane"],
                        "native_octant_counts": counts,
                        "reference_octant_counts": expected["octant_occupied_counts"],
                        "changed_points": len(flipped),
                    }
                )
                for point in flipped:
                    native_nearest = int(native_diffs[:, point].argmin())
                    ref_nearest = int(reference_diffs[:, point].argmin())
                    involved = sorted({native_nearest, ref_nearest})
                    record = {
                        "id": req["id"],
                        "plane": frame["plane"],
                        "point_index": int(point),
                        "native_point": native_points[point].tolist(),
                        "reference_point": cloud[point].tolist(),
                        "native_occupied": bool(native_mask[point]),
                        "reference_occupied": bool(independent_mask[point]),
                        "witness_atoms": [
                            {
                                "atom_index": int(indices[index]),
                                "native_aligned": native_atoms[index],
                                "reference_aligned_position": aligned[
                                    indices[index]
                                ].tolist(),
                                "reference_radius": float(radii[indices[index]]),
                                "native_squared_margin": float(
                                    native_diffs[index, point]
                                ),
                                "reference_squared_margin": float(
                                    reference_diffs[index, point]
                                ),
                                "margin_change": float(
                                    native_diffs[index, point]
                                    - reference_diffs[index, point]
                                ),
                            }
                            for index in involved
                        ],
                    }
                    output.write(json.dumps(record, allow_nan=False) + "\n")
    if any(sha(Path(p)) != digest for p, digest in inputs.items()):
        raise ValueError("Changed input during diagnostic")
    summary = {
        "kind": "minimized_point_atom_boundary_witnesses",
        "global_scientific_pass": False,
        "threshold_added": False,
        "inputs": inputs,
        "full_corpus_scanned": n_all,
        "selected_same_table_rows": len(selected),
        "frames_reconstructed": frames_reconstructed,
        "changed_point_decisions": points_changed,
        "exact_native_total_and_octant_bits_reconstructed": True,
        "unchanged_independent_occupied_counts_reconstructed": True,
        "frames": witnesses,
        "interpretation": (
            "Discrete f32/f64 membership differences reproduced from actual "
            "captured grid/aligned atoms versus unchanged independent grid/frame/"
            "radii. Point/atom margins retain each sign flip. These are arithmetic "
            "boundary witnesses, not an empirical global error allowance."
        ),
        "artifacts": {p.name: sha(p) for p in args.output.iterdir() if p.is_file()},
    }
    (args.output / "complete.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["observations", "references", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    main(parser.parse_args())
