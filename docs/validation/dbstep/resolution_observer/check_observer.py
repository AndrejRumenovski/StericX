#!/usr/bin/env python3
"""Verify observer integrity on preserved 401/0 R=2.5; not a chemistry oracle."""

import argparse
import hashlib
import json
import struct
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
FIELDS = ("buried_volume", "quadrants", "octants", "near_vbur", "far_vbur")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def f32(value):
    return struct.unpack("<f", struct.pack("<f", value))[0]


def bits(value):
    if isinstance(value, list):
        return [bits(x) for x in value]
    return struct.unpack("<I", struct.pack("<f", value))[0]


def invoke(binary, requests, output, label):
    payload = "".join(json.dumps(r, separators=(",", ":")) + "\n" for r in requests)
    (output / (label + "_requests.jsonl")).write_text(payload)
    run = subprocess.run(
        [str(binary.resolve())],
        input=payload,
        text=True,
        capture_output=True,
        check=False,
    )
    (output / (label + "_stdout.jsonl")).write_text(run.stdout)
    (output / (label + "_stderr.txt")).write_text(run.stderr)
    assert run.returncode == 0, (label, run.stderr)
    results = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(results) == len(requests), label
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-binary", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    original = ROOT / "docs/validation/dbstep/volume_v2"
    case_id = "SIG-NIHDA-401/0/sphere_radius=2.5"
    request = next(
        json.loads(line)
        for line in (original / "requests.jsonl").read_text().splitlines()
        if json.loads(line)["id"] == case_id
    )
    historical = next(
        json.loads(line)
        for line in (original / "native_stdout.jsonl").read_text().splitlines()
        if json.loads(line)["id"] == case_id
    )
    request["include_points"] = True
    base = invoke(args.baseline_binary, [request], output, "baseline_geometry")[0]
    observed = invoke(args.binary, [request], output, "resolution_geometry")[0]
    assert base == observed, "legacy geometry operation changed"
    for key in historical:
        if key != "dump":
            assert historical[key] == base[key], key
    assert historical["dump"]["orientations"] == base["dump"]["orientations"]
    points = base["dump"]["points"]
    requests = []
    raw_requests = []
    for frame in base["dump"]["orientations"]:
        common = {
            "atoms": frame["aligned_atoms"],
            "config": request["config"],
            "include_points": True,
        }
        requests.extend(
            [
                {**common, "id": str(frame["plane"]) + "/grid", "op": "grid_occupancy"},
                {
                    **common,
                    "id": str(frame["plane"]) + "/raw",
                    "op": "raw_point_decisions",
                    "points": points,
                },
            ]
        )
        raw_requests.append(
            {
                "id": str(frame["plane"]),
                "op": "raw_occupancy",
                "points": points,
                "atoms": frame["aligned_atoms"],
                "sphere_radius": request["config"]["sphere_radius"],
            }
        )
    results = invoke(args.binary, requests, output, "mask")
    old_raw = invoke(args.baseline_binary, raw_requests, output, "baseline_raw")
    frames = []
    for index, frame in enumerate(base["dump"]["orientations"]):
        grid, raw = results[2 * index : 2 * index + 2]
        assert grid["point_source"] == "native_integration_grid"
        assert raw["point_source"] == "supplied_points"
        for key in grid:
            if key not in ("point_source", "id"):
                assert grid[key] == raw[key], key
        assert grid["points"] == points
        assert grid["point_bits"] == [bits(p) for p in points]
        aggregate = grid["aggregate"]
        for field in FIELDS:
            assert (
                bits(aggregate[field])
                == bits(frame[field])
                == bits(old_raw[index][field])
            ), field
            assert bits(aggregate[field]) == aggregate["_bits"][field], field
        assert grid["point_bins"] == old_raw[index]["point_bins"]
        assert grid["grid_count"] == len(points) == base["dump"]["grid_count"]
        assert grid["occupied_count"] == sum(grid["occupied_mask"])
        for label, count in (("quadrant", 4), ("octant", 8)):
            column = 0 if label == "quadrant" else 1
            total = [0] * count
            occupied = [0] * count
            for mask, region in zip(
                grid["occupied_mask"], grid["point_bins"], strict=True
            ):
                total[region[column]] += 1
                occupied[region[column]] += mask
            assert total == grid[label + "_total"]
            assert occupied == grid[label + "_occupied"]
        volume = grid["sphere_volume"]

        def normalized(occupied, total, divisor=1, volume=volume):
            ratio = f32(f32(occupied) / f32(total)) if total else 0.0
            return f32(f32(ratio * volume) / divisor)

        reconstructed = {
            "buried_volume": normalized(grid["occupied_count"], len(points))
        }
        for label, count in (("quadrant", 4), ("octant", 8)):
            reconstructed[label + "s"] = [
                normalized(o, t, count)
                for o, t in zip(
                    grid[label + "_occupied"], grid[label + "_total"], strict=True
                )
            ]
        for name, octants in (
            ("far_vbur", reconstructed["octants"][:4]),
            ("near_vbur", reconstructed["octants"][4:]),
        ):
            value = 0.0
            for octant in octants:
                value = f32(value + octant)
            reconstructed[name] = value
        assert {key: bits(value) for key, value in reconstructed.items()} == aggregate[
            "_bits"
        ]
        frames.append(
            {
                "plane": frame["plane"],
                "grid_count": len(points),
                "occupied_count": grid["occupied_count"],
                "all_aggregate_bits_reconstructed": True,
                "native_masks_identical_for_generated_and_supplied_points": True,
            }
        )
    invalid = invoke(
        args.binary,
        [
            {"id": "empty", "op": "raw_point_decisions", "points": [], "atoms": []},
            {
                "id": "underflow",
                "op": "raw_point_decisions",
                "points": [[0, 0, 0]],
                "atoms": [],
                "config": {"sphere_radius": 1e-20},
            },
            {
                "id": "invalid_density",
                "op": "grid_occupancy",
                "atoms": [],
                "config": {"density": 0},
            },
        ],
        output,
        "invalid_domain",
    )
    assert all(
        "error" in result and "occupied_mask" not in result for result in invalid
    )
    result = {
        "purpose": (
            "Observation adapter integrity; independent science is evaluated separately"
        ),
        "case": case_id,
        "historical_outputs_preserved": True,
        "legacy_geometry_unchanged": True,
        "legacy_raw_aggregate_unchanged": True,
        "frames": frames,
        "invalid_domain_rejected": True,
        "baseline_binary": {
            "path": str(args.baseline_binary.resolve()),
            "sha256": sha(args.baseline_binary),
        },
        "observer_binary": {
            "path": str(args.binary.resolve()),
            "sha256": sha(args.binary),
        },
        "controller_sha256": sha(__file__),
        "inputs": [
            {"path": str(original / name), "sha256": sha(original / name)}
            for name in ("requests.jsonl", "native_stdout.jsonl")
        ],
        "evidence": [
            {"path": str(path), "sha256": sha(path)}
            for path in sorted(output.iterdir())
        ],
    }
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
