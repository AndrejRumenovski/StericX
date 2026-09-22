"""Independent arithmetic of fixed-probe coordinate-rounding clearance intervals.

This is conditional on the represented trial normal. It does not enclose the
uncertainty of the inferred normal or establish experimental geometry errors.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def down(value):
    return math.nextafter(value, -math.inf)


def up(value):
    return math.nextafter(value, math.inf)


def width(value):
    magnitude = np.abs(value)
    adjacent = np.nextafter(magnitude, np.float32(np.inf))
    if np.isinf(adjacent):
        return (float(magnitude) - float(np.nextafter(magnitude, np.float32(0)))) / 2
    return (float(adjacent) - float(magnitude)) / 2


def interval(positions, donor, probe):
    lows, highs = [], []
    for index, atom in enumerate(positions):
        if index == donor:
            continue
        lower = upper = 0.0
        for a, d, p in zip(atom, positions[donor], probe, strict=True):
            da, dd = width(a), width(d)
            a, d, p = float(a), float(d), float(p)
            lo = down(down(down(a - da) - up(d + dd)) - p)
            hi = up(up(up(a + da) - down(d - dd)) - p)
            near = 0.0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
            far = max(abs(lo), abs(hi))
            lower = max(0.0, down(lower + max(0.0, down(near * near))))
            upper = up(upper + up(far * far))
        lows.append(lower)
        highs.append(upper)
    return min(lows), min(highs)


def calculate(request):
    xyz = np.array([a["position"] for a in request["atoms"]], dtype=np.float32)
    donor = request["donor"]
    unit = xyz[request["neighbors"]] - xyz[donor]
    unit /= np.sqrt(np.sum(unit * unit, axis=1))[:, None]
    unit = np.array(sorted(unit.tolist()), dtype=np.float32)
    normal = np.cross(unit[0], unit[1])
    normal /= np.sqrt(np.dot(normal, normal))
    # These witnesses contain no hydrogens; fail instead of silently changing
    # which atoms participate if the diagnostic input scope changes.
    assert all(a["element"].upper() != "H" for a in request["atoms"])
    positive = interval(xyz, donor, normal)
    negative = interval(xyz, donor, -normal)
    selected = (
        1 if positive[0] > negative[1] else -1 if negative[0] > positive[1] else 0
    )
    return {
        "id": request["id"],
        "fixed_represented_normal": normal.tolist(),
        "positive_clearance_squared_interval": positive,
        "negative_clearance_squared_interval": negative,
        "selected_sign": selected,
    }


def main(args):
    rows = [
        calculate(json.loads(line)) for line in args.requests.read_text().splitlines()
    ]
    assert len(rows) == 124
    assert all(row["selected_sign"] == 0 for row in rows)
    control = {
        "id": "clearly_off_plane_obstruction_control",
        "donor": 0,
        "neighbors": [1, 2, 3],
        "atoms": [
            {"element": "P", "position": [0, 0, 0]},
            {"element": "C", "position": [1.8, 0, 0]},
            {"element": "C", "position": [-0.9, 1.5588458, 0]},
            {"element": "C", "position": [-0.9, -1.5588458, 0]},
            {"element": "C", "position": [0, 0, 2.5]},
        ],
    }
    controlled = calculate(control)
    assert controlled["selected_sign"] != 0
    chosen = (
        np.array(controlled["fixed_represented_normal"]) * controlled["selected_sign"]
    )
    assert chosen[2] < 0
    result = {
        "scope": "coordinate rounding only, conditional on fixed trial normal",
        "not_an_uncertain_normal_or_experimental_error_enclosure": True,
        "inputs": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [args.requests, Path(__file__)]
        },
        "rows": rows,
        "control_request": control,
        "control": controlled,
    }
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--requests", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args())
