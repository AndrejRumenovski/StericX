# ruff: noqa: E501
"""Independent geometry mathematics; no imports/calls to StericX kernels.

Sterimol B1 uses exact analytic candidates of the support-envelope minimum,
not StericX's angular scan. BV uses vectorized union-of-balls membership on
independently constructed float64 cubature and explicit per-octant masks.
"""

import itertools
import math
import warnings

import numpy as np
from morfeus import BuriedVolume, Pyramidalization, Sterimol
from morfeus.utils import get_radii
from scipy.spatial import cKDTree


def radii_for(case):
    # Reference Bondi implementation independent of StericX's radius table.
    vals = get_radii([a["element"] for a in case["atoms"]], radii_type="bondi")
    return np.array(
        [a.get("radius", r) for a, r in zip(case["atoms"], vals, strict=True)]
    )


def unit(v):
    n = np.linalg.norm(v)
    if n == 0 or not np.isfinite(n):
        raise ValueError("undefined direction")
    return v / n


def center_for(case, xyz):
    if "center" in case:
        return np.array(case["center"])
    donor = case["donor"]
    vec = xyz[case["neighbors"]] - xyz[donor]
    vec = vec / np.linalg.norm(vec, axis=1)[:, None]
    s = -vec.sum(axis=0)
    if np.linalg.norm(s) <= 1e-12:
        raise ValueError("planar donor has no unique opposing-bond direction")
    return xyz[donor] + case["config"]["center_distance"] * unit(s)


def sterimol_exact(xyz, radii, dummy, base, excluded=None):
    direction = unit(xyz[base] - dummy)
    axis_candidates = np.eye(3)
    seed = axis_candidates[np.argmin(abs(axis_candidates @ direction))]
    x = unit(np.cross(direction, seed))
    y = np.cross(direction, x)
    keep = np.ones(len(xyz), dtype=bool)
    if excluded is not None:
        keep[excluded] = False
    relative = xyz[keep] - dummy
    rs = radii[keep]
    axial = relative @ direction
    transverse = relative @ np.array([x, y]).T
    angles = [0.0]
    angles.extend((np.arctan2(transverse[:, 1], transverse[:, 0]) + np.pi).tolist())
    for i, j in itertools.combinations(range(len(rs)), 2):
        dx, dy = transverse[i] - transverse[j]
        amplitude = np.hypot(dx, dy)
        if amplitude == 0:
            continue
        ratio = (rs[j] - rs[i]) / amplitude
        if abs(ratio) <= 1:
            phase = np.arctan2(dy, dx)
            delta = np.arccos(ratio)
            angles.extend([phase + delta, phase - delta])
    directions = np.column_stack([np.cos(angles), np.sin(angles)])
    support = (transverse @ directions.T + rs[:, None]).max(axis=0)
    return {
        "l": float(max(axial + rs)),
        "b1": float(min(support)),
        "b5": float(max(np.linalg.norm(transverse, axis=1) + rs)),
    }


def sterimol_refs(case, xyz, controlled_center=None):
    radii = radii_for(case)
    d = case["donor"]
    attach = case.get("attach", d)
    base = case.get("sterimol_neighbor", case["reference"])
    exact = sterimol_exact(xyz, radii, xyz[attach], base, attach)
    m = Sterimol(
        [a["element"] for a in case["atoms"]],
        xyz.copy(),
        attach + 1,
        base + 1,
        radii=radii,
    )
    out = {
        "bond": {
            "independent": exact,
            "morfeus": {
                "l": float(m.L_value_uncorrected),
                "b1": float(m.B_1_value),
                "b5": float(m.B_5_value),
            },
            "morfeus_corrected_l": float(m.L_value),
        }
    }
    if not case.get("sterimol_only"):
        try:
            center = (
                center_for(case, xyz)
                if controlled_center is None
                else np.asarray(controlled_center)
            )
            exact = sterimol_exact(xyz, radii, center, d)
            mm = Sterimol(
                ["H"] + [a["element"] for a in case["atoms"]],
                np.vstack([center, xyz]),
                1,
                d + 2,
                radii=np.r_[0.0, radii],
            )
            out["coordination"] = {
                "center": center.tolist(),
                "independent": exact,
                "morfeus": {
                    "l": float(mm.L_value_uncorrected),
                    "b1": float(mm.B_1_value),
                    "b5": float(mm.B_5_value),
                },
            }
        except (ValueError, FloatingPointError) as e:
            out["coordination"] = {"error": str(e)}
    return out


def pyramidalization_ref(case, xyz):
    v = xyz[case["neighbors"]] - xyz[case["donor"]]
    if v.shape != (3, 3):
        raise ValueError("requires three neighbors")
    norms = np.linalg.norm(v, axis=1)
    if np.any(norms == 0):
        raise ValueError("coincident donor and substituent")
    v = v / norms[:, None]
    determinant = abs(np.linalg.det(v))
    angles = []
    for k in range(3):
        other = np.delete(v, k, axis=0)
        normal = np.cross(*other)
        length = np.linalg.norm(normal)
        if length == 0:
            raise ValueError("collinear substituents")
        cosine = np.clip(abs(v[k] @ normal) / length, 0, 1)
        angle = math.acos(cosine)
        if np.dot(other.sum(axis=0), v[k]) > 0:
            angle = -angle
        angles.append(angle)
    mean = np.mean(angles)
    independent = {
        "pyr_p": float(2 - determinant if mean < 0 else determinant),
        "pyr_alpha": float(np.rad2deg(mean)),
    }
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = Pyramidalization(
            xyz.copy(),
            case["donor"] + 1,
            neighbor_indices=[i + 1 for i in case["neighbors"]],
        )
    return {
        "independent": independent,
        "morfeus": {"pyr_p": float(m.P), "pyr_alpha": float(m.alpha)},
        "morfeus_warnings": [str(w.message) for w in caught],
    }


def grid(radius, density, dtype=np.float64):
    n = round((8 * radius**3 / density) ** (1 / 3))
    axes = np.linspace(-radius, radius, n, dtype=dtype)
    p = np.stack(np.meshgrid(axes, axes, axes, indexing="ij"), axis=-1).reshape(-1, 3)
    return p[np.sum(p * p, axis=1) <= radius * radius]


def volume_on_grid(points, centers, radii, R):
    occupied = np.zeros(len(points), bool)
    tree = cKDTree(points)
    for c, r in zip(centers, radii, strict=True):
        occupied[tree.query_ball_point(c, r)] = True
    sphere = 4 * np.pi * R**3 / 3
    # Explicit coordinate sign masks; canonical IDs (++),( -+),(--),(+-), z+ then z-.
    octants = []
    counts = []
    occupied_counts = []
    for z in [1, -1]:
        for x, y in [(1, 1), (-1, 1), (-1, -1), (1, -1)]:
            mask = (
                (points[:, 0] * x > 0) & (points[:, 1] * y > 0) & (points[:, 2] * z > 0)
            )
            n = mask.sum()
            k = (occupied & mask).sum()
            counts.append(int(n))
            occupied_counts.append(int(k))
            octants.append(float(k / n * sphere / 8) if n else None)
    quadrants = (
        [octants[i] + octants[i + 4] for i in range(4)]
        if all(x is not None for x in octants)
        else [None] * 4
    )
    v = float(occupied.mean() * sphere)
    return {
        "buried_volume": v,
        "percent_buried_volume": float(occupied.mean() * 100),
        "quadrants": quadrants,
        "octants": octants,
        "near_vbur": sum(octants[4:]),
        "far_vbur": sum(octants[:4]),
        "grid_count": len(points),
        "octant_grid_counts": counts,
        "octant_occupied_counts": occupied_counts,
        "occupied_total": int(occupied.sum()),
    }


def bv_refs(case, xyz, center=None):
    cfg = case["config"]
    R = cfg["sphere_radius"]
    density = cfg["density"]
    center = center_for(case, xyz) if center is None else np.asarray(center)
    radii = radii_for(case) * cfg["radii_scale"]
    elements = [a["element"] for a in case["atoms"]]
    d = case["donor"]
    orientations = []
    points = grid(R, density)
    for plane in case["neighbors"]:
        z = unit(center - xyz[d])
        x = unit((xyz[plane] - center) - np.dot(xyz[plane] - center, z) * z)
        y = np.cross(z, x)
        aligned = (xyz - center) @ np.array([x, y, z]).T
        include = np.array([cfg["include_hydrogens"] or e != "H" for e in elements])
        v = volume_on_grid(points, aligned[include], radii[include], R)
        # Present zero metal at origin explicitly to avoid center aliasing differences in reference tool.
        m = BuriedVolume(
            ["H", *elements],
            np.vstack([[0.0, 0.0, 0.0], aligned]),
            1,
            radii=np.r_[0.0, radii],
            radius=R,
            density=density,
            include_hs=cfg["include_hydrogens"],
        )
        m.octant_analysis()
        order = [0, 1, 2, 3, 7, 6, 5, 4]
        mo = [m.octants["buried_volume"][i] for i in order]
        orientations.append(
            {
                "plane": plane,
                "independent": v,
                "morfeus": {
                    "buried_volume": float(m.buried_volume),
                    "percent_buried_volume": float(m.fraction_buried_volume * 100),
                    "quadrants": [float(mo[i] + mo[i + 4]) for i in range(4)],
                    "octants": [float(t) for t in mo],
                    "near_vbur": float(sum(mo[4:])),
                    "far_vbur": float(sum(mo[:4])),
                },
            }
        )

    def aggregate(key):
        rows = [o[key] for o in orientations]
        q = np.array([r["quadrants"] for r in rows])
        octs = np.array([r["octants"] for r in rows])
        first = rows[0]
        return {
            k: first[k]
            for k in ["buried_volume", "percent_buried_volume", "near_vbur", "far_vbur"]
        } | {
            "qvbur_min": float(q.min()),
            "qvbur_max": float(q.max()),
            "ovbur_min": float(octs.min()),
            "ovbur_max": float(octs.max()),
            "max_delta_qvbur": float(abs(q - np.roll(q, 1, axis=1)).max()),
        }

    return {
        "center": center.tolist(),
        "orientations": orientations,
        "independent": aggregate("independent"),
        "morfeus": aggregate("morfeus"),
    }
