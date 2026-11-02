"""Independent research references for a union of balls inside a sampling sphere.

This module never imports StericX, DBSTEP or Morfeus. It is not a production
descriptor implementation. Point witnesses use rational arithmetic; continuous
volume bounds use adaptive box inclusion, not either tool's point lattice.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np


def _fraction(value):
    if isinstance(value, Fraction):
        return value
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Reference inputs must be finite")
    return Fraction.from_float(value)


def _float_enclosure(value):
    """Adjacent binary64 bounds of an exactly represented rational number."""
    value = _fraction(value)
    rounded = float(value)
    if not math.isfinite(rounded):
        raise ValueError("Reference quantity is outside finite binary64 range")
    represented = Fraction.from_float(rounded)
    lower = math.nextafter(rounded, -math.inf) if represented > value else rounded
    upper = math.nextafter(rounded, math.inf) if represented < value else rounded
    return lower, upper


def sphere_margin_exact(point, center, radius):
    """Exact sum((point-center)^2)-radius^2 for supplied binary values.

    A nonpositive result means inclusion for a closed ball. In particular this
    predicate does not round a positive distance margin back onto the boundary.
    """
    if len(point) != 3 or len(center) != 3:
        raise ValueError("Point and center must have three coordinates")
    radius = _fraction(radius)
    if radius < 0:
        raise ValueError("Radius must be nonnegative")
    return (
        sum(
            (_fraction(p) - _fraction(c)) ** 2
            for p, c in zip(point, center, strict=True)
        )
        - radius**2
    )


def f32_predicate_certificate(point, center, radius):
    """A priori roundoff certificate, without evaluating a binary32 kernel.

    The documented subtract/square/two-add graph has relative error at most
    gamma_5 on its exact squared distance. A separate radius square contributes
    u*r^2, where u=2^-24. The returned ambiguous region is not an acceptance
    tolerance: every disputed decision must still retain its exact margin and
    an operation-order explanation. Sufficient domain checks exclude underflow,
    overflow and values that are not the actual accepted binary32 inputs.
    """
    if len(point) != 3 or len(center) != 3:
        raise ValueError("Point and center must have three coordinates")
    values = [*map(_fraction, point), *map(_fraction, center), _fraction(radius)]
    largest = Fraction(2**128 - 2**104)
    for value in values:
        if abs(value) > largest or _fraction(np.float32(float(value))) != value:
            raise ValueError("Roundoff certificate requires exact binary32 inputs")
    radius = values[-1]
    if radius <= 0:
        raise ValueError("Roundoff certificate requires a positive radius")
    squares = [(p - c) ** 2 for p, c in zip(values[:3], values[3:6], strict=True)]
    exact_squared = sum(squares)
    radius_squared = radius**2
    unit_roundoff = Fraction(1, 2**24)
    smallest_normal = Fraction(1, 2**126)
    gamma = 5 * unit_roundoff / (1 - 5 * unit_roundoff)
    if (
        any(s and s < smallest_normal / (1 - unit_roundoff) ** 3 for s in squares)
        or exact_squared * (1 + gamma) > largest
        or radius_squared < smallest_normal / (1 - unit_roundoff)
        or radius_squared > largest / (1 + unit_roundoff)
    ):
        raise ValueError("Predicate is outside the certified finite-normal domain")
    margin = exact_squared - radius_squared
    bound = gamma * exact_squared + unit_roundoff * radius_squared
    classification = (
        "stable_inside"
        if margin < -bound
        else "stable_outside"
        if margin > bound
        else "ambiguous"
    )
    return {
        "classification": classification,
        "exact_inside": margin <= 0,
        "margin_exact": str(margin),
        "error_bound_exact": str(bound),
        "error_bound": _float_enclosure(bound)[1],
        "unit_roundoff": "1/16777216",
        "gamma": str(gamma),
        "finite_normal_domain": True,
    }


def _balls(centers, radii):
    centers = np.asarray(centers, dtype=np.float64).reshape((-1, 3))
    radii = np.asarray(radii, dtype=np.float64)
    if radii.shape != (len(centers),):
        raise ValueError("One radius is required per center")
    if not np.isfinite(centers).all() or not np.isfinite(radii).all():
        raise ValueError("Centers and radii must be finite")
    if np.any(radii < 0):
        raise ValueError("Atomic radii must be nonnegative")
    return centers, radii


def union_mask_f64(points, centers, radii):
    """Direct binary64 union membership; this is not an exact predicate.

    Zero-radius entries are explicitly ghost atoms and occupy no nodes. Each
    atom's mask is combined by logical OR, so overlaps are counted only once.
    """
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("Points must be a finite N by 3 array")
    centers, radii = _balls(centers, radii)
    result = np.zeros(len(points), dtype=bool)
    for center, radius in zip(centers, radii, strict=True):
        if radius:
            offsets = points - center
            result |= np.einsum("ij,ij->i", offsets, offsets) <= radius * radius
    return result


def union_mask_exact(points, centers, radii):
    """Rational union predicate intended for small adversarial point sets."""
    centers, radii = _balls(centers, radii)
    return np.asarray(
        [
            any(
                radius > 0 and sphere_margin_exact(point, center, radius) <= 0
                for center, radius in zip(centers, radii, strict=True)
            )
            for point in points
        ],
        dtype=bool,
    )


@lru_cache(maxsize=4)
def pi_interval(terms=64):
    """Rational enclosure from Machin's identity and alternating-series bounds."""
    if not isinstance(terms, int) or terms < 1:
        raise ValueError("A positive integer series length is required")

    def arctangent(reciprocal):
        total = sum(
            Fraction((-1) ** k, (2 * k + 1) * reciprocal ** (2 * k + 1))
            for k in range(terms)
        )
        adjacent = total + Fraction(
            (-1) ** terms, (2 * terms + 1) * reciprocal ** (2 * terms + 1)
        )
        return min(total, adjacent), max(total, adjacent)

    a_lower, a_upper = arctangent(5)
    b_lower, b_upper = arctangent(239)
    return 16 * a_lower - 4 * b_upper, 16 * a_upper - 4 * b_lower


def sphere_intersection_interval(radius, atom_radius, separation):
    """Certified analytic two-ball intersection volume, in input length cubed.

    Inputs are interpreted exactly as their supplied binary values. Separation
    is an explicit scalar distance; uncertainty in deriving it is not included.
    The result brackets the true mathematical pi, not binary64's pi constant.
    """
    radius, atom_radius, separation = map(_fraction, (radius, atom_radius, separation))
    if min(radius, atom_radius, separation) < 0:
        raise ValueError("Radii and separation must be nonnegative")
    if separation >= radius + atom_radius:
        coefficient = Fraction(0)
    elif separation <= abs(radius - atom_radius):
        coefficient = Fraction(4, 3) * min(radius, atom_radius) ** 3
    else:
        coefficient = (
            (radius + atom_radius - separation) ** 2
            * (
                separation**2
                + 2 * separation * (radius + atom_radius)
                - 3 * (radius - atom_radius) ** 2
            )
            / (12 * separation)
        )
    pi_lower, pi_upper = pi_interval()
    lower = _float_enclosure(coefficient * pi_lower)[0]
    upper = _float_enclosure(coefficient * pi_upper)[1]
    return lower, upper


def sphere_intersection_volume(radius, atom_radius, separation):
    lower, upper = sphere_intersection_interval(radius, atom_radius, separation)
    return lower + (upper - lower) / 2


def box_ball_relation_exact(lower, upper, center, radius):
    """Return -1 outside, +1 contained, or 0 unresolved for an exact box."""
    if not all(len(v) == 3 for v in (lower, upper, center)):
        raise ValueError("Box and center require three coordinates")
    low, high, origin = ([*map(_fraction, v)] for v in (lower, upper, center))
    if any(lo > hi for lo, hi in zip(low, high, strict=True)):
        raise ValueError("Box lower boundary exceeds upper boundary")
    radius = _fraction(radius)
    if radius < 0:
        raise ValueError("Radius must be nonnegative")
    r2 = radius**2
    minimum = sum(
        max(lo - c, c - hi, Fraction(0)) ** 2
        for lo, hi, c in zip(low, high, origin, strict=True)
    )
    maximum = sum(
        max(abs(lo - c), abs(hi - c)) ** 2
        for lo, hi, c in zip(low, high, origin, strict=True)
    )
    if minimum > r2:
        return -1
    return 1 if maximum <= r2 else 0


def _box_ball_relation(lower, upper, center, radius):
    """Conservative binary64 interval classifications of enclosed exact boxes."""
    # Bounding boxes enclose their exact rational endpoints. Each arithmetic
    # operation is expanded by one adjacent float in the required direction.
    near = np.maximum(
        0,
        np.maximum(
            np.nextafter(lower - center, -np.inf),
            np.nextafter(center - upper, -np.inf),
        ),
    )
    far = np.maximum(
        abs(np.nextafter(lower - center, -np.inf)),
        abs(np.nextafter(upper - center, np.inf)),
    )
    near2 = np.nextafter(near * near, -np.inf)
    far2 = np.nextafter(far * far, np.inf)
    minimum = np.nextafter(
        np.nextafter(near2[:, 0] + near2[:, 1], -np.inf) + near2[:, 2], -np.inf
    )
    maximum = np.nextafter(
        np.nextafter(far2[:, 0] + far2[:, 1], np.inf) + far2[:, 2], np.inf
    )
    r2_lower = np.nextafter(radius * radius, -np.inf)
    r2_upper = np.nextafter(radius * radius, np.inf)
    return minimum > r2_upper, maximum <= r2_lower


def certified_union_volume(
    centers, radii, radius, center=(0.0, 0.0, 0.0), max_depth=7, batch_size=16384
):
    """Enclose continuous sphere-intersected union volume with dyadic boxes.

    Boxes outside either component are discarded. A box contributes to the
    lower bound only if wholly inside the sampling sphere and at least one
    atomic ball. Unclassified boxes contribute only to the upper bound. This
    sufficient inclusion test is conservative when several balls jointly cover
    a box; overlap is never double counted. Refinement reduces uncertainty.
    """
    centers, radii = _balls(centers, radii)
    center = np.asarray(center, dtype=np.float64)
    radius = float(radius)
    if (
        center.shape != (3,)
        or not np.isfinite(center).all()
        or not math.isfinite(radius)
        or radius <= 0
        or max(abs(center).max(initial=0), abs(centers).max(initial=0), radius) > 1e100
        or radii.max(initial=0) > 1e100
    ):
        raise ValueError(
            "Finite molecular geometry and positive sphere radius required"
        )
    if not isinstance(max_depth, int) or not 0 <= max_depth <= 20:
        raise ValueError("max_depth must be an integer from 0 through 20")
    if not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    root_low = np.asarray(
        [_float_enclosure(_fraction(c) - _fraction(radius))[0] for c in center]
    )
    root_high = np.asarray(
        [_float_enclosure(_fraction(c) + _fraction(radius))[1] for c in center]
    )
    retained = np.asarray(
        [
            atom_radius > 0
            and not _box_ball_relation(
                root_low[None, :], root_high[None, :], atom_center, atom_radius
            )[0][0]
            for atom_center, atom_radius in zip(centers, radii, strict=True)
        ],
        dtype=bool,
    )
    centers, radii = centers[retained], radii[retained]
    root_volume = (2 * _fraction(radius)) ** 3
    final_units = 8**max_depth
    inside_units = 0
    active = np.zeros((1, 3), dtype=np.int64)
    child_offsets = np.asarray(
        list(itertools.product((0, 1), repeat=3)), dtype=np.int64
    )
    history = []
    for depth in range(max_depth + 1):
        axis_enclosures = []
        divisions = 2**depth
        for value in center:
            axis_enclosures.append(
                np.asarray(
                    [
                        _float_enclosure(
                            _fraction(value)
                            + _fraction(radius) * Fraction(2 * k - divisions, divisions)
                        )
                        for k in range(divisions + 1)
                    ]
                )
            )
        uncertain_batches = []
        inside_count = outside_count = 0
        for start in range(0, len(active), batch_size):
            indices = active[start : start + batch_size]
            low = np.column_stack(
                [axis[indices[:, j], 0] for j, axis in enumerate(axis_enclosures)]
            )
            high = np.column_stack(
                [axis[indices[:, j] + 1, 1] for j, axis in enumerate(axis_enclosures)]
            )
            sphere_out, sphere_in = _box_ball_relation(low, high, center, radius)
            union_out = np.ones(len(indices), dtype=bool)
            union_in = np.zeros(len(indices), dtype=bool)
            for atom_center, atom_radius in zip(centers, radii, strict=True):
                atom_out, atom_in = _box_ball_relation(
                    low, high, atom_center, atom_radius
                )
                union_out &= atom_out
                union_in |= atom_in
            definitely_out = sphere_out | union_out
            definitely_in = sphere_in & union_in
            if np.any(definitely_out & definitely_in):
                raise ArithmeticError("Contradictory interval box classification")
            inside_count += int(np.count_nonzero(definitely_in))
            outside_count += int(np.count_nonzero(definitely_out))
            uncertain_batches.append(indices[~(definitely_out | definitely_in)])
        unknown = (
            np.concatenate(uncertain_batches)
            if uncertain_batches
            else np.empty((0, 3), dtype=np.int64)
        )
        unit_multiplier = 8 ** (max_depth - depth)
        inside_units += inside_count * unit_multiplier
        upper_units = inside_units + len(unknown) * unit_multiplier
        lower_exact = root_volume * Fraction(inside_units, final_units)
        upper_exact = root_volume * Fraction(upper_units, final_units)
        lower = _float_enclosure(lower_exact)[0]
        upper = _float_enclosure(upper_exact)[1]
        history.append(
            {
                "depth": depth,
                "boxes_examined": len(active),
                "new_inside_boxes": inside_count,
                "new_outside_boxes": outside_count,
                "unresolved_boxes": len(unknown),
                "lower": lower,
                "upper": upper,
                "width": (
                    math.nextafter(upper - lower, math.inf) if upper > lower else 0.0
                ),
            }
        )
        if depth == max_depth or not len(unknown):
            break
        active = (2 * unknown[:, None, :] + child_offsets[None, :, :]).reshape((-1, 3))
    return {
        "lower": lower,
        "upper": upper,
        "width": history[-1]["width"],
        "lower_exact": str(lower_exact),
        "upper_exact": str(upper_exact),
        "history": history,
        "maximum_depth": max_depth,
        "scope": (
            "Continuous union of supplied exact binary-coordinate balls; "
            "no source-geometry or physical-radius uncertainty."
        ),
    }


def ratio_volume_error_bound(
    occupied_estimate, sphere_estimate, occupied_error, sphere_error
):
    """Bound analytic-sphere-normalized error from two grid-volume bounds.

    If |Qocc-Vocc| <= Eocc and |Qsphere-Vsphere| <= Esphere, then
    |Vsphere*Qocc/Qsphere - Vocc| <= Eocc + (Qocc/Qsphere)*Esphere.
    Numerator and denominator must each have independently justified bounds.
    No clamp to [0, 100%] is performed or validated here.
    """
    q_occ, q_sphere, e_occ, e_sphere = map(
        _fraction, (occupied_estimate, sphere_estimate, occupied_error, sphere_error)
    )
    if q_sphere <= 0 or min(q_occ, e_occ, e_sphere) < 0:
        raise ValueError(
            "Positive denominator and nonnegative estimates/errors required"
        )
    return _float_enclosure(e_occ + q_occ / q_sphere * e_sphere)[1]
