"""Independent analytic and interval checks of the research volume reference."""

import math
import sys
import unittest
from fractions import Fraction
from itertools import pairwise
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import dbstep_resolution_reference as reference


class RationalGeometry(unittest.TestCase):
    def test_roundoff_certificate_uses_an_independent_error_bound(self):
        inside = reference.f32_predicate_certificate([0.5, 0, 0], [0, 0, 0], 1)
        outside = reference.f32_predicate_certificate([2, 0, 0], [0, 0, 0], 1)
        boundary = reference.f32_predicate_certificate([1, 0, 0], [0, 0, 0], 1)
        self.assertEqual(inside["classification"], "stable_inside")
        self.assertEqual(outside["classification"], "stable_outside")
        self.assertEqual(boundary["classification"], "ambiguous")
        point = [-2.045454502105713, -0.4545454978942871, 1.3636364936828613]
        witness = reference.f32_predicate_certificate(point, [0, 0, 0], 2.5)
        self.assertEqual(witness["classification"], "ambiguous")
        self.assertFalse(witness["exact_inside"])
        self.assertLess(
            abs(Fraction(witness["margin_exact"])),
            Fraction(witness["error_bound_exact"]),
        )

    def test_roundoff_bound_refuses_unjustified_domains(self):
        for point, radius in [
            ([0.6, 0, 0], 1),  # binary64 0.6 is not the accepted binary32 value
            ([float(np.float32(2.0**-70)), 0, 0], 1),  # squared underflow
            ([float(np.float32(2.0**70)), 0, 0], 1),  # squared overflow
            ([0, 0, 0], float(np.float32(2.0**-70))),
        ]:
            with (
                self.subTest(point=point, radius=radius),
                self.assertRaises(ValueError),
            ):
                reference.f32_predicate_certificate(point, [0, 0, 0], radius)

    def test_binary_boundary_is_not_decimal_boundary(self):
        self.assertGreater(
            reference.sphere_margin_exact([0.6, 0.8, 0], [0, 0, 0], 1), 0
        )
        self.assertEqual(reference.sphere_margin_exact([1, 0, 0], [0, 0, 0], 1), 0)
        point = [-2.045454502105713, -0.4545454978942871, 1.3636364936828613]
        margin = reference.sphere_margin_exact(point, [0, 0, 0], 2.5)
        self.assertEqual(
            margin, Fraction("0.000000216744183489936403930187225341796875")
        )

    def test_exact_union_overlaps_ghosts_and_boundary(self):
        points = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [4, 0, 0]]
        centers, radii = [[0, 0, 0], [1, 0, 0], [4, 0, 0]], [1, 1, 0]
        expected = [True, True, True, False]
        self.assertEqual(
            reference.union_mask_exact(points, centers, radii).tolist(), expected
        )
        self.assertEqual(
            reference.union_mask_f64(points, centers, radii).tolist(), expected
        )

    def test_machin_pi_enclosure(self):
        lower, upper = reference.pi_interval()
        self.assertGreater(
            lower, Fraction("3.14159265358979323846264338327950288419716939937510")
        )
        self.assertLess(
            upper, Fraction("3.14159265358979323846264338327950288419716939937511")
        )
        self.assertLess(upper - lower, Fraction(1, 10**80))

    def test_analytic_lenses_cover_all_branches(self):
        pi_lower, pi_upper = reference.pi_interval()
        for r, s, d, coefficient in [
            (1, 1, 0, Fraction(4, 3)),
            (2, 1, 0.5, Fraction(4, 3)),
            (1, 1, 1, Fraction(5, 12)),
            (1, 1, 2, Fraction(0)),
            (1, 1, 3, Fraction(0)),
            (1, 0, 0, Fraction(0)),
        ]:
            with self.subTest(r=r, s=s, d=d):
                lower, upper = reference.sphere_intersection_interval(r, s, d)
                self.assertLessEqual(Fraction(lower), coefficient * pi_lower)
                self.assertGreaterEqual(Fraction(upper), coefficient * pi_upper)
                self.assertEqual(
                    (lower, upper), reference.sphere_intersection_interval(s, r, d)
                )
        near_tangent = reference.sphere_intersection_interval(
            1, 1, math.nextafter(2, 0)
        )
        self.assertGreater(near_tangent[0], 0)

    def test_direct_membership_is_order_independent(self):
        rng = np.random.default_rng(20261102)
        points = rng.uniform(-2, 2, (100, 3))
        centers = rng.uniform(-1, 1, (8, 3))
        radii = rng.uniform(0.2, 1.2, 8)
        exact = reference.union_mask_exact(points, centers, radii)
        np.testing.assert_array_equal(
            exact, reference.union_mask_f64(points, centers, radii)
        )
        permutation = rng.permutation(8)
        np.testing.assert_array_equal(
            exact,
            reference.union_mask_f64(points, centers[permutation], radii[permutation]),
        )


class CertifiedBoxes(unittest.TestCase):
    def test_outward_classifications_against_exact_rational_boxes(self):
        rng = np.random.default_rng(4729)
        lows = rng.uniform(-3, 3, (2048, 3))
        highs = lows + rng.uniform(0, 2, (2048, 3))
        origin = np.array([0.6, 0.8, -0.1])
        outside, inside = reference._box_ball_relation(lows, highs, origin, 1.3)
        self.assertFalse(np.any(outside & inside))
        self.assertGreater(int(outside.sum()), 0)
        self.assertGreater(int(inside.sum()), 0)
        for i, (low, high) in enumerate(zip(lows, highs, strict=True)):
            exact = reference.box_ball_relation_exact(low, high, origin, 1.3)
            if outside[i]:
                self.assertEqual(exact, -1)
            if inside[i]:
                self.assertEqual(exact, 1)

    def test_tangent_and_adjacent_float_boxes_do_not_get_false_certificates(self):
        for x in (math.nextafter(1, 0), 1.0, math.nextafter(1, math.inf)):
            low = np.array([[x, 0, 0]])
            high = low.copy()
            outside, inside = reference._box_ball_relation(low, high, np.zeros(3), 1)
            exact = reference.box_ball_relation_exact(low[0], high[0], [0, 0, 0], 1)
            if outside[0]:
                self.assertEqual(exact, -1)
            if inside[0]:
                self.assertEqual(exact, 1)

    def assert_contains_analytic(self, result, analytical):
        self.assertLessEqual(result["lower"], analytical[0])
        self.assertGreaterEqual(result["upper"], analytical[1])
        history = result["history"]
        for old, new in pairwise(history):
            self.assertLessEqual(old["lower"], new["lower"])
            self.assertGreaterEqual(old["upper"], new["upper"])

    def test_single_atom_bounds_shrink_and_contain_analytic_volume(self):
        analytical = reference.sphere_intersection_interval(1, 0.75, 0.5)
        result = reference.certified_union_volume([[0.5, 0, 0]], [0.75], 1, max_depth=7)
        self.assert_contains_analytic(result, analytical)
        self.assertLess(result["upper"] - result["lower"], 0.25)
        self.assertGreater(len(result["history"]), 3)

    def test_overlapping_spheres_never_double_count(self):
        result = reference.certified_union_volume(
            [[-0.5, 0, 0], [0.5, 0, 0]], [1, 1], 2, max_depth=6
        )
        pi_lower, pi_upper = reference.pi_interval()
        analytical = (
            reference._float_enclosure(Fraction(9, 4) * pi_lower)[0],
            reference._float_enclosure(Fraction(9, 4) * pi_upper)[1],
        )
        self.assert_contains_analytic(result, analytical)
        self.assertLess(result["upper"], 8 * math.pi / 3)

    def test_empty_and_completely_covering_union(self):
        empty = reference.certified_union_volume([], [], 1)
        self.assertEqual((empty["lower"], empty["upper"]), (0, 0))
        full = reference.certified_union_volume([[0, 0, 0]], [3], 1, max_depth=5)
        self.assert_contains_analytic(
            full, reference.sphere_intersection_interval(1, 3, 0)
        )

    def test_exact_translation_permutation_and_batch_size(self):
        centers = np.array([[-0.5, 0, 0], [0.5, 0, 0]])
        base = reference.certified_union_volume(centers, [1, 1], 2, max_depth=4)
        translated = reference.certified_union_volume(
            centers[::-1] + np.array([8, -4, 2]),
            [1, 1],
            2,
            center=[8, -4, 2],
            max_depth=4,
            batch_size=11,
        )
        self.assertEqual(base, translated)

    def test_ratio_error_bound_includes_denominator(self):
        # Vocc=1, Vs=4, Qocc=1.1, Qs=3.8; numerator alone underbounds error.
        bound = reference.ratio_volume_error_bound(
            1.1, 3.8, 0.100000000000001, 0.200000000000001
        )
        self.assertGreater(bound, 0.1)
        self.assertGreaterEqual(bound, abs(4 * 1.1 / 3.8 - 1))
        self.assertEqual(reference.ratio_volume_error_bound(1, 2, 0, 0), 0)

    def test_invalid_inputs_are_errors(self):
        with self.assertRaises(ValueError):
            reference.sphere_intersection_interval(1, -1, 0)
        with self.assertRaises(ValueError):
            reference.certified_union_volume([[0, 0, 0]], [1], 0)
        with self.assertRaises(ValueError):
            reference.ratio_volume_error_bound(1, 0, 0, 0)
        with self.assertRaises(ValueError):
            reference.box_ball_relation_exact([0, 0, 0], [1, 1, 1], [0, 0, 0], -1)


if __name__ == "__main__":
    unittest.main()
