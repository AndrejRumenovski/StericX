"""Independent analytical checks of the diagnostic cell-bound derivation."""

import importlib.util
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

REFERENCE_AVAILABLE = all(
    importlib.util.find_spec(name) is not None for name in ("dbstep", "morfeus")
)

if REFERENCE_AVAILABLE:
    from analyze_dbstep_volume_convergence import (
        analytical_interval,
        cell_geometry,
        estimator_enclosure,
        native_ball_trace,
        primitive_cases,
        shell_upper,
        tensor_points,
    )
    from dbstep_resolution_reference import union_mask_exact


@unittest.skipUnless(
    REFERENCE_AVAILABLE,
    "requires pinned DBSTEP/Morfeus; seven tests mandatory in reference environment",
)
class ConvergenceCellBounds(unittest.TestCase):
    def test_equal_sphere_nonuniform_axis_with_exact_point_reference(self):
        axis = np.array([-1.25, -0.8, -0.3, 0.0, 0.2, 0.7, 1.25])
        axes = [axis] * 3
        points = tensor_points(axes)
        mask = union_mask_exact(points, [[0, 0, 0]], [1.0])
        cells = cell_geometry(axes, 1.0)
        assert cells["weight_error_per_selected_node"] > 0
        observed = 4 * math.pi / 3
        bound, _ = estimator_enclosure(
            axes, [[0, 0, 0]], [1.0], 1.0, mask, mask, observed
        )
        assert (
            bound["certified_volume_lower"]
            <= observed
            <= bound["certified_volume_upper"]
        )

    def test_arbitrary_mask_defect_is_charged_and_truth_still_enclosed(self):
        axis = np.linspace(-1.2, 1.2, 13)
        axes = [axis] * 3
        points = tensor_points(axes)
        sphere = union_mask_exact(points, [[0, 0, 0]], [1.0])
        occupied = union_mask_exact(points, [[0, 0, 0]], [0.5])
        occupied[0] = True  # Deliberate exterior decision defect.
        value = (4 * math.pi / 3) * occupied.sum() / sphere.sum()
        bound, differing = estimator_enclosure(
            axes, [[0, 0, 0]], [0.5], 1.0, sphere, occupied, value
        )
        assert 0 in differing["numerator_disagreement_indices"]
        truth = math.pi / 6
        assert (
            bound["certified_volume_lower"] <= truth <= bound["certified_volume_upper"]
        )

    def test_finite_endpoint_coverage_is_mandatory(self):
        with self.assertRaisesRegex(ValueError, "do not cover"):
            cell_geometry([np.array([-0.5, 0, 0.5])] * 3, 1.0)

    def test_shell_contains_exact_inner_outer_ball_difference(self):
        radius, delta = Fraction(3, 2), Fraction(1, 8)
        expected = (4 * math.pi / 3) * float(
            (radius + delta) ** 3 - (radius - delta) ** 3
        )
        assert float(shell_upper(radius, delta)) >= expected

    def test_analytical_primitive_limits_and_overlap_are_independent(self):
        cases = {case["id"]: case for case in primitive_cases()}
        truth = {name: analytical_interval(case) for name, case in cases.items()}
        assert truth["empty"] == [0, 0]
        assert truth["disjoint"] == [0, 0]
        assert truth["external_tangent"] == [0, 0]
        assert truth["coincident"][0] <= 4 * math.pi / 3 <= truth["coincident"][1]
        assert truth["separated"][0] <= math.pi / 3 <= truth["separated"][1]
        # Two unit balls separated by one: lens volume is 5*pi/12.
        union = 9 * math.pi / 4
        assert truth["overlapping"][0] <= union <= truth["overlapping"][1]

    def test_empty_union_normalization_is_exact_zero(self):
        axis = np.linspace(-1.25, 1.25, 11)
        axes = [axis] * 3
        points = tensor_points(axes)
        sphere = union_mask_exact(points, [[0, 0, 0]], [1])
        empty = np.zeros(len(points), dtype=bool)
        bound, _ = estimator_enclosure(axes, [], [], 1, sphere, empty, 0)
        assert bound["total_absolute_error_upper"] == 0
        assert bound["certified_volume_lower"] == 0
        assert bound["certified_volume_upper"] == 0

    def test_native_circle_boundary_trace_has_rational_explanation(self):
        point = np.asarray([0.6, 0.8, 0], dtype=np.float32).astype(float)
        witness = native_ball_trace(point, [0, 0, 0], 1)
        assert witness["trace_inside"]
        assert not witness["certificate"]["exact_inside"]
        assert witness["certificate"]["classification"] == "ambiguous"
        assert witness["trace_exact_dispute_within_bound"]


if __name__ == "__main__":
    unittest.main()
