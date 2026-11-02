"""Independent-reference integration regressions for finite-domain comparisons.

Run explicitly in the frozen DBSTEP environment. General research environments
without DBSTEP or the source-bound observer report a skip, never a reference pass.
"""

import importlib.util
import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
EXECUTABLE = Path(
    os.environ.get(
        "STERICX_RESOLUTION_OBSERVER",
        ROOT / ".stericx/dbstep_resolution/mask_build_v2/bin/stericx-audit-observer",
    )
)
AVAILABLE = importlib.util.find_spec("dbstep") is not None and EXECUTABLE.is_file()


@unittest.skipUnless(AVAILABLE, "requires frozen DBSTEP and source-bound Rust observer")
class PointDomainRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from analyze_dbstep_point_domains import analyze_mask, verify_counts
        from dbstep_resolution_reference import (
            f32_predicate_certificate,
            sphere_margin_exact,
        )

        cls.analyze_mask = staticmethod(analyze_mask)
        cls.verify_counts = staticmethod(verify_counts)
        cls.certificate = staticmethod(f32_predicate_certificate)
        cls.margin = staticmethod(sphere_margin_exact)

    def observe(self, points=None, atoms=None, radius=1):
        payload = {
            "id": "integration-regression",
            "op": "raw_point_decisions" if points is not None else "grid_occupancy",
            "atoms": atoms or [],
            "config": {"sphere_radius": radius, "density": 0.01},
            "include_points": True,
        }
        if points is not None:
            payload["points"] = points
        process = subprocess.run(
            [str(EXECUTABLE)],
            input=json.dumps(payload) + "\n",
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(process.stderr, "")
        result = json.loads(process.stdout)
        self.assertNotIn("error", result)
        return result

    def test_original_failure_is_a_second_clip_not_atomic_occupation(self):
        frozen = ROOT / "docs/validation/dbstep/resolution_baseline/original_failure"
        response = json.loads((frozen / "native_stdout.jsonl").read_text())
        request = json.loads((frozen / "requests.jsonl").read_text())
        points = np.float32(response["dump"]["points"]).astype(float)
        radii = np.float32(
            [
                np.float32(a["vdw_radius"])
                * np.float32(request["config"]["radii_scale"])
                for a in response["atoms"]
                if a["element"].upper() != "H"
            ]
        ).astype(float)
        for frame in response["dump"]["orientations"]:
            observed = self.observe(points.tolist(), frame["aligned_atoms"], radius=2.5)
            centers = np.float32(
                [a["position"] for a in frame["aligned_atoms"]]
            ).astype(float)
            analysis, tree, mask = self.analyze_mask(points, centers, radii, observed)
            self.assertEqual(analysis["native_dbstep_mask_differences"], 0)
            self.assertEqual(analysis["native_direct_mask_differences"], 0)
            self.assertEqual(self.verify_counts(observed), [])
            retained = tree.query_ball_point(np.zeros(3), 2.5)
            self.assertEqual((len(points), len(retained)), (5565, 5533))
            removed = sorted(set(range(len(points))) - set(retained))
            self.assertEqual(int(mask[removed].sum()), 6)
            for index in removed:
                self.assertGreater(self.margin(points[index], [0, 0, 0], 2.5), 0)
                self.assertEqual(
                    self.certificate(points[index], [0, 0, 0], 2.5)["classification"],
                    "ambiguous",
                )

    def test_exact_and_adjacent_float_atom_boundaries(self):
        points = np.float32(
            [
                [1, 0, 0],
                [np.nextafter(np.float32(1), np.float32(2)), 0, 0],
                [0.6, 0.8, 0],
            ]
        )
        observed = self.observe(
            points.tolist(), [{"position": [0, 0, 0], "radius_squared": 1}]
        )
        self.assertEqual(observed["occupied_mask"], [True, False, True])
        analysis, _, _ = self.analyze_mask(
            points.astype(float), np.zeros((1, 3)), np.ones(1), observed
        )
        self.assertEqual(analysis["native_dbstep_mask_differences"], 1)
        self.assertEqual(analysis["unexplained_count"], 0)
        self.assertFalse(analysis["disputed_points"][0]["exact_rational"])

    def test_checker_rejects_a_stable_wrong_kernel_decision(self):
        observed = {"occupied_mask": [True]}
        analysis, _, _ = self.analyze_mask(
            np.array([[2.0, 0, 0]]), np.zeros((1, 3)), np.ones(1), observed
        )
        self.assertEqual(analysis["unexplained_count"], 1)
        self.assertTrue(analysis["disputed_points"][0]["stable_disagreement"])

    def test_region_zero_planes_have_declared_positive_assignment(self):
        observed = self.observe(
            [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
            [{"position": [0, 0, 0], "radius_squared": 4}],
        )
        self.assertEqual(
            observed["point_bins"], [[0, 0], [1, 1], [0, 0], [3, 3], [0, 0], [0, 4]]
        )
        self.assertEqual(self.verify_counts(observed), [])

    def test_row_cache_handles_point_and_atom_permutations(self):
        atoms = [
            {"position": [-0.5, 0, 0], "radius_squared": 0.5625},
            {"position": [0.5, 0.25, 0], "radius_squared": 0.25},
        ]
        original = self.observe(atoms=atoms)
        points = np.float32(original["points"])
        order = np.random.default_rng(20261102).permutation(len(points))
        changed = self.observe(points[order].tolist(), atoms[::-1])
        np.testing.assert_array_equal(
            np.array(original["occupied_mask"])[order], changed["occupied_mask"]
        )
        self.assertEqual(original["aggregate"]["_bits"], changed["aggregate"]["_bits"])
        self.assertEqual(self.verify_counts(changed), [])

    def test_empty_union_is_zero_without_removing_nodes(self):
        observed = self.observe()
        self.assertGreater(observed["grid_count"], 0)
        self.assertEqual(observed["occupied_count"], 0)
        self.assertEqual(observed["aggregate"]["buried_volume"], 0)


if __name__ == "__main__":
    unittest.main()
