"""Independent thermodynamic equations and retained C-audit failure witnesses."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from decimal import Decimal, localcontext
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_data as prep  # noqa: E402
import stericx_quantum as quantum  # noqa: E402

AUDIT_INPUTS = ROOT / "docs/scientific_accuracy_audit/kinetics/inputs"


def populations(energies, temperature):
    """Independent Decimal partition sum using exact SI kB and NA."""
    with localcontext() as context:
        context.prec = 70
        gas = Decimal("1.380649e-23") * Decimal("6.02214076e23") / Decimal(4184)
        values = [Decimal(str(value)) for value in energies]
        low = min(values)
        raw = [
            (-(value - low) / (gas * Decimal(str(temperature)))).exp()
            for value in values
        ]
        return [float(value / sum(raw)) for value in raw]


def backend(temperature):
    result = object.__new__(quantum.QuantumBackend)
    result.config = quantum.QuantumConfig(
        cache_dir=Path("unused_cache"), temperature_k=temperature
    )
    return result


def frames():
    return [
        quantum.XyzFrame(
            ("C",), np.zeros((1, 3)), "analytic", -100 + value / 627.5094740628975
        )
        for value in (0, 1)
    ]


class ThermodynamicsRemediationTests(unittest.TestCase):
    def test_crest_temperature_witness_and_explicit_provenance(self):
        for temperature in (298.15, 500):
            rows = backend(temperature)._conformer_thermodynamics(
                frames(), AUDIT_INPUTS / "crest_default.log"
            )
            np.testing.assert_allclose(
                [row["boltzmann_weight"] for row in rows],
                populations([0, 1], temperature),
                rtol=2e-14,
            )
            self.assertIn("population_method", rows[0])

    def test_negative_populations_rejected_at_any_temperature(self):
        for temperature in (298.15, 500):
            with self.assertRaisesRegex(quantum.QuantumBackendError, "non-negative"):
                backend(temperature)._conformer_thermodynamics(
                    frames(), AUDIT_INPUTS / "crest_negative.log"
                )

    def test_no_table_fallback_declares_unit_degeneracy(self):
        rows = backend(500)._conformer_thermodynamics(
            frames(), AUDIT_INPUTS / "crest_missing.log"
        )
        np.testing.assert_allclose(
            [row["boltzmann_weight"] for row in rows],
            populations([0, 1], 500),
            atol=1e-12,
        )
        self.assertTrue(all(row["degeneracy"] == 1 for row in rows))
        self.assertEqual(
            rows[0]["population_method"], "xyz_electronic_energies_unit_degeneracy"
        )

    def test_complete_rotamer_group_is_sum_not_degeneracy_times_representative(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "crest.log"
            path.write_text(
                "Erel/kcal Etot weight/tot conformer set degen\n"
                "1 0.0 -100 0.4 0.7 1 2\n"
                "2 0.5 -99.999 0.3\n"
                "3 1.0 -99.998 0.3 0.3 2 1\nT /K : 298.15\n"
            )
            rows = backend(500)._conformer_thermodynamics(frames(), path)
            reference = populations([0, 0.5, 1], 500)
            np.testing.assert_allclose(
                [row["boltzmann_weight"] for row in rows],
                [reference[0] + reference[1], reference[2]],
                rtol=2e-14,
            )
            path.write_text(path.read_text().replace("2 0.5 -99.999 0.3\n", ""))
            with self.assertRaisesRegex(
                quantum.QuantumBackendError, "complete rotamer"
            ):
                backend(500)._conformer_thermodynamics(frames(), path)
            # At the documented matching temperature the group table is sufficient.
            self.assertEqual(
                backend(298.15)._conformer_thermodynamics(frames(), path)[0][
                    "boltzmann_weight"
                ],
                0.7,
            )

    def test_invalid_table_does_not_fall_back_to_apparently_valid_xyz(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "crest.log"
            for table in (
                "1 bad 2\n",
                "",
                "1 nan -100 1 1 1 1\n",
                "1 0 -100 1 0 1 1\n",
            ):
                path.write_text(
                    "Erel/kcal Etot weight/tot conformer set degen\n"
                    + table
                    + "T /K : 500\n"
                )
                with self.assertRaises(quantum.QuantumBackendError):
                    backend(500)._conformer_thermodynamics(frames(), path)

    def test_population_normalization_scale_and_domain(self):
        for values in ([1, 3], [1e-300, 3e-300], [1e308, 1e308]):
            expected = [0.5, 0.5] if values[0] == values[1] else [0.25, 0.75]
            np.testing.assert_allclose(
                quantum._normalized_populations(values), expected, rtol=0, atol=0
            )
        for values in ([], [0, 0], [-1, 2], [math.nan, 1], [math.inf, 1]):
            with self.assertRaises(quantum.QuantumBackendError):
                quantum._normalized_populations(values)

    def test_ee_domain_temperature_and_missing_values(self):
        for value in (-101, -100, 100, 101, math.inf, -math.inf):
            with self.assertRaisesRegex(ValueError, "ee measurements"):
                prep.ee_to_ddg(pd.Series([value]), 298.15)
        for temperature in (0, -1, math.nan, math.inf):
            with self.assertRaisesRegex(ValueError, "temperature"):
                prep.ee_to_ddg(pd.Series([50]), temperature)
        self.assertTrue(
            math.isnan(prep.ee_to_ddg(pd.Series([math.nan]), 298.15).iloc[0])
        )

    def test_ee_equation_magnitude_and_tiny_values(self):
        with localcontext() as context:
            context.prec = 70
            gas = Decimal("1.380649e-23") * Decimal("6.02214076e23") / Decimal(4184)
            for value in (0, 1e-12, 4, -88, 99.999):
                x = Decimal(str(abs(value))) / 100
                expected = float(gas * Decimal("353.15") * ((1 + x) / (1 - x)).ln())
                actual = prep.ee_to_ddg(pd.Series([value]), 353.15).iloc[0]
                self.assertAlmostEqual(
                    actual, expected, delta=max(1e-27, abs(expected) * 2e-12)
                )

    def test_ni_hda_metadata_preserves_supplied_target_and_only_derives_missing(self):
        source = pd.DataFrame(
            {
                "id": [2064, 1],
                "smiles": ["P", "P"],
                "ddG_abs": [0.028072, math.nan],
                "ee": [100, 88],
                "nbo_P_boltz": [0, 0],
            }
        )
        result = prep.normalize_public_sigman(source)
        self.assertEqual(result["Temp_K"].tolist(), [353.15, 353.15])
        self.assertEqual(result["Experimental_ddG"].iloc[0], 0.028072)
        self.assertAlmostEqual(
            result["Experimental_ddG"].iloc[1],
            prep.ee_to_ddg(pd.Series([88]), 353.15).iloc[0],
        )

    def test_mmff_retains_status_one_and_validates_temperature(self):
        with (
            patch.object(prep.AllChem, "EmbedMultipleConfs", return_value=(0, 1, 2)),
            patch.object(prep.AllChem, "MMFFHasAllMoleculeParams", return_value=True),
            patch.object(
                prep.AllChem,
                "MMFFOptimizeMoleculeConfs",
                return_value=[(0, 0), (1, 1), (-1, -10)],
            ),
        ):
            result = prep.embed_and_optimize("CP(C)C", 1, 3, 0.1, 6, 500)
            self.assertEqual(result.conformer_ids, (0, 1))
            self.assertEqual(result.mmff_statuses, (0, 1))
            np.testing.assert_allclose(
                result.boltzmann_weights, populations([0, 1], 500), rtol=2e-14
            )
        for temperature in (0, -1, math.nan, math.inf):
            with self.assertRaisesRegex(ValueError, "temperature"):
                prep.embed_and_optimize("P", 1, 1, 0.1, 6, temperature)


if __name__ == "__main__":
    unittest.main()
