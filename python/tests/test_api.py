import json
import math
import struct
from pathlib import Path

import pytest
import stericx as sx

XYZ = "4\nligand\nP 0 0 0\nC 1.4 0 .45\nC -.7 1.212 .45\nC -.7 -1.212 .45\n"
SDF = (
    "ligand\nStericX\n\n  4  3  0  0  0  0            999 V2000\n"
    "0 0 0 P\n1.4 0 .45 C\n-.7 1.212 .45 C\n-.7 -1.212 .45 C\n"
    "  1  2  1  0\n  1  3  1  0\n  1  4  1  0\nM  END\n$$$$\n"
)
COORD = sx.DescriptorOptions(axis="coordination")


def f32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def test_single_molecule_units_and_exact_native_float():
    m = sx.Molecule.from_xyz_string("2\naxis\nC 0 0 0\nC 0 0 1.5\n")
    d = m.descriptors(sx.DescriptorOptions(donor_index=0, reference_index=1))
    assert d.sterimol.require().l == f32(3.2)
    assert d.sterimol.require().b1 == f32(1.7)
    assert d.sterimol.value.unit == "angstrom"
    assert d.buried_volume.error.code == "descriptor_error"
    with pytest.raises(sx.CalculationError):
        d.buried_volume.require()
    with pytest.raises(sx.InputError):
        sx.Molecule.from_xyz_string(XYZ + XYZ)


def test_order_failed_middle_and_missing_population():
    e = sx.Ensemble.from_sdf_string(SDF + "bad\n$$$$\n" + SDF)
    result = e.descriptors(COORD)
    assert [c.index for c in result.conformers] == [0, 1, 2]
    assert result.conformers[0].identifier == result.conformers[2].identifier
    assert result.conformers[1].geometry.error
    assert result.sterimol_b5.mean is None
    assert result.population is None
    assert result.sterimol_b5.boltzmann_mean is None
    available = e.descriptors(COORD, missing="available_cases")
    assert available.sterimol_b5.successful_count == 2
    assert [c.index for c in available.sterimol_b5.missing_conformers] == [1]
    assert available.regional_unavailable_reason


def test_thermodynamic_metadata_and_analytical_degeneracy():
    e = sx.Ensemble.from_sdf_string(SDF * 2)
    with pytest.raises(sx.PopulationError):
        e.descriptors(COORD, population=sx.Boltzmann(298.15, "enumerated states"))
    energies = [
        sx.Energy(0, "kcal_mol", "electronic", "test method", "calculation 1"),
        sx.Energy(
            0, "kcal_mol", "electronic", "test method", "calculation 2", degeneracy=3
        ),
    ]
    e.set_energies(energies)
    r = e.descriptors(
        COORD,
        population=sx.Boltzmann(
            298.15, "two states, second represents three microstates"
        ),
        ligand_id="test",
    )
    assert r.population.weights == pytest.approx((0.25, 0.75), abs=1e-15)
    assert r.population.origin["energy_kind"] == "electronic"
    assert r.sterimol_b5.boltzmann_mean is not None
    assert len(r.regions) == 36
    assert r.sensitivity["ranges"]["sterimol_b5"] == 0
    assert "standard_deviation" in r.sensitivity_csv
    before = e.records
    with pytest.raises(sx.PopulationError):
        e.set_energies([energies[0]])
    assert e.records == before
    r = e.descriptors(
        COORD, population=sx.SuppliedWeights((1, 3), "external empirical frequencies")
    )
    assert r.sterimol_b5.weighted_mean is not None
    assert r.sterimol_b5.boltzmann_mean is None


def test_structural_failure_keeps_per_conformer_descriptors():
    e = sx.Ensemble.from_xyz_string(XYZ + XYZ.replace("P 0", "N 0"))
    r = e.descriptors(COORD)
    assert len(r.conformers) == 2
    assert r.statistics_error.code == "ensemble_correspondence_error"
    assert r.conformers[0].sterimol.successful


def test_tables_keep_errors_and_float_values():
    pytest.importorskip("numpy")
    pytest.importorskip("pandas")
    r = sx.Ensemble.from_sdf_string(SDF + "bad\n$$$$\n").descriptors(COORD)
    n = r.to_numpy()
    p = r.to_pandas()
    assert len(n) == len(p) == 2
    assert n["sterimol_b5"][0] == r.conformers[0].sterimol.value.b5
    assert math.isnan(n["sterimol_b5"][1])
    assert n["geometry_error"][1]
    assert p.attrs["units"]["sterimol_b5"] == "angstrom"
    assert p.loc[1, "geometry_error"]


def test_batch_is_ordered_and_preserves_unreadable_files(tmp_path):
    p = tmp_path / "ligand.xyz"
    p.write_text(XYZ)
    rows = sx.calculate_files([p, tmp_path / "missing.xyz", p], COORD, threads=2)
    assert [r.path for r in rows] == [str(p), str(tmp_path / "missing.xyz"), str(p)]
    assert rows[0].conformers[0].sterimol.value == rows[2].conformers[0].sterimol.value
    assert rows[1].error.code == "input_error"
    with pytest.raises(sx.InputError):
        sx.calculate_files([p], threads=0)


def test_rdkit_preserves_coordinates_atom_order_bonds_and_conformer_ids():
    Chem = pytest.importorskip("rdkit.Chem")
    m = Chem.MolFromSmiles("P(C)(C)C")
    for identity, x in [(7, 1.234567891), (19, 1.456789123)]:
        c = Chem.Conformer(4)
        c.SetId(identity)
        for i, pos in enumerate(
            [(0, 0, 0), (x, 0, 0.45), (-0.7, 1.212, 0.45), (-0.7, -1.212, 0.45)]
        ):
            c.SetAtomPosition(i, pos)
        m.AddConformer(c, assignId=False)
    e = sx.Ensemble.from_rdkit(m)
    r = e.descriptors(COORD)
    assert [c.properties[0]["value"] for c in r.conformers] == ["7", "19"]
    assert r.conformers[0].geometry.value.coordinates_angstrom[1][0] == f32(1.234567891)
    assert r.conformers[0].geometry.value.bonds == ((0, 1), (0, 2), (0, 3))
    assert [c.GetId() for c in m.GetConformers()] == [7, 19]
    assert (
        sx.Molecule.from_rdkit(m, conformer_id=19).descriptors(COORD).identifier
        == "rdkit:19"
    )
    with pytest.raises(sx.InputError):
        sx.Molecule.from_rdkit(m)


def test_invalid_options_and_unknown_fields_fail():
    e = sx.Ensemble.from_xyz_string(XYZ)
    with pytest.raises(sx.InputError):
        e.descriptors({"axis": "coordination"})
    with pytest.raises(sx.InputError):
        e.descriptors(sx.DescriptorOptions(center=(math.nan, 0, 0)))
    with pytest.raises(sx.InputError):
        e.descriptors(missing="silently_drop")
    with pytest.raises(TypeError):
        sx.VolumeOptions(radius=3.5)


GOLDENS = json.loads((Path(__file__).parent / "data/native_goldens.json").read_text())


@pytest.mark.parametrize("golden", GOLDENS, ids=[g["id"] for g in GOLDENS])
def test_binding_preserves_frozen_native_descriptor_bits(golden):
    options = dict(golden["options"])
    options["volume"] = sx.VolumeOptions(**options["volume"])
    r = sx.Molecule.from_xyz_string(golden["xyz"]).descriptors(
        sx.DescriptorOptions(**options)
    )
    expected = golden["expected"]
    for name in ("l", "b1", "b5"):
        assert getattr(r.sterimol.require(), name) == f32(
            expected["sterimol"]["value"][name]
        )
    for name, value in expected["buried_volume"]["value"]["summary"].items():
        assert getattr(r.buried_volume.require().summary, name) == f32(value)
    for name, value in expected["pyramidalization"]["value"].items():
        assert getattr(r.pyramidalization.require(), name) == f32(value)
