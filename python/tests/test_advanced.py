import json
from pathlib import Path

import pytest
import stericx as sx

DATA = Path(__file__).parent / "data"
DENSITY = json.loads((DATA / "density_goldens.json").read_text())


def grid(path=None):
    return sx.DensityGrid.from_cube(
        path or DATA / "Ne_medium.cube",
        field="total_electron_density",
        density_unit="electrons_per_bohr_cubed",
        provenance="official DBSTEP total SCF fixture; method not inferred",
    )


@pytest.mark.parametrize("expected", DENSITY, ids=lambda r: str(r["radius_angstrom"]))
def test_density_binding_matches_independent_reference_goldens(expected):
    g = grid()
    r = g.buried_volume(
        center_angstrom=(0, 0, 0),
        radius_angstrom=expected["radius_angstrom"],
        threshold_electrons_per_bohr3=0.0016,
    )
    assert r.sphere_nodes == expected["sphere_nodes"]
    assert r.occupied_nodes == expected["occupied_nodes"]
    assert r.buried_volume_angstrom3 == expected["buried_volume_angstrom3"]
    assert r.percent_buried_volume == expected["percent_buried_volume"]
    assert r.grid["bohr_to_angstrom"] == 0.529177210544
    assert r.method == "electron_density_node_quadrature"


def test_density_requires_explicit_field_units_threshold_and_coverage():
    with pytest.raises(TypeError):
        sx.DensityGrid.from_cube(DATA / "Ne_medium.cube")
    with pytest.raises(sx.InputError):
        sx.DensityGrid.from_cube(
            DATA / "Ne_medium.cube",
            field="orbital",
            density_unit="electrons_per_bohr_cubed",
            provenance="x",
        )
    with pytest.raises(sx.DensityError):
        grid().buried_volume(
            center_angstrom=(0, 0, 0),
            radius_angstrom=3.5,
            threshold_electrons_per_bohr3=0.0016,
        )
    with pytest.raises(sx.DensityError):
        grid().buried_volume(
            center_angstrom=(0, 0, 0),
            radius_angstrom=1.0,
            threshold_electrons_per_bohr3=0.0,
        )


def test_density_profile_preserves_every_position_and_duplicate():
    p = grid().buried_volume_scan(
        [1.0, 1.5, 2.0, 10.0, 1.0],
        center_angstrom=(0, 0, 0),
        threshold_electrons_per_bohr3=0.0016,
    )
    assert p.definition == "sampling_sphere_radius_at_fixed_center"
    assert [s.requested_radius_angstrom for s in p.samples] == [
        1.0,
        1.5,
        2.0,
        10.0,
        1.0,
    ]
    assert p.samples[3].outcome.error
    assert p.samples[0].outcome.value == p.samples[4].outcome.value
    for sample, expected in zip(p.samples[:3], DENSITY, strict=True):
        assert sample.outcome.value.occupied_nodes == expected["occupied_nodes"]
        assert sample.outcome.value.sphere_nodes == expected["sphere_nodes"]


def test_hard_sphere_profiles_use_existing_kernel_and_fixed_center():
    golden = json.loads((DATA / "native_goldens.json").read_text())[0]
    opts = dict(golden["options"])
    opts["volume"] = sx.VolumeOptions(**opts["volume"])
    opts = sx.DescriptorOptions(**opts)
    m = sx.Molecule.from_xyz_string(golden["xyz"])
    p = m.buried_volume_scan([3.5, -1.0, 2.5, 4.5, 3.5], options=opts)
    assert p.samples[0].outcome.value == m.descriptors(opts).buried_volume.value
    assert p.samples[1].outcome.error
    assert p.samples[0].outcome.value == p.samples[-1].outcome.value
    assert all(
        s.outcome.value.center == p.center_angstrom
        for s in p.samples
        if s.outcome.successful
    )
    assert p.samples[2].outcome.value.summary.buried_volume == pytest.approx(
        16.43008804321289, abs=0.0
    )
    with pytest.raises(sx.CalculationError):
        m.buried_volume_scan([float("nan")], options=opts)


def test_free_energy_temperature_mismatch_is_not_electronic_weighting():
    xyz = "2\ntest\nC 0 0 0\nC 0 0 1.5\n"
    e = sx.Ensemble.from_xyz_string(xyz + xyz)
    e.set_energies(
        [
            sx.Energy(
                0,
                "kcal_mol",
                "gibbs_free_energy",
                "same method",
                "source",
                temperature_k=298.15,
            )
        ]
        * 2
    )
    with pytest.raises(sx.PopulationError):
        e.descriptors(population=sx.Boltzmann(300, "explicit states"))
