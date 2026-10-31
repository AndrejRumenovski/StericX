# Advanced Python descriptor methods

The public Python API now exposes the validated density mode and fixed-center
radius profiles, alongside ordered conformer descriptors, Rust ensemble
statistics, population calculations and sensitivity reports. All scientific
calculations call the shared Rust engine. Python adapts inputs and typed outputs.

```python
from stericx import Ensemble, DescriptorOptions, Energy, Boltzmann

ensemble = Ensemble.from_sdf("ligands.sdf")
options = DescriptorOptions(axis="coordination", donor_index=0)
results = ensemble.descriptors(options, ligand_id="ligand")
print(results.sterimol_b5.mean, results.sterimol_b5.minimum)
print(results.sterimol_b5.maximum, results.sterimol_b5.boltzmann_mean)
# Last value is None until complete, justified energy metadata is supplied.
# Never infer energies or units from SDF property names or titles.
print(results.sensitivity["ranges"]["sterimol_b5"])
results.to_pandas().to_csv("conformers.csv", index=False)
```

`Energy` requires value, unit, kind, model and provenance. Gibbs free energies
also require their evaluation temperature. `Boltzmann(T, assumptions)` validates
all energies and assumptions in Rust; `SuppliedWeights` remains distinguishable
from energy-derived populations. Failed conformers remain in the source order;
`missing="available_cases"` explicitly opts into conditional statistics and
reports retained population mass. The default requires completeness per field.

```python
from stericx import DensityGrid

grid = DensityGrid.from_cube(
    "density.cube",
    field="total_electron_density",
    density_unit="electrons_per_bohr_cubed",
    provenance="method/basis, charge, calculation and file identifier",
)
profile = grid.buried_volume_scan(
    [1.0, 1.5, 2.0],
    center_angstrom=(0.0, 0.0, 0.0),
    threshold_electrons_per_bohr3=0.0016,
)
for sample in profile.samples:
    print(sample.requested_radius_angstrom, sample.outcome)
```

The threshold above is a documented convention, not a universal physical
constant. The cube must actually contain total electron density in the declared
units and cover every requested sphere. No field, threshold, missing density or
center is inferred. Density results report integer counts, cell volume, grid
buried volume, sampled and analytical sphere volumes, sampling error and percent.
Hard-sphere and density absolute volumes use different finite normalizations.
See [methodology](../../methods/electron_density_buried_volume.md).

`Molecule.buried_volume_scan(radii_angstrom, options=...)` and
`Ensemble.buried_volume_scan(..., conformer_index=...)` vary the sphere radius,
not center distance. They preserve order, duplicates and individual failures,
and report actual f32 radii for native hard-sphere calculations. The known
R=2.5 Å boundary difference with DBSTEP remains documented; the wrapper does not
change the scientific algorithm.

**71 Python tests pass.** Advanced tests include three independently validated
Ne density radii, exact node counts and descriptor values, explicit units and
threshold failures, incomplete domains, profile identity/order/failures, and
Gibbs-temperature mismatch rejection. The original DBSTEP Ne cube and its MIT
license accompany the tests. Wheel/OS/version distribution admission is separate.
