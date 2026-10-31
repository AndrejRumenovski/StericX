# Conformer descriptor contract

The Rust API is `steric_x::ensemble::{Ensemble, DescriptorOptions}`:

```rust
use steric_x::ensemble::{DescriptorOptions, Ensemble, SterimolAxis};

let ensemble = Ensemble::from_sdf("ligands.sdf")?;
let options = DescriptorOptions {
    axis: SterimolAxis::Coordination,
    ..Default::default()
};
let results = ensemble.descriptors(&options);
for conformer in &results.conformers {
    println!("{} {:?}", conformer.index, conformer.sterimol);
}
# Ok::<(), String>(())
```

Each source record produces exactly one ordered result. `index` is the unique
zero-based positional key; `identifier` preserves the original title, including
duplicates and blanks. Original record text, source path/label, format, parsed
coordinates/radii/bonds, and raw SDF properties are retained. Duplicate property
names are not overwritten. XYZ supports concatenated frames. An invalid XYZ
count prevents trustworthy recovery of subsequent boundaries: the remaining
text becomes one explicit error record, never a guessed set of conformers.

Each descriptor family has a tagged success/error outcome. A failed parse does
not remove the record. An invalid buried-volume configuration does not erase
valid pyramidalization or explicit bond-axis Sterimol. Explicit donor/reference
indices also allow bond-axis Sterimol for molecules lacking a trivalent donor;
unavailable regional volume and pyramidalization remain errors.

Units: Sterimol Å; total/near/far/quadrant/octant volume Å³; percent volume %;
pyramidalization P dimensionless and alpha degrees. Native coordinates and
descriptors retain f32 precision. No scientific constants or default arithmetic
were changed. Defaults match the descriptor CLI (bond axis, 2.28 Å virtual-center
distance); the lower-level volume config's separate 2.1 Å default is unchanged.
Bond-axis L is raw support; coordination-axis L adds the documented 0.40 Å
correction unless explicitly overridden. Optional endpoint-inclusive angular
sampling reuses the existing coordination-axis kernel.

`buried_volume.value.summary` contains all nine existing scalar fields.
`orientations` contains all three original plane estimates, their plane atom
indices and source-coordinate X/Y/Z vectors. Quadrants are (++), (-+), (--),
(+-); octants follow those four at z ≥ 0 then at z < 0. Negative Z faces the
donor. Zero coordinates belong to positive regions. Regions have independent
finite-grid normalizations and need not sum to the separately estimated total
on odd grids. These are molecular frames, not fixed laboratory quadrants;
cross-conformer regional comparisons require corresponding atom identities.

Loading does not assert that arbitrary SDF records are the same chemical species
or infer atom mappings. All rows remain inspectable. Ensemble statistics and
thermodynamic populations have additional contracts in the following steps.

Energies are optional typed metadata attached using `set_energy`. The caller
declares units, electronic/force-field/Gibbs interpretation, common model,
source provenance, degeneracy and free-energy evaluation temperature. Raw SDF
properties are preserved even when their units are unknown; they are never
silently converted into energies. No populations are inferred at this step.

Validation: `tests/conformer_descriptors.rs` adds seven ingestion, analytical
sphere-support, independence-of-failures, regional-identity and metadata tests.
All 329 Rust tests pass. Existing per-file calculations now share their donor
selection routine with the library API. The complete scientific replay remains
a separate gate; agreement with a native API test is not an independent physical
validation of its descriptor convention.
# Ensemble statistics

`EnsembleDescriptors::statistics(&StatisticsOptions::default(), population)`
returns 14 named scalar descriptors and, when ordered explicit connectivity is
available for every record, all molecular-frame quadrant and octant statistics.
Each field carries its unit, total and successful counts, failed conformer IDs,
and all conformer IDs tied at each extreme. Atom element order and any supplied
adjacency must agree. This check does not establish bond orders, charge,
stereochemistry, search completeness or chemical identity; the caller must
provide a corresponding ensemble. Regional keys identify the original plane
atom and zero-based bin. They do not represent a laboratory-fixed direction.

For finite descriptor values \(x_1,\ldots,x_n\), the arithmetic mean is
\(\sum_i x_i/n\), minimum and maximum are the order extrema, and standard
deviation is \(\sqrt{\sum_i(x_i-\bar x)^2/n}\) (population standard deviation,
`ddof=0`). The median is the middle sorted value for odd n, or the arithmetic
mean of the two middle values for even n. A singleton has standard deviation
zero. These are descriptive statistics of the supplied states, not uncertainty
estimates or a correction for uneven conformer sampling. Reductions use f64
over the unchanged native f32 descriptor values.

No weights are inferred from an SDF property, missing energies, or the number
of conformers. Without an explicitly supplied `Population`, both weighted
fields are unavailable. `Population::supplied(weights, provenance)` accepts
finite nonnegative weights with positive total and normalizes them in input
order. Its weighted mean is \(\sum_i w_i x_i\); this is never labelled a
Boltzmann mean. The caller must associate supplied weights with the complete
ordered record list. Energy-derived populations additionally check that their
ordered energy metadata matches the descriptor records.

The default `RequireComplete` policy leaves **all** statistics for a field
unavailable if any conformer lacks that field. Other descriptor families can
remain available. Opting into `AvailableCases` explicitly conditions on the
successful subset S: arithmetic statistics use only S, and weighted means use
\(\sum_{i\in S}w_i x_i/\sum_{i\in S}w_i\). Results retain failed IDs,
`complete=false` and the original retained population mass. Zero retained mass
leaves the weighted mean unavailable. This conditional result must not be
reported as the full ensemble. Nonfinite values count as failed, not zeros.

`Ensemble::boltzmann_population(T, assumptions)` requires declared, valid energy
metadata for every record. It uses
\(w_i=g_i\exp[-(E_i-E_{min})/(RT)]/\sum_j g_j\exp[-(E_j-E_{min})/(RT)]\)
with log shifting for stable normalization. All energies must share unit, kind
and model. Temperature is positive finite kelvin; Gibbs free energies must have
been evaluated at that exact temperature. Units are kcal/mol, kJ/mol, hartree
per molecule, or eV per molecule. The existing gas constant is unchanged;
conversions use exact SI Avogadro/elementary-charge constants, thermochemical
1 kcal = 4184 J, and CODATA 2022 hartree energy (see source citation in code).

Electronic and force-field energy weights are explicitly approximations, not
free-energy equilibrium populations. The supplied assumptions must address
state membership and degeneracy. Equal-energy records remain separate and
normally each have degeneracy one. An explicitly aggregated state can have
integer multiplicity g; do not count multiplicity both in enumerated records
and g, or again if its entropy is already included in the free energy. This API
does not invent historical energies, establish equilibration, or infer missing
thermochemical/solvent contributions. Independent thermodynamic validation is
recorded separately before scientific release.

