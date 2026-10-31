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
