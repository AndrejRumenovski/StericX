# StericX

**A Rust molecular steric-descriptor engine with explicit conformer analysis and
Python bindings.** StericX calculates Sterimol, hard-sphere buried volume and
pyramidalization, and provides a separately defined electron-density-grid
buried-volume mode. Independent comparisons, failed cases and reproducible
inputs accompany the implementation.

[![CI](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license)

## Install and calculate

The Python package is **0.4.0.dev0**, backed by the existing 0.3.0 Rust engine.
It has not been published to PyPI. From a checkout, with Rust installed:

```bash
python -m pip install ./python
# Optional table and RDKit adapters:
python -m pip install './python[tables,rdkit]'
```

A portable Linux x86-64 wheel has been built and installed in isolated Python
3.11, 3.12, 3.13 and 3.14 environments. Each installation passes 71 package tests;
a separate source-archive installation also passes. Actual CI also passes all 12 Python-version/OS combinations on Linux x86-64,
Windows x86-64 and macOS arm64, plus source installation. The
[fresh installation and CI replay](docs/validation/dbstep/resolution_advanced/RESULTS.md)
preserves exact artifacts and source identities.

```python
from stericx import Molecule, DescriptorOptions

mol = Molecule.from_file("ligand.xyz")
options = DescriptorOptions(axis="coordination", donor_index=0)
result = mol.descriptors(options)

print(result.sterimol.require())  # L, B1, B5 in Å
print(result.buried_volume.require().summary)  # volume in Å³, fraction in %
print(result.pyramidalization.require())
```

Choose the actual donor index for your geometry; indices are zero-based.
Coordination-axis Sterimol uses the inferred virtual center and the documented
+0.40 Å L correction. Bond-axis calculations require a meaningful attachment
axis. Undefined descriptors return structured errors independently, so a
failure does not silently remove a conformer.

For the standalone CLI:

```bash
cargo build --release --locked
./target/release/stericx descriptors ligand.xyz --sterimol-axis coordination --format json
./target/release/stericx --help
```

The [retained CLI guide](LEGACY_WORKFLOWS.md) covers database construction,
comparison, batch files and legacy research workflows.

## Conformers are explicit data

```python
from stericx import Ensemble, DescriptorOptions

ensemble = Ensemble.from_sdf("ligands.sdf")
results = ensemble.descriptors(DescriptorOptions(axis="coordination"))

print(results.sterimol_b5.mean)
print(results.sterimol_b5.minimum, results.sterimol_b5.maximum)
print(results.sterimol_b5.standard_deviation)  # population SD, divisor n
print(results.sterimol_b5.median)
print(results.sterimol_b5.boltzmann_mean)  # None without valid energies

for conformer in results.conformers:
    print(conformer.index, conformer.identifier, conformer.sterimol)

results.to_pandas().to_csv("conformers.csv", index=False)
```

Every input record keeps its order, identity, source geometry, properties and
failure status. Supported conformers expose L/B1/B5, total and percent buried
volume, quadrant and octant volumes, and pyramidalization. Regional results
carry their molecular frame; ensemble regional statistics require consistent
explicit atom mapping/connectivity. XYZ files alone do not establish that map.

Statistics describe the supplied finite ensemble. The default requires complete
coverage for each descriptor; an explicitly selected available-case policy
reports the retained subset and population mass. Sensitivity reports include
range, standard deviation and all conformers tied at each extreme. No universal
flexibility score is imposed.

Boltzmann weighting requires explicit comparable energies, units, provenance,
energy kind, temperature, state degeneracy and assumptions. Electronic or
force-field energy weighting is labelled accordingly; it is not silently called
a free-energy population. Missing energies never produce invented uniform
Boltzmann populations. The
[API and mathematical definitions](docs/paton_recommendations/CONFORMER_API.md)
and [executed notebook](python/examples/stericx_research.ipynb) show the complete
workflow, including recorded MMFF energies, RDKit, NumPy/pandas, batch errors and
sensitivity exports.

## Density grids and radius profiles

```python
from stericx import DensityGrid

grid = DensityGrid.from_cube(
    "density.cube",
    field="total_electron_density",
    density_unit="electrons_per_bohr_cubed",
    provenance="Your calculation method, basis, state and source file",
)
profile = grid.buried_volume_scan(
    [2.0, 2.5, 3.0],
    center_angstrom=(0.0, 0.0, 0.0),
    threshold_electrons_per_bohr3=0.0016,
)
```

This scan changes the **sphere radius at a fixed center**. Molecular
`buried_volume_scan` follows the same radius definition. Density occupation is
strictly density > the explicit threshold; volume is occupied-node count times
cell volume. The threshold shown is DBSTEP's documented convention, not a
universal physical constant. Density and hard-sphere volumes are distinct
quantities.

The cube reader preserves origin, vectors, dimensions, atom coordinates and
field declaration. The validated analysis supports positive, axis-aligned,
isotropic scalar density grids. Unsupported orbital/multifield/unit conventions
and spheres extending beyond the supplied grid fail explicitly. No padding or
interpolation is hidden. Read the
[method and limitations](docs/methods/electron_density_buried_volume.md) and
[independent density validation](docs/validation/dbstep/DENSITY_RESULTS.md).

## What independent validation establishes

| Evidence | Result and limit |
|---|---|
| Legacy numerical preservation | All 128,021 observations in 13 lanes match the task's frozen baseline byte for byte. This proves unchanged outputs, not physical validity. |
| DBSTEP Sterimol, supplied common radii | 61 cases; common signed-L subset maximum error 1.10e-6 Å; B5 maximum 4.27e-7 Å. Default finite-angle B1 maximum difference is 0.01891 Å. Radii and signed/absolute L conventions remain distinct. |
| DBSTEP hard-sphere volume | All 69 cases / 207 frames have identical atomic occupancy masks on the same supplied points; all 138 public volume/percentage comparisons meet the original arithmetic limits. The historical failed lane clipped the sample domain again using a different sphere predicate; its failures remain preserved. |
| Continuous-volume accuracy | 320 comparisons on 20 fixed geometries separate finite-grid semantics from integration error. For ligand 401 at R=2.5 Å, default-grid errors are provably at least 0.034–0.058 Å³ across its three frames. Resolution must be checked for the intended use. |
| DBSTEP density mode | All 12 real-grid/analytic comparisons have identical occupied and sphere node counts; explicit unit conversion and volume checks pass. Grid-resolution error remains separate. |
| Morfeus conformer study | 784 direct conformer fields and 770 ensemble reductions on the fixed 11-ligand/56-conformer study; all comparisons retained. Ensemble treatment is not uniformly superior against historical targets. |
| Thermodynamic populations | Independent high-precision analytical cases cover units, temperature, degeneracy, missing information and extreme energy separation. |

Full residuals, MAE/RMSE/max/median errors, slopes/intercepts and acceptance
policies are in the [DBSTEP reports](docs/validation/dbstep/),
[conformer study](docs/paton_recommendations/) and
[legacy replay receipt](docs/paton_recommendations/validation/final_replay_comparison.json).
The [scientific audit](docs/scientific_accuracy_audit/SCIENTIFIC_ACCURACY_AUDIT.md)
and [remediation](docs/scientific_remediation/SCIENTIFIC_REMEDIATION.md) retain
previous defects, corrections and unfavorable results.

The [DBSTEP root-cause study](docs/validation/dbstep/ROOT_CAUSE_ANALYSIS.md)
documents the corrected comparison, an independently confirmed DBSTEP
numerator/denominator inconsistency, and the approved separate validation lanes.
No Rust scientific output, radius or default setting changed to obtain agreement.

**Direct descriptor agreement and historical Kraken reproduction answer different
questions.** The [Kraken investigation](docs/kraken_residual_forensics/KRAKEN_RESIDUAL_FORENSICS.md)
compares 1,541 ligands and 31,611 available conformers, with complete residuals.
Exact identity with the original historical geometries is established for zero
ligands; historical conformer energies are unavailable. High correlation on
available inputs does not resolve those provenance gaps.

Native geometry kernels retain f32 arithmetic and the default 360-direction B1
scan; exact rotational invariance is not claimed for this finite scan. The
lower-level buried-volume library retains its 2.1 Å center-distance default;
the CLI and new descriptor API use 2.28 Å. Regional finite grids and analytic
sphere normalization have documented conventions. No scientific constants were
changed to improve reference agreement.

For hard-sphere volume, `VolumeOptions(density=...)` controls nominal grid-cell
volume in Å³; smaller values generally refine the lattice. It is not an accuracy
tolerance. The [fixed convergence study](docs/validation/dbstep/resolution_analysis/convergence_v2/CONVERGENCE_REPORT.md)
includes analytical errors, independent continuous-volume bounds and all
nonmonotonic refinement steps. Arithmetic agreement on identical nodes does
not establish continuum accuracy.

## Performance and reproducibility

Rust supplies the shared scientific engine; PyO3 binds it without duplicating
chemistry in Python. Independent files use ordered Rayon processing. Performance
claims require fixed geometries/settings, numerical admission and retained raw
samples. The [Paton performance report](docs/paton_recommendations/performance/PERFORMANCE.md)
compares the extended engine with this task's frozen baseline, including
allocation counts, Python overhead and a separately retained I/O-contention study.

The earlier [Morfeus benchmark](docs/benchmarks/morfeus_current/BENCHMARK.md)
is scoped to nine finite-grid volume outputs on 56 repeated geometries. Its
reported speedups do not establish an all-descriptor or universal speed ratio.
The [optimization report](docs/performance/SCIENTIFICALLY_EXACT_OPTIMIZATION.md)
retains accepted and rejected changes.

Start with the [final implementation report](docs/paton_recommendations/FINAL_REPORT.md),
[complete scientific validation](docs/paton_recommendations/SCIENTIFIC_VALIDATION.md),
[research examples](python/examples/README.md), and
[scientific replay instructions](docs/scientific_remediation/REPRODUCE.md).
Legacy reaction modelling remains experimental: unfavorable validation, failed
folds and missing prospective outcomes are preserved. Descriptor validation does
not establish catalyst-prediction accuracy.

## Citation and license

Cite original reference methods and datasets when used: the
[Kraken platform](https://doi.org/10.1021/jacs.1c09718),
[Morfeus](https://github.com/digital-chemistry-laboratory/morfeus), and
[DBSTEP](https://github.com/patonlab/DBSTEP). Density-method primary references
are listed in the [method document](docs/methods/electron_density_buried_volume.md).
StericX citation metadata is in [CITATION.cff](CITATION.cff); the existing
[0.3.0 archive](https://doi.org/10.5281/zenodo.22309199) predates these development
extensions. No new archived release is implied.

Licensed under [Apache-2.0](LICENSE-APACHE) or [MIT](LICENSE-MIT), at your option.
Third-party fixtures retain their source licenses and provenance.
