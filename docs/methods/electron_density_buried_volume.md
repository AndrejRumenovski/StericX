# Electron-density buried volume: definition and validation contract

This is a separate descriptor from the union of scaled atomic hard spheres.
Implementation is scientifically feasible for a restricted scalar cube pathway:
the official DBSTEP source, real density fixtures and independent mathematical
checks are available. Public scientific support is conditional on the validation
results, not on this methodology document alone.

## Scientific interpretation and threshold

For total electron density ρ(x), threshold τ, fixed center c and radius R, define
occupation as `ρ(x) > τ`, and the sampling region as `||x-c|| <= R`. The continuous
target is the integral of this occupation indicator over that sphere. This
measures the volume of the chosen density superlevel set inside a declared
region. It is not integrated electron count, an atomic-radius volume, a solvent
accessible surface, or an experimentally established catalyst-selectivity model.

Alibakhshi and Schäfer (2024) compared iso-density surfaces against thermodynamic
effective surfaces for 104 molecules. They report approximately 1.6% mean unsigned
percentage deviation and correlation 0.995 at 0.0016 atomic units, using their
specified quantum-chemical protocol. That supports a documented threshold
convention; it does not establish a universal threshold or validate every
metal-ligand buried-volume application. Other thresholds define other surfaces.
The API therefore requires an explicit positive threshold in electrons/Bohr³,
an assertion that the scalar field is total electron density, and provenance for
the quantum calculation. It never guesses field identity from a cube filename.
[Primary paper](https://doi.org/10.1038/s41467-024-50408-8),
[open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11271626/).

Bader et al.'s earlier atomic-volume work is bibliographically relevant:
[Properties of atoms in molecules: atomic volumes (1987)](https://doi.org/10.1021/ja00260a006).
Its full text was not available in this session; no implementation detail is
attributed to an unread paper.

## Reference implementation

Use unmodified DBSTEP 2.2.0 at
[`4365353e0a2f24d05336e619759e350bc630a3cb`](https://github.com/patonlab/DBSTEP/tree/4365353e0a2f24d05336e619759e350bc630a3cb).
`parse_data.CubeParser` reads scalar cube values in x-slow/y-middle/z-fast order.
`sterics.occupied_dens` applies strict `density > isoval`; `buried_vol` counts
occupied nodes and all sampling nodes inside the sphere, returning their ratio
as a percentage. The public API's actual default isoval is 0.0016, despite an
older 0.002 function docstring. Explicit settings avoid relying on that mismatch.
[Parser](https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/parse_data.py),
[occupation/integration](https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/sterics.py),
[API settings](https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/Dbstep.py).

The pinned density pathway treats x spacing as a common spacing and constructs
axis-aligned grids. It cannot independently validate general oblique or
anisotropic integration. It also resizes the denominator grid if the sphere
plus 10% margin extends outside the original grid; this can alter lattice
spacing. Native verified analysis will reject incomplete domains and never
extrapolate, zero-pad, resize or interpolate density. Reference fixtures must
cover the 10% margin so public DBSTEP does not resize them.

## Input and coordinate contract

Read formatted Gaussian-compatible scalar cube files with positive atom count,
positive grid counts, one scalar per node and all coordinates in Bohr. Preserve
two comment lines, source identifier, atomic numbers/charges/coordinates, origin,
three full grid increment vectors, dimensions and original density values.
An optional scalar-count header must equal one. Reject orbital cubes, multiple
fields, ambiguous signed-count unit extensions, incomplete/excess data,
nonfinite numbers and negative total density. Fortran D exponents may be read
without changing numerical meaning. The caller must explicitly declare total
density in electrons/Bohr³. Other fields can have the same file format.
[Official VMD cube-format documentation](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html).

Initial analysis supports only positive axis-aligned isotropic grids; the reader
may preserve other nonsingular vectors but analysis rejects them explicitly.
There is no silent flattening of off-diagonal vectors or conversion of ambiguous
units. Center coordinates are explicit in the source coordinate system, in Å;
a caller may select a stored atom position. No molecular rotation is needed.

Convert Bohr to Å using CODATA 2022 `a0 = 0.529177210544 Å`, explicitly recorded
in results. DBSTEP uses `0.529177249 Å`; retain and quantify that convention
difference rather than changing either implementation to fit. Density values
and threshold remain in electrons/Bohr³, so coordinate conversion never changes
the occupation criterion.
[NIST CODATA 2022 constants](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).

## Numerical definition

For accepted grid vectors a,b,c in Å, cell volume is `|det(a,b,c)|`. Sample
coordinates are `origin + i*a + j*b + k*c`, with k fastest. Count nodes with
squared distance `<= R²`, and among these those with density strictly `> τ`.
Report both integer counts and:

* grid buried volume = `N_occupied * cell_volume` in Å³;
* sampled sphere volume = `N_sphere * cell_volume` in Å³;
* percent buried volume = `100 * N_occupied/N_sphere`;
* analytic sphere volume = `4*pi*R³/3` and the relative sampling-volume error.

This deterministic node quadrature has boundary/discretization error. No partial
voxel correction or hidden analytic-sphere normalization is applied to the grid
buried volume. A sphere must fit within the node-domain bounds on all axes;
empty sampling sets are errors. All coordinates, radii and thresholds must be
finite. Derived overflows or zero cell volume are errors. Refinement studies
must keep the physical density, center, threshold and region fixed.

## Radius profiles

DBSTEP's scan varies the sampling-sphere **radius** about a fixed center. It does
not move the center along a bond. Name the public argument `radii_angstrom` and
record this definition with every profile. Each requested position is an
independent calculation with its own grid/counts and result/error; preserve
caller order, duplicate radii and failures. Native hard-sphere profiles retain
the existing integration algorithm and defaults at each radius. Density profiles
reuse the supplied density grid and reject radii outside its domain. Shell
profiles and center-displacement scans are outside this initial API.

## Admission tests fixed before implementation

Use all three declared official scalar-density fixtures: `Ne_medium.cube`,
`benzene_coarse.cube`, `CH2CMe3_100.cube`, first-atom centers, R=1.0, 1.5, 2.0 Å,
and explicit τ=0.0016 electrons/Bohr³. Check domain/margin sufficiency first;
retain unavailable comparisons. Freeze fixture and reference hashes. Compare
unmodified public DBSTEP and its independent occupied-density/support functions.
Every common-domain node count must agree exactly; arithmetic percentages and
volumes have a fixed `64*epsilon_f64*max(1,|value|)` comparison allowance after
explicitly accounting for the recorded Bohr conversion. Any membership difference
must be investigated, not absorbed into an empirical tolerance.

Also test analytically enumerable grids, strict threshold equality, all-empty
and all-full occupation, translated centers/origins, x/y/z indexing, malformed
and incompatible inputs, domain failures, determinism and convergence against a
known spherical superlevel set. Use independently generated analytic densities
at several spacings; an analytic geometric boundary layer provides an a priori
quadrature error bound. The real cubes establish algorithm agreement, not
quantum-density accuracy or universal threshold validity. Only validated modes
may be exposed in Python or presented as release features.
