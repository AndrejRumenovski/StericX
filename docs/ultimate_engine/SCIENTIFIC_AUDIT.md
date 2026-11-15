# Scientific audit

Status: source audit and independent calculation replays are complete; final
performance admission is separate. This document distinguishes confirmed defects,
tested behavior and reference conventions. Existing
scientific-remediation and DBSTEP failure evidence remains unchanged.

## Confirmed correction S1: R² depends on target units

`model::fit::metrics` discarded target variance below 1e-12;
`model::evaluation::score_frozen_predictions` discarded variance below machine
epsilon. Both thresholds have absolute squared target units. They incorrectly
returned `None` for small, nonconstant targets.

Independent witness: actual targets (-s,+s), predictions (-s/2,+s/2).
The mean is zero, RSS=s²/2 and TSS=2s², so R²=1-RSS/TSS=3/4. Using dyadic
s=2^-30 makes every operation exactly representable. Both original functions
failed this test; see [r2_baseline_failure.log](r2_baseline_failure.log).

The correction tests TSS>0 and evaluates the ratio lazily. Constant targets
still return `None`; poor predictions can still have negative R². No fit
regularization, physical constant, prediction or scientific tolerance changes.
Regression coverage checks s=1, 2^-30 and 2^-100, and the fit helper also checks
constant targets and the exact score -3. All 58 model unit tests and six public
training API tests pass, including the unchanged golden training report and
bit-for-bit comparison against the former CLI composition. Raw validation logs
are in `.stericx/ultimate_engine/r2_{model_tests,training_api}.log`.

## Audit coverage and conventions

The table summarizes fresh evidence; detailed counts and limits are in
[VALIDATION_REPORT.md](VALIDATION_REPORT.md).

| Calculation family | Definition and audit focus | Current evidence/status |
|---|---|---|
| Sterimol | Raw L versus external +0.40 Å convention; signed support B1 on 360 directions; radial B5; f64 projection with f32 inputs/results; explicit dummy and axis | 61 legacy DBSTEP cases; 39 high-precision continuous-angle cases and matched external comparisons |
| Hard-sphere buried volume | Union membership on the exact f32 sphere lattice, radii scaling, H exclusion, three donor planes, independently normalized regions | 69 common-node cases/207 frames; 31,721 independent full-volume references; analytical convergence |
| Quadrants/octants/hemispheres | Zero belongs to the positive side; each region uses its own finite-grid denominator; near is negative local Z | All bins in baseline replay captured; regional sums need not equal total on odd grids |
| Coordination geometry | Negative sum of unit bonds, trivalent topology, planar sign ambiguity requires explicit center | Existing ambiguity/degeneracy gates preserved; geometric lone pair is not an electronic orbital |
| Pyramidalization | Plane geometry and angular descriptors, proper degeneracy rejection | Fresh geometry and fixed 2,517-variant Kraken union replay, with three reference neighbor/precision conventions |
| Ensembles/statistics | Identity/order, minima/maxima/ties, weighted means/distributions, missing-data provenance | 11 ligands/56 conformers plus independent aggregation equations and exact legacy replay |
| Populations/thermodynamics | Declared energy units, reference-energy subtraction, temperature, stable exponential normalization | 42 valid/11 invalid high-precision population cases; missing-energy policy preserved |
| Electron density | Explicit total density, Bohr conversion, grid transforms, strict density > threshold, cell-volume integral | 12 fresh independent DBSTEP/analytic Gaussian cases pass; convergence remains discretization-dependent |
| Radius profiles/sensitivity | Fixed center, requested radius order/duplicates and per-radius failures | 20 positions, 16 reference comparisons and four expected errors; historical strict gate remains false |
| Eyring/ee/reaction calculations | Energy/temperature units, signed conventions, stable extreme limits and uncertainty | 744 scalar and 51 aggregation/rejection checks, 14 CLI cases; scoped pass |
| Predictive models | Feature conventions, AVX2/scalar reduction differences, nested versus fixed-feature validation, intervals, extrapolation | S1 corrected; 7,175 model plus 449 supplement comparisons pass |
| Input/radii/storage | Coordinate precision, element/radius policies, parser errors, record serialization | Frozen source/data hashes and exact parser replay retained; public radius fallbacks are conventions, not fitted constants |

## Calculation-by-calculation review

### Sterimol and coordinate construction

The legacy scalar quantity is forward axial support L, minimum sampled signed
transverse support B1, and maximum radial support B5. The attachment atom is
excluded in the real-atom axis API; the dummy-axis API has its own explicit
origin. A maximum with zero belongs to B1's support definition. It is not the
minimum radial atom distance or the smallest centered enclosing cylinder.
The shortest-arc transverse frame and one-degree phase are part of the finite
algorithm. A rigid rotation can change sampled B1 through angular phase; such
a difference must not be called a failed continuous rotational invariant.

Axis/plane construction uses wider intermediates where required for very small
or large represented coordinates. Zero-length axes, malformed indices,
nonfinite positions and invalid radii return errors. The source coordinates
and ordinary descriptor outputs remain binary32. Translation cannot restore
coordinate differences already lost during binary32 input conversion. Exact
atom permutations with correctly remapped identities are a stronger discrete
test than rotating, rounding, and then expecting bitwise equality.

The continuous f64 methods in [NUMERICAL_ACCURACY.md](NUMERICAL_ACCURACY.md)
remove angular sampling error explicitly. They do not change the legacy scan,
radii, axis or the optional external 0.40 Å L convention. Morfeus samples both
B1 and B5; DBSTEP's default alignment rounds to eight decimal places and its
classic L uses an absolute axial convention. These differences were checked
against pinned implementation source and preserved failed comparisons.

### Hard spheres, regional measures and center assumptions

The original lattice generator, point order, f32 sphere inclusion and atomic
predicate `(dx*dx + dy*dy) + dz*dz <= radius_squared` remain authoritative for
the default method. SIMD uses separate multiply/add operations, ordered
comparisons and integer OR/popcount reductions. It does not replace the test
with `dz² <= r² - dx² - dy²`, which can differ near floating-point absorption.
An exact `(radius_bits,density_bits)` cache key prevents a nearby grid from
being substituted. Cache size affects memory and construction cost only.

The geometric coordination center opposes the sum of unit donor-bond vectors.
It is a geometric convention, not a calculated lone-pair orbital or the sum
of raw bond displacements used in some historical datasets. SDF topology is
authoritative; XYZ inference uses the existing covalent-radius cutoff. Planar
sign ambiguity requires an explicit center. The historical planar clearance
interval encloses coordinate rounding conditional on its represented normal;
it does not enclose uncertainty in that inferred normal.

Total volume averages the same three donor-plane estimates used for regional
extrema. Each region has its own finite-grid population denominator. Zero
coordinates are assigned to the positive region and near is negative local Z.
Consequently, summing regional values need not recover the total on odd grids.
This documented estimator differs from a single common voxel-volume partition.
Radius profiles preserve requested order, duplicates, fixed-center semantics
and per-radius failures; changing radius also changes the discrete grid.

Atomic hard spheres, element radius fallbacks, scale 1.17 and hydrogen exclusion
are model choices. They do not establish quantum electron exclusion, solvent
accessibility or experimental selectivity. No constant was fitted to improve
agreement. The exact-node DBSTEP comparison is deliberately separated from
continuous-volume error and from DBSTEP's additional native sphere clipping.

### Pyramidalization

Three distinct donor-to-neighbor directions are normalized in f64, then sorted
to make the reduction deterministic under neighbor permutations. P is the
absolute scalar triple product, with the retained `2-P` inverted-pyramid
correction. Alpha is the mean signed angle to the other two directions' plane
normal. The sign is determined by the apex projection on their bisector.
Dot products entering acos are clamped to its mathematical domain. Coincident
or collinear defining directions return errors; valid planar results remain
valid numbers. Analytical planar, tetrahedral, orthogonal, nearly collinear
and permutation witnesses complement the external comparisons.

The reported `pyr_alpha` is not automatically Morfeus's separate `P_angle`
quantity. Neighbor mapping must be supplied or matched: Morfeus can use its
nearest-three or connectivity convention. See the
[official definition documentation](https://digital-chemistry-laboratory.github.io/morfeus/pyramidalization.html).

### Conformer identity, statistics and populations

Ensemble input retains frame identity and order, including unsuccessful
records. Minima and maxima retain all exactly tied conformer identities.
Arithmetic mean, median and population standard deviation use the available
f32 descriptor values with f64 accumulation; the standard deviation divides
by N, not N-1. This is ensemble dispersion, not an uncertainty in its mean.
The default missing-conformer policy requires completeness. Explicit available
subsets report missing identities and retained population mass; weighted means
then normalize by that retained mass. Failure and renormalization are visible.

Explicit population APIs require finite nonnegative supplied weights and
provenance, or complete energy metadata, a common unit/kind/model, positive
temperature and stated degeneracy/state-membership assumptions. Boltzmann
weights are proportional to `g_i exp[-(E_i-E_min)/(RT)]`. Log shifting, staged
division and compensated normalization protect extreme cases. Gibbs energies
must have been evaluated at the population temperature. Electronic/force-field
energies are labelled approximations to equilibrium free-energy populations.
Conformer generation, duplicate-state detection and equilibration are caller
responsibilities; numerical normalization cannot establish them. Legacy record
aggregation conventions remain distinct from these explicit population APIs.

Energy conversions use SI exact charge and Avogadro constants, the
thermochemical 4184 J/kcal, and the stated CODATA 2022 Hartree energy. The
existing gas constant equals `N_A*k_B/4184` in the declared precision. These
values agree with the [NIST constants table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt)
checked during this campaign; they were not changed.

### Density grids

Cube ingestion validates dimensions, numeric counts, finite scalar fields,
coordinate metadata and explicitly declared total electron density. Fortran
D/d exponents retain the same conversion and failure behavior after U1.
Analysis accepts positive, axis-aligned isotropic grids only; other parsed
vectors remain metadata and are rejected for this validated analysis route.
The sphere must fit within the sampled grid. No padding, extrapolation or
interpolation is introduced.

The quantity is the node quadrature of `density > threshold` inside the
explicit sphere. Absolute volume is occupied count times the determinant-based
cell volume; percentage uses the sampled sphere count. The analytic sphere
volume and its sampling discrepancy are reported separately. Equality at the
density threshold is unoccupied. Cell spacing, threshold and origin jointly
affect discretization error. Gaussian/DBSTEP agreement validates those tests,
not the electronic structure calculation or a universal isodensity threshold.
The new hard-sphere reference API does not apply to density superlevel sets.

### Eyring rates, selectivity and temperature scans

The absolute unimolecular Eyring rate assumes unit transmission coefficient:
`k=(k_B*T/h)*exp(-DeltaG_dagger/(R*T))`, in s^-1 with kcal/mol and kelvin.
The logarithmic API avoids absolute-rate underflow; the checked f32 API reports
underflow/overflow, while the legacy convenience API retains zero/infinity or
NaN behavior. A barrier difference alone cannot supply an absolute rate.

Selectivity uses a stable logistic ratio and `100*tanh(abs(DeltaDeltaG)/(2RT))`
for absolute ee. The declared signed difference is G‡(S)-G‡(R), so positive
values favor R. Equal barriers are racemic. Temperature scans are conditional
on the supplied energy model; holding a free-energy difference fixed across
temperature is an assumption, not a prediction of entropy/enthalpy variation.
The low-level range helper is not an unbounded-resource API: astronomically
many requested steps cannot be represented as an in-memory vector.

### Regression, intervals, applicability and empirical studies

Feature construction, target units, descriptor aggregation and prediction
serialization are part of the model contract. The AVX2 eight-term predictor
and scalar fallback have different floating-point reduction orders; hardware
dispatch is not a universal bitwise-reproducibility guarantee. S1 repairs a
demonstrated dimensionful R² cutoff. It does not provide arbitrary-range
statistics when f64 variance itself underflows or overflows.

The fitting code standardizes selected features and solves normal equations;
even the nominal zero-ridge route adds the retained 1e-10 diagonal floor to
nonintercept terms. Near-singular problems therefore have a stabilizing
regularization convention, not an exact arbitrary-precision OLS solution.
Feature-scale cutoffs, pivot checks and LASSO convergence remain explicit
numerical limitations. Replacing them with a different estimator would require
new model semantics and independent validation; this campaign does not silently
retrain historical models.

Validation distinguishes fixed-feature leave-one-out, grouped folds, and
regularization tuning conditional on an already selected feature set. The
compatibility key `nested_loo` does not imply nested feature selection.
Coefficient intervals assume the stated fixed linear model and residual
distribution. Marginal coefficient envelopes do not guarantee joint coverage;
bootstrap propagation retains whole sampled vectors. Ordinary Mahalanobis
distance requires a nonsingular ordinary covariance and valid stored training
points; no ridge distance is relabelled ordinary. Range and leverage checks
are diagnostics, not proof of prospective reliability.

Historical empirical coverage failures, repeated-ligand panels, failed outer
folds, inaccessible primary supporting information, and conflicting published
ligand values remain limitations. Inclusive classifier threshold corrections
and existing provenance/geometry-refusal rules are replayed unchanged. Storage,
search, hashes and serialization were reviewed as preservation/order boundaries;
they do not create additional physical descriptors.

Finite-grid agreement is not a certificate of continuous-volume accuracy. The
fresh eight-geometry study shows nonmonotone convergence and an external near
contact that can be missed entirely. The current default remains the validated
finite-grid definition. Density is a separate scalar-field descriptor; it does
not validate the hard-sphere physical model or a universal density threshold.

Only x86_64 hardware is available. Cross-architecture bitwise equivalence is
not established by the current runs. The AVX2 model predictor and its scalar
fallback use different reductions, so identical real-arithmetic formulas alone
cannot establish identical floating-point bits.

## Engineering defect S2: declared minimum Rust version did not build

The manifest declares Rust 1.85, but the frozen baseline uses later let-chain
syntax. A fresh official Rust 1.85.0 toolchain reproduced the compiler failures;
[the original log](msrv_baseline_failure.log) is preserved. Six existing
conditionals in model validation/fitting, screening and input-path resolution
were rewritten using equivalent Option/Result predicates or nested checks.
AVX2 unsafe blocks were extended to include arithmetic intrinsics that were
unsafe APIs on that older compiler; dispatch and operation order are unchanged.

All 364 current Rust unit/integration/example tests pass on actual Rust 1.85.0,
with no test tolerance changes. The current 1.97 toolchain also passes the
suite and strict Clippy. This establishes compiler-version compatibility on
x86_64, not numerical validation on different hardware architectures. Original
failed builds and the passing replay are retained under
`.stericx/ultimate_engine/msrv_*`. The declared minimum was not raised.
