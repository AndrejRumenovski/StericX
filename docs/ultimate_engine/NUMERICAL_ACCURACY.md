# Numerical accuracy investigation

Status: baseline convergence and explicit N1/N2 methods validated in the
scopes below. Final whole-engine replay and controlled performance audit remain
pending. The method proposals are retained to explain their admission criteria.

## Quantities and baseline

Hard-sphere volume is the Lebesgue volume of S ∩ (union of atomic balls), where
S is the specified integration sphere. Existing native results instead report
a finite-lattice occupied fraction times 4πR³/3, averaged over three molecular
frames. Region denominators are separate. A lattice's agreement with another
program using precisely the same points cannot bound continuous-volume error.
The eight analytical cases, raw runs and rational/π reference enclosures are
recorded in [baseline_measurements.json](baseline_measurements.json).

Increasing lattice density does not monotonically improve error: the overlapping
pair has absolute errors 0.827624, 0.0532848, 0.0803326 and 0.0148628 Å³ over the
four tested spacings. A near-tangent lens of volume approximately 3.31e-10 Å³
receives a discrete nonzero value at one resolution and zero at a finer one.
Any precise/reference API must expose its method and limitations separately.

## Continuous Sterimol support envelope: planned finite-event method

After an orthonormal axis projection, each atomic ball has transverse disc
center p_i and radius r_i. Define h(n)=max(0,max_i(p_i·n+r_i)) for unit transverse
n. Continuous B1 is min h(n); B5 is max_i(||p_i||+r_i), and raw L is
max_i(z_i+r_i), without the optional historical 0.40 Å correction.

On any open angular segment with a unique active disc, a local minimum occurs
at a stationary point of that disc's support. Otherwise a minimum occurs at
an equality of two supports. Enumerate the negative-center stationary
directions and both unit solutions of (p_i-p_j)·n=r_j-r_i. Include a zero-radius
disc at the origin to represent the outer max with zero; coincident centers
need no pair event, and a default direction covers constant envelopes. Evaluate
the entire support envelope at every candidate and keep the minimum.

This finite candidate set contains a global minimizer in exact real arithmetic.
A straightforward implementation costs O(N³) time and O(N) storage. It removes
angular-grid bias, but floating-point projection, hypot, square roots and
near-tangent equality decisions still need numerical validation. It is not a
certified floating-point error bound. Proposed implementation uses f64 throughout
and returns f64, while retaining the original f32 input precision. Validate with
analytical disc arrangements and an independent high-precision angular search
with a Lipschitz enclosure; measure cost against dense angular scans.

## Continuous volume: planned independent approaches

Precise candidate: Fubini reduction. At each XY integration location, clip every
atomic sphere's exact Z chord to the integration sphere's Z chord, merge the
intervals, and integrate their union length over XY. This reduces a 3-D
occupancy grid to 2-D quadrature. Refinement differences are empirical error
estimates, not certificates, because circle boundaries and interval topology
produce nonsmooth integrands. Analytical single-sphere intersections supply
independent references and potential explicitly identified hybrid cases.

Reference candidate: adaptive box enclosure. Bound squared distance from each
axis-aligned box to each sphere using outward-rounded interval operations.
Discard boxes certainly outside S or every atomic sphere. A box certainly
inside S and at least one atomic sphere contributes to both lower and upper
bounds; unresolved boxes contribute only to the upper bound. Bisect unresolved
boxes and accumulate conservatively. The union of disjoint leaves covers the
initial integration-sphere bounding box, so these are geometric volume bounds
if every arithmetic step is outward-rounded. A work limit must return unresolved
bounds and a false convergence flag, never an invented guarantee. Test exact
dyadic classification with an independent rational oracle and analytical volumes.

Higher precision and interval methods address arithmetic and integration error,
not uncertain molecular conformations, radii or the hard-sphere physical model.
Changing radii to match an external result is excluded. Existing default calls
remain the fast finite-grid method; new methods require explicit selection.

## N1 validation and admission

The continuous-event Sterimol method is implemented as explicit
`SterimolCalculator::compute_continuous` and `compute_continuous_with_dummy`
methods returning `ContinuousSterimolParams` with f64 fields. Existing entry
points and their 360-direction f32 arithmetic are unchanged. The new methods
also retain f64 axial intermediates through the L reduction. They do not invent
precision absent from the f32 source coordinates/radii.

Four analytical/geometry/validation unit tests pass. The independent angular
interval search uses 80 decimal digits and a per-disc lower bound on each
angular interval; it does not enumerate the Rust stationary/pair events. It
converged to width <=1e-12 times geometric scale on 39 cases: 27 analytical,
random and adversarial arrangements, and 12 frozen molecular geometries. Every
continuous result lies within the reference bracket plus a separately declared
128-epsilon arithmetic allowance. The numerical bracket is a high-precision
validation reference, not a formal mpmath rounding certificate. L and B5 also
agree with independent 80-digit projections/formulas.

[Full accuracy/runtime samples](n1_accuracy_measurements.json) include five
native timing observations per method per case. For the 32-disc random witness,
the event method uses about 49 microseconds versus about 789 microseconds for a
23,041-direction scan; the latter retains about 1.57e-4 Å of B1 error, whereas
the event value lies in the ~2e-12 Å reference bracket. The default scan takes
about 11 microseconds but has about 7.87e-3 Å error. This is an accuracy/cost
comparison, not a claim that equal-accuracy legacy calls are 16× faster.
Numbers are exploratory and final controlled measurement remains pending.

### External convention investigation and preserved failures

The first external check incorrectly assumed Morfeus's B5 was the analytic
radial maximum. Its inspected `calculate` method actually maximizes the finite
angular support grid. Preserve that failed comparison at
`.stericx/ultimate_engine/n1_external`. For N endpoint-inclusive directions,
continuous B1 <= sampled B1 <= continuous B1 + rho*pi/(N-1); the analogous B5
shortfall is at most rho*(1-cos(pi/(N-1))). The revised check uses these derived
discretization bounds and the original arithmetic allowance.

A second lane exposed DBSTEP rotation's explicit rounding to eight decimal
places. Preserve `.stericx/ultimate_engine/n1_external_corrected`; its 12 real
molecule B5 differences are not defects in the continuous Rust formula. The
matched lane calls DBSTEP's unchanged descriptor kernel on an unrounded
orthonormal projection and retains the upstream rounded-alignment values in
separate fields. DBSTEP's absolute-axis L is recorded, not conflated with signed
forward L. All 39 matched comparisons pass; the float allowance was not enlarged.
See [external comparisons](n1_external_comparisons.json). Original/revised
controllers and inspected upstream source are retained with the raw evidence.

Memcheck on single-disc, 32-disc and molecular cases reports zero errors and
zero lost bytes (544 runtime bytes reachable). Repeated exact-order/invariance
unit checks, complete final replay and final controlled cost measurements remain
part of campaign acceptance. No globally certified floating-point error claim
is made for N1.

## N2: explicit continuous-volume methods

The Rust core exposes `sphere_union_volume` for f64 balls and
`BuriedVolumeCalculator::compute_continuous_with_center` for molecules. The
latter uses an explicit center, radius, radii scale and hydrogen policy. Its
bounds refer to the rounded f64 scaled radii. Both return total volume in Å³;
no new regional, density or thermodynamic interpretation is introduced.
Existing default APIs continue to use the exact native f32 finite grid.

`ContinuousVolumeMethod::MidpointSlices { subdivisions }` merges clipped Z
chords at each XY midpoint, then integrates their union length with compensated
summation. It reports no error bound or target-met flag. A fixed resolution is
not itself an accuracy guarantee. For the centered-sphere case, the measured
absolute error improves to about 4.65e-6 Å³ at 128 subdivisions, then worsens to
about 2.50e-4 at 256. Richardson extrapolation or a monotone-resolution stopping
rule is therefore not admitted for these nonsmooth integrands.

`ContinuousVolumeMethod::Reference { absolute_tolerance, max_boxes }` uses
outward-rounded analytical ball/lens formulas when applicable, including
inclusion–exclusion for two balls wholly inside the integration sphere. Other
unions use conservative adaptive boxes. Returned lower/upper bounds and
`absolute_error_bound_angstrom3` refer to the continuous mathematical volume.
`target_met` is explicitly false when the requested tolerance is not achieved.
A work limit never becomes a claimed error bound.

### Bound justification and numerical contract

Inputs are exact supplied binary64 values, finite with coordinates/radii at
most 1e100 in magnitude; atomic radii may be zero, integration radius must be
positive. IEEE round-to-nearest and gradual underflow are required. Each
arithmetic operation used in an enclosure widens toward adjacent floats.
Adjacent-float operations are implemented by IEEE bit stepping to support Rust
1.85. Squaring handles zero-crossing intervals, square roots are widened, and
uncertain products return conservative unbounded intermediates before physical
volume clipping. Zero-radius balls have zero Lebesgue volume.

The π interval is enclosed by neighboring binary64 floats around π. The
independent Python oracle checks this against rational alternating-series
bounds for Machin's formula, π=16 atan(1/5)-4 atan(1/239), using 64 terms.
The exact rational certificate is retained in
`.stericx/ultimate_engine/n2_validation_v4/pi_certificate.json`.

For a single sphere intersection, center distance is interval-enclosed.
Intersection volume decreases with separation, so endpoint evaluations give
valid bounds. Mathematical containment/disjoint/partial-overlap branches are
selected conservatively; ambiguous branches fall back to boxes. Partial lens
volume is evaluated as π(R+r-d)²[d+2(R+r)-3(R-r)²/d]/12, with every operation
outward-rounded. Two-contained-ball unions use V1+V2-Vintersection with interval
subtraction. No general many-ball union is approximated as a sum of pair terms.

For general unions, the initial box is [-R,R]³ around the integration center.
Relative atomic centers are intervals enclosing exact coordinate subtraction.
Minimum and maximum squared distance bound each box's relation to a ball.
A box is discarded only if outside the integration ball or every atomic ball;
it contributes to the lower bound only if wholly inside the integration ball
and at least one atomic ball. All unresolved boxes contribute to the upper
bound. This sufficient condition remains conservative when several balls
jointly cover a box.

Binary splits partition the real box at a representable interior point. Each
cell's volume is interval-enclosed. The lower sum accumulates confirmed lower
volumes downward. Refinement subtracts the stored parent upper volume and adds
retained child upper volumes upward; previous rounding slack is retained.
Unsplittable cells keep their upper contribution. The reported estimate lies
inside the final interval, and its maximum endpoint distance is rounded upward.
These facts establish the reported error bound under the numerical contract;
they do not establish a bound for uncertain source structures or radii.

### Independent checks and measured limitations

[Accuracy and cost data](n2_accuracy_measurements.json) contain 244 observations:
77 midpoint integrations and 167 reference integrations. All 167 reference
intervals contain independent rational/analytical volume enclosures. Twenty-two
atom-permutation checks are exact. Cases include empty/full volumes, tangent
and near-tangent lenses, overlapping unions, exact rigid transformations,
100 randomized lenses, and scales from 1e-200 to 1e99. Seven Rust tests include
10,000 independent integer/dyadic checks of box-distance and box-volume bounds.

Twenty-three requested reference targets remain unmet: 15 general-box cases
and eight extreme-underflow cases. All are reported as false. The 12 molecular
cases also fail to attain ±0.05 Å³ within 262,143 classified boxes. Their bounds
overlap the independently implemented prior octree oracle, which uses rational
cell-volume counting and different subdivision/accumulation. See
[molecular comparisons](n2_molecular_comparisons.json). This corroborates the
implementation; it is not an exact molecular-volume oracle or proof that the
midpoint estimates have a small error. These molecular bounds remain broad.

The analytical cases typically enclose their volumes within several binary64
ULPs. The general box method is expensive for tight tolerances and is retained
as an explicit conservative reference, not substituted into ordinary workflows.
Single-call nanosecond timings for analytical cases are exploratory and include
clock overhead; final cost measurements must batch these short calls.

Native and allocation-instrumented outputs match exactly for all 236 runs of
the earlier, non-underflow-expanded matrix. Memcheck on the analytical lens,
midpoint integration and general-box path finds no errors or lost bytes (544
runtime bytes reachable). [Allocation measurements](n2_allocations.json) include
zero net live-byte change. A 65,535-box overlapping-chain run requests about
4.72 MB across heap growth and peaks around 2.42 MB including the loaded input
matrix; the analytical single-ball path uses one 128-byte allocation in this
adapter. This is allocation traffic, not retained memory per calculation.

An initial stress fixture exceeded the declared coordinate domain by placing
an atom at 2e100. The engine correctly rejected it; the failed run is preserved
at `.stericx/ultimate_engine/n2_validation`. Subsequent fixtures stay inside the
published domain and use exact rational inclusion–exclusion in the reference.
No production tolerance or domain was expanded to make that run pass.
