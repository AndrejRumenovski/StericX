# Numerical accuracy investigation

Status: baseline convergence measured; new algorithms below are proposed and
are not yet admitted or exposed as validated production modes.

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
