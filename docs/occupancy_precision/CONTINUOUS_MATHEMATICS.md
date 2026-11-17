# Continuous target and candidate certificate

The target is
\(V=\int_{\mathbb R^3}1_{B(c,R)}(x)\,1_{\cup_i B(c_i,r_i)}(x)\,dx\).
All input binary64 coordinates and configured scaled radii are interpreted as
exact real values. The finite-grid estimator is a separate quantity. Chemical
radius and coordinate uncertainty are outside this numerical certificate.

The existing reference encloses unknown boundary boxes by [0, box volume].
For a regular surface of area A and cell side h, unresolved volume scales like
A h while boundary-cell count scales like A/h². Thus the basic enclosure width
can decay only like the inverse square root of boundary work in that regime.
This is an asymptotic explanation, not a measured hardware or runtime law.
The fresh twelve-case replay has bounds of ±1.176–3.065 Å³ at 262,143 boxes;
merely accelerating its classification does not tighten those mathematical bounds.

## C1: chord bounds and certified smooth patches

By Fubini, integrate over the XY square [-R,R]² the length of the union of
atomic Z intervals, clipped to the integration sphere's interval. A sphere's
endpoints are \(c_z\pm\sqrt{q}\), with
\(q=r²-(x-c_x)²-(y-c_y)²\), wherever q is nonnegative.

For each rectangle, outward-rounded interval arithmetic encloses q and both
endpoints. A nonpositive upper q excludes positive chord length. If lower q
is positive, the interval formed by the upper bound of the lower endpoint and
the lower bound of the upper endpoint is contained in every chord throughout
the rectangle. The converse outer interval contains every chord. Clip these
intervals to the analogous integration-sphere bounds, merge their unions and
multiply their lengths by rectangle area. This supplies conservative cell
volume bounds even where circles or endpoint changes make the integrand nonsmooth.

For tighter bounds, require all contributing chords to exist throughout the
rectangle and prove which endpoint wins each clipping operation. Connected
components of the **inner** intervals prove connectivity everywhere. Require
the outer spans of distinct components to be disjoint. Within each component,
interval comparisons must prove a single minimum lower and maximum upper
endpoint (identical sphere endpoints are recognized). Only then is union length
a fixed sum of smooth sphere endpoint functions over the whole rectangle.

For a hemisphere endpoint, both diagonal second derivatives have magnitude at
most \(r²/q_{min}^{3/2}\). Applying the one-dimensional midpoint remainder in
each axis gives the rectangle integral error bound

\[\frac{A}{24}\left(w_x²\sup|f_{xx}|+w_y²\sup|f_{yy}|\right).\]

Sum these bounds for every selected endpoint. Evaluate the exact rectangle
midpoint through an interval enclosure, and outward-round its value, area,
derivative bounds and remainder. Intersect the resulting integral interval
with the coarse chord interval. If derivative denominators underflow or the
required conditions cannot be proved, keep the coarse bound and subdivide.

This adapts the standard derivative-remainder verification approach discussed
in [Rump, *Verification methods*, §12](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf).
The endpoint/topology argument above is specific to this implementation and
requires independent tests; citing general quadrature theory does not validate
it. Classical quadrature identities and their smoothness conditions are also
documented in [NIST DLMF §3.5](https://dlmf.nist.gov/3.5).

Refine the cell with the widest volume interval, bisecting its longest XY side.
Parent-excluded sphere candidates remain excluded from every child because
children are subsets. All aggregate additions and parent replacements are
outward-rounded. Keep the intersection of successive global enclosures
separately from raw accumulated bounds to ensure nesting without corrupting
subtraction bookkeeping. A work/representability limit returns a false target
flag and the actual remaining error bound. Relative error is bounded by absolute
error divided by the certified positive lower volume; otherwise it is undefined.

The public API passes 167 original analytical enclosures (18 strict stress
targets remain unmet), 96 additional triangle-free many-sphere/rigid-motion
observations, 84 molecular work/accuracy observations and 6,000 independent
rational predicate checks. All twelve original molecules meet 0.05 Å³ within
262,143 cells and 0.01 Å³ within 1,048,575 cells. The independent octree intervals
also overlap the new intervals; this corroborates consistency and is not an
exact molecular-volume reference. See `c1_validation.json`.


## Explicit Rust interface

```rust
use steric_x::geometry::{AtomicSphere, CertifiedVolumeOptions, sphere_union_volume_certified};
let region = AtomicSphere { center: [0.0; 3], radius: 3.5 };
let atoms = [AtomicSphere { center: [1.0, 0.0, 0.0], radius: 1.7 }];
let result = sphere_union_volume_certified(region, &atoms, CertifiedVolumeOptions {
    absolute_tolerance: 0.05, max_cells: 262143,
})?;
// Inspect result.target_met and result.status before accepting a requested tolerance.
```

`BuriedVolumeCalculator::compute_certified_with_center` uses the existing f64
continuous scaling convention and an explicit center. It returns total volume;
finite-grid regional calculations and orientation averaging retain their existing
meaning. Runtime includes preparation and integration (and atom conversion for
the molecular convenience call). Geometric evaluations count sphere/chord
queries and preprocessing pair queries, not arithmetic instructions. Cell budgets
count all classified XY rectangles; the original Reference budget counts 3D
boxes, so equal counts are not equal costs. Analytical shortcuts use no XY cells.

The certificate assumes IEEE round-to-nearest, gradual underflow and no unsafe
fast-math/FMA rewrite. Coordinates/radii are exact supplied binary64 values.
Unsupported/nonfinite inputs are rejected; budget exhaustion and arithmetic
resolution limits return false target flags. Contradictory independently formed
bounds return an error rather than a fabricated interval. Relative error uses
the positive certified lower volume; it is undefined when that lower bound is
zero. Analytical and adaptive bounds include arithmetic rounding. Runtime is
excluded from deterministic scientific-output comparisons.


## Exact positive-interval specialization

For finite nonnegative endpoint intervals, multiplication's minimum is the
product of the lower endpoints and its maximum is the product of the upper
endpoints. Round-to-nearest multiplication and adjacent-float widening are
monotone, so selecting these two products yields the same bound bits as the
general four-product/min/max implementation. Exact-zero handling is retained.
Signed or infinite intervals use the original general implementation, including
its NaN-product fallback. Positive division uses the identical outward reciprocal
followed by this multiplication.

A deterministic 100,000-case test spans binary64 exponents, underflow, overflow,
signed zero and fallback intervals against the unchanged four-product reference.
All 24 molecular tolerance cases preserve every certificate and work decision in
eight alternating experimental pairs. This changes the cost of certified
arithmetic, not the certificate or its acceptance threshold.


## Per-integration scratch reuse

Span, range, ordering, component and endpoint vectors are cleared and reused
within one integration. A component stores an index range into the same sorted
order rather than a newly allocated index vector. The endpoint comparisons,
iteration order, floating-point operations, heap priorities and counters are
unchanged. Every retained cell owns its candidate indices; no heap cell refers
to reusable scratch storage. No scratch or geometry survives an API call.

Paired research trials preserve all bounds, work counts and stopping decisions
on the original twelve molecules at 0.05 and 0.01. Full optimized-certificate
replay is an additional admission gate, not replaced by these performance trials.
