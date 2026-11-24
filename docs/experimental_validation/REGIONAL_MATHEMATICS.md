# Certified continuous regional buried volume

## Quantity and explicit interfaces

For the integration ball S=B(c,R) and atomic union A, let Q0=(x>=cx,y>=cy),
Q1=(x<cx,y>=cy), Q2=(x<cx,y<cy), Q3=(x>=cx,y<cy). Octants 0–3 also have
z>=cz; octants 4–7 have z<cz. Vj=∫1S 1A 1Qj dx (or the corresponding octant).
Boundary planes have zero Lebesgue measure; their inclusive side is consistent
with the existing grid labels. Regional percentage is 100 Vj/(4πR³/(3k)),
k=4 or 8, and total percentage uses k=1. These are continuous quantities,
separate from independently normalized finite lattice ratios.

`sphere_union_volume_certified_regions` accepts exact represented f64 balls and
returns total, four quadrants and eight octants, each with estimate, lower/upper
volume bounds, absolute/relative error, percent enclosure, target flag and status.
The enclosing result reports the shared computational budget, all classified cells,
subdivisions, geometric evaluations, smooth/range/empty pieces, heap peak and runtime.
A work item contributes to both Z halves; no fictitious independent per-bin runtime
is assigned. The shared runtime/counters apply to every member integral.

`BuriedVolumeCalculator::prepare_certified_frame` uses the existing donor/plane
basis and f32 coordinate subtraction/dot convention, exposes axes and local balls,
and uses the same f64 radius scaling as the explicit continuous APIs. Its
`compute_regions` method certifies these **represented local balls**. It does not
certify the error of the molecular frame rotation, source f32 rounding, empirical
radii, or chemistry. The caller supplies center and plane; repeat for the three
existing donor-plane choices when desired. Default finite-grid routines are unchanged.
The molecular default center distance is 2.1 Å; the Kraken-matched experimental
study separately fixes 2.28 Å and records that choice.

## Shared subdivision and plane clipping

Four disjoint quadrant rectangles partition [-R,R]² at x=y=0. Each adaptive cell
stores a quadrant ID and bounds for both positive and negative Z volume. Geometry
and candidate queries are performed once. First intersect each atomic chord with
the integration-ball chord. Then intersect with [0,+∞) and (-∞,0] in local Z.
The plane endpoint is the exact constant zero, whose first/second derivatives vanish.
Interval min/max enclose clipping. If interval comparisons cannot prove which
endpoint wins, the endpoint identity is unknown and smooth integration is disabled.

As in the preceding certified chord method, union the inner chords to bound the
minimum occupied length and outer chords to bound the maximum. Multiply by an
outward area enclosure. Missing chords contribute zero lower length. All sphere,
plane tangency, emerging chord and ambiguous overlapping-boundary cases retain this
coarse enclosure, without assuming smoothness at a singularity.

A smooth bound is admitted only when every contributing chord exists over the
rectangle, inner spans prove connectivity, outer spans prove separation of
components, and interval ordering proves each active minimum/maximum endpoint.
Endpoint functions are hemispheres cz±sqrt(r²-(x-cx)²-(y-cy)²), or constant zero.
For each curved endpoint, |fxx|,|fyy|<=r²/q_min^(3/2). The tensor midpoint remainder
is A(wx² sup|fxx|+wy² sup|fyy|)/24, summed over selected endpoints. Constant plane
endpoints contribute zero. All midpoint, area, curvature and remainder operations
are interval-enclosed. Intersect this bound with the coarse bound; contradictory
bounds return an error. Underflow/nonfinite curvature falls back to range bounds.

Refine the cell with greatest summed Z-half uncertainty and bisect its longest XY
side. A child inherits only parent candidates: geometric exclusion on a superset
proves exclusion on its subsets. The binary cells partition the original domain;
no point sampling or pairwise-overlap approximation defines the certificate.

## Global error guarantee and consistency

For each octant, downward lower sums and upward upper sums enclose its integral.
Replacing a parent uses outward subtraction/addition and retains rounding slack.
Maintain an independent intersection of successive enclosures to guarantee nested
bounds as the cell budget grows. Restrict volume lower bounds to zero and octant
upper bounds to the known sphere/8 volume. Quadrant intervals are outward sums of
their two octants. Total bounds are outward sums of the four quadrants, intersected
with [0,V(S)]. Thus every reported global/regional integral contains the real
quantity if the local predicates and IEEE arithmetic contract hold. Estimates are
midpoints; the absolute error is the maximum endpoint distance, enclosed using
TwoSum. Relative error divides by a strictly positive certified lower volume,
otherwise it is undefined. Percentage division uses a positive interval denominator
and intersects with the physical [0,100] domain; if its lower denominator is zero,
percentage is unavailable rather than fabricated.

Stopping requires the total **and all twelve regions** to meet the original
requested absolute tolerance. The one global budget includes four roots and hence
must be >=4. Exhaustion returns current certificates and false flags for unmet
members. Unsplittable cells retain their uncertainty; an empty refinable heap
returns floating_point_limit. Different members may individually meet tolerance
when the shared calculation fails. Independently calculated estimates need not
sum bitwise; containment of their common mathematical quantities is the requirement.

## Independent evidence

`validate_regions.py` builds exact Fraction quadratic cross-section envelopes,
splits them at the regional plane, and integrates each polynomial exactly. The
only transcendental factor is an independently bounded Machin-formula π. This
reference supports genuine triple/higher overlaps and never calls the production
chord integrator. Exact sphere caps and eight separated balls supplement it.
Percentage references cancel π algebraically. Axis permutations and exact translated
centers exercise all region labels and local-coordinate subtraction.

The first admitted mathematical run contains **96 geometries / 1,248 regional
integrals, zero invalid volume or percentage enclosures**. All 92 requests with
0.01/0.02 Å³ targets succeed; four deliberately strict 1e-8 requests at
4/32/1024/16384 cells remain unmet and visible. All per-region flags and shared
cell counts pass; exact rational sums verify quadrant/octant consistency. Five
Rust integration tests additionally cover nesting, duplicate/permutation invariance,
full/empty volume, explicit input failures and reuse of existing molecular axes.

Exploratory runtime median is 18.82 ms, maximum 195.97 ms; maximum classified cells
194,858. Whole observer max RSS is 22,356 KiB. These shared-workstation observations
are diagnostics, not controlled performance claims. Further molecular, radius and
full regression gates are separate evidence.

An initial regional test exposed a negative-subnormal lower bound for an empty
quadrant from unconditional interval addition. Although it was a conservative
mathematical enclosure, it violated the nonnegative volume contract. Intersecting
regional sums with [0,+∞) fixes that contract without weakening any tolerance.
The failed source and test log remain in `regional_v1/`; passing evidence is in
`regional_v2/` and `regional_independent_v1/` under the campaign evidence directory.
