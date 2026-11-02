# Independent mathematical reference for buried-volume resolution

This document specifies the research reference and its acceptance assumptions
before it is used to judge a changed comparison policy. The helper is
`scripts/dbstep_resolution_reference.py`; it imports neither StericX nor DBSTEP
nor Morfeus. It is not a second production chemistry engine. The production
descriptor implementation remains Rust.

## Mathematical quantity

Let S be the closed evaluation ball with center c and radius R, and let U be the
union of the selected positive-radius atomic balls. The continuous quantity is

```text
V = volume(S intersect U)
percentage = 100 V / ((4/3) pi R^3).
```

Atomic overlaps count once. Zero-radius entries are explicitly treated as ghost
atoms with no occupied nodes. This convention agrees with their zero continuum
volume but is separately stated because a finite point exactly at a zero-radius
center would otherwise have nonzero grid weight. The reference uses supplied
coordinates, radii, scaling and center; it does not infer bonds, centers or radii.

The helper interprets each supplied binary float as its exact rational value.
This separates the calculation from unavailable pre-rounding coordinates and
from physical uncertainty in an atomic-radius model. No tool is physical ground
truth merely because its implementation uses a wider float type.

## Exact finite-point witnesses

For point p and ball (a,r), the signed squared-distance margin is

```text
m = sum_j (p_j - a_j)^2 - r^2.
```

`sphere_margin_exact` evaluates m as a `Fraction`, using the exact binary
representations. Inclusion in the closed ball means m <= 0. This supplies an
independent answer even where two float operation orders disagree. For example,
the binary64 values written as 0.6 and 0.8 do not form an exact 3–4–5 triangle:
their exact squared-distance margin relative to radius 1 is positive.

`union_mask_exact` combines these predicates by set union for small witnesses.
`union_mask_f64` is a transparent, direct vectorized calculation for larger point
sets. Its result is explicitly a binary64 computation, not an exact predicate;
all disputed or roundoff-sensitive decisions require the rational witness.
Agreement of total counts alone cannot conceal compensating point differences.

`f32_predicate_certificate` supplies a separate forward-error bound without
evaluating a binary32 copy of the production kernel. For actual accepted
binary32 inputs let d_j=p_j-a_j, s=sum(d_j^2), m=s-r^2 and u=2^-24. A rounded
subtraction contributes one relative factor; squaring repeats that factor and
adds another; at most two positive additions give at most five factors per
term. Consequently `abs(s_computed-s) <= gamma_5*s`, with
`gamma_5=5u/(1-5u)`. Radius squaring contributes at most u*r^2. A changed point
decision therefore requires `abs(m) <= gamma_5*s + u*r^2`. Points strictly
outside that band have a stable inclusion/exclusion decision; any disagreement
there fails. Points inside remain ambiguous and retain their exact rational
margins and actual operation-order traces; ambiguity alone is not a pass.
The band is derived from the operation graph, not from a measured discrepancy.

The certificate checks the sufficient finite-normal domain explicitly using
exact rational inputs. Each nonzero d_j^2 must be at least
`N/(1-u)^3`, `s*(1+gamma_5)` must not exceed F, and r^2 must lie between
`N/(1-u)` and `F/(1+u)`, where `N=2^-126` and `F=2^128-2^104`. Zero coordinate
differences are exact and allowed. These checks exclude intermediate underflow
and overflow; unsupported cases cannot silently receive a certificate. Inputs
must already be exactly representable in binary32. In particular, the common
radius is the explicitly declared scaled radius; a different scaling or radius
table remains a distinct input convention. The sphere test without subtraction
would permit gamma_3; using gamma_5 is a stated conservative bound, not a fit.

An exact finite-set comparison must specify the actual set of sample points.
If a reference clips a supplied already-selected set using a different sphere
predicate, it has changed that finite set. The resulting denominator must be
recorded separately. A comparison with identical atomic occupancy but different
sphere membership is not an arithmetic-only comparison of the same estimator.

## Analytic sphere intersections

For ball radii R and r with separation d:

* d >= R+r gives zero intersection volume.
* d <= |R-r| gives the complete smaller-ball volume.
* Otherwise, joining the two spherical caps gives

```text
V = pi (R+r-d)^2 [d^2 + 2d(R+r) - 3(R-r)^2] / (12d).
```

The coefficient of pi is calculated exactly as a rational number, avoiding
catastrophic subtraction in nearly tangent witnesses. A certified pi interval
comes from Machin's identity `pi = 16 atan(1/5) - 4 atan(1/239)` and the two
successive partial sums of each alternating arctangent series. The default
64 terms enclose pi more tightly than 1e-80. Final volume bounds are rounded
outward to adjacent binary64 values. No fitted numerical tolerance enters the
analytic interval.

This formula also validates two overlapping atomic spheres wholly inside S:
their union volume is the sum of individual volumes minus their pairwise lens.
Axis-aligned test separations are supplied exactly; the analytic helper does not
silently certify a separately rounded square-root calculation of d.

## Independent continuous-volume enclosure

The adaptive reference partitions the cube enclosing S into dyadic boxes. This
is a different integration method from center-point counting on either tool's
lattice. For box [l,h] and ball center a, the exact squared-distance extrema are

```text
d_min^2 = sum_j max(l_j-a_j, a_j-h_j, 0)^2
d_max^2 = sum_j max(abs(l_j-a_j), abs(h_j-a_j))^2.
```

A box is outside a ball if d_min^2 > r^2, and wholly inside if d_max^2 <= r^2.
It is outside S intersect U if it is outside S or every atomic ball. It is
wholly inside the intersection if it is inside S and at least one atomic ball.
Joint covering by several balls without one covering the entire box remains
unresolved; this is conservative rather than an incorrect rejection.

The exact dyadic endpoints are enclosed in binary64 intervals. Every subtraction,
square and addition in the min/max predicate is expanded outward by one adjacent
binary64 number, including the radius square. An inconclusive interval remains
unclassified. Independent rational box predicates check the direction of every
definite classification in the randomized/adversarial unit tests.

Definite interior boxes contribute to both lower and upper volume bounds.
Unresolved boxes contribute only to the upper bound and are subdivided at the
next depth. Volumes accumulate as exact rational multiples of the root cube's
volume; only the final reported endpoints are converted outward to binary64.
This avoids accumulating a large uncertain floating sum of voxel volumes.

The method produces a certified interval conditional on the supplied geometry
and standard IEEE round-to-nearest binary64 arithmetic. Each depth records the
number of examined, inside, outside and unresolved boxes. Nested bounds should
increase/decrease monotonically, respectively. The interval width, runtime and
maximum depth must be published. A wide interval is not evidence of high
precision; no observed tool residual is used to select the width or stopping rule.

## Discretization and ratio normalization are separate

For an exact regular center lattice with cubic cell side h, let delta be half
the cell diagonal, `sqrt(3) h / 2`. Assigning every cell according to its center
can disagree with the continuous indicator only within delta of a boundary.
For a ball radius r, a conservative full boundary-layer volume is

```text
B(r,delta) = (4/3) pi [(r+delta)^3 - max(r-delta,0)^3].
```

Thus a conservative occupied grid-volume error bound is B(R,delta) plus the sum
of B(r_i,delta) over relevant atomic balls, and the sampled sphere-volume error
is at most B(R,delta). This union bound can be very loose; overlapping atomic
surfaces are counted repeatedly. It must not be advertised as an observed
accuracy or used alone to establish a precise default resolution. The adaptive
box interval supplies a more informative independent enclosure.

The boundary-layer argument assumes a complete cell tiling and exact center
predicates. Rounded/nonuniform lattice positions, finite endpoints and ambiguous
float predicates require explicit extra accounting; the ideal-lattice formula
cannot simply be applied to arbitrary captured coordinates.

Absolute grid-count volume and an analytic-sphere-normalized fraction are
different finite estimators. Write Qocc=Nocc*h^3, Qsphere=Nsphere*h^3, and let
Eocc and Esphere independently bound their errors from continuous Vocc and
Vsphere. For Qsphere>0, exact algebra gives

```text
Vsphere Qocc/Qsphere - Vocc
  = (Qocc - Vocc) + (Vsphere - Qsphere) Qocc/Qsphere

absolute error <= Eocc + (Qocc/Qsphere) Esphere.
```

`ratio_volume_error_bound` evaluates that bound with rational arithmetic and
outward conversion. Both denominator and numerator must have valid, separately
established error bounds. A denominator formed using a different predicate is
not silently replaced with the numerator's sphere count. Clamping a computed
percentage to [0,100] is not justified by this bound and can hide an inconsistent
finite estimator; raw counts and the unclamped value remain required evidence.

## Prospective validation lanes

1. **Finite predicate and normalization:** independently compare exact supplied
   geometries, radii, every relevant point decision, denominator, and algebraic
   normalization. Record exact margins for all disagreements; distinguish
   documented rounding uncertainty from an incompatible call/comparison contract.
2. **Continuous quantity:** analytic cases and adaptive certified intervals test
   the underlying volume definition independently of the sampled lattice.
   Report each estimator's distance to the certified interval and possible
   continuous error range; interval overlap is not exact descriptor agreement.
3. **Resolution and cost:** use a predeclared refinement ladder and cohort,
   retaining every result, failure, elapsed time and nonmonotonic step. Compare
   native density controls with DBSTEP spacing controls only after recording
   their actual point counts and lattice positions. Last-step differences are
   sensitivity measurements, not certified error bars.

This methodology does not itself approve a replacement release gate. The complete
comparison policy, including what constitutes an explained rounding difference
and what quantitative accuracy is claimed, must be documented and approved before
the changed comparison is scored. The original strict failed gate remains part
of the historical record.

## Reference-development check

The first helper test run had 12 passes and one unmet enclosure-width request:
the R=1, r=0.75, d=0.5 analytic lens was correctly enclosed at depth 6, but its
interval width was 0.292724609375 rather than below the test's 0.25 volume-unit
resolution target. The formula and target were retained and the test was refined
to depth 7. This was a request for more reference resolution, not a failed
analytic inclusion or a widened scientific acceptance tolerance. Subsequent
focused tests additionally exercise the forward-error certificate's valid and
rejected domains.
