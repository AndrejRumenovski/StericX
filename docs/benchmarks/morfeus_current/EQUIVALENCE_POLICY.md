# Prospective scientific equivalence policy, revision 1

Policy scope: the frozen implementations and 56 XYZ geometries in `freeze.json`
and `input_manifest.json`. This document is frozen before applying its rules to a
new comparison. Its SHA-256 and timestamp are recorded in
`revision_policy/policy_freeze.json`. Amendments require a new version, explicit
reason, a new hash, and a new evaluation; a failed case cannot change this policy.

## Provenance and prior exposure

No approved numerical pass/fail tolerance was found in the independent audit.
Its observed maxima are descriptive, not acceptance limits. The earlier run
already exposed all 56 residuals and timings. That history cannot be undone:
this revision is **not blinded**, and the earlier gate and performance claims are
withdrawn. The complete prior publication is archived in
`revision_policy/prior_publication.tar.gz`. A new evaluation is prospective with
respect to this frozen policy, not with respect to first exposure to this corpus.
No earlier timing is relabelled as having followed this policy.

The rules below follow definitions, exact arithmetic and certified interval
predicates, not the candidate residual distribution. There is no fitted absolute
or relative tolerance, and no use of an audit maximum. Correlation, regression,
small MAE and coincident rounded results do not independently establish admission.

## Descriptor decisions

| Descriptor | Expected relationship | Numerical source of disagreement | Independent reference | Acceptance criterion | Justification |
| ---------- | --------------------- | -------------------------------- | --------------------- | -------------------- | ------------- |
| Sterimol L | Same continuous axial support plus 0.4 Å, under matched axis and radii; cross-tool bit identity is not expected | Decimal to f32 vs f64 input, inferred center/axis, dot products, output rounding | Existing analytic support diagnostic is useful but lacks a certified propagated input/axis error bound | **Unresolved; excluded** until an independently derived interval/error bound includes input and axis construction | Observed microångström residuals are not tolerances; API couples L to the other Sterimol outputs |
| Sterimol B1 | Same continuous minimum of transverse support functions; different discrete angular estimators | Native 360 unique angles in shortest-arc +Z frame; reference 3600 linspace directions including duplicate endpoint in Kabsch +X frame; phase and precision | Existing independent analytic stationary/crossover solution of the support envelope; useful diagnostic, not a certified full error budget | **Unresolved; excluded**; no grid-error tolerance adopted | Angular-count matching does not match phase. A future bound can use the support function Lipschitz constant and angular covering radius, but must also bound axis/roundoff error independently |
| Sterimol B5 | Native analytic maximum radial support vs reference sampled angular support | Finite angular underestimation plus input/axis/radius rounding | Analytic radial maximum diagnostic, not a certified propagated error budget | **Unresolved; excluded** | Similar outputs do not make the distinct finite algorithms identical; future angular bound needs independently bounded geometry error |
| Buried volume | Same explicitly defined finite lattice estimator, not a certified continuum integral | f32/f64 grid coordinates, frame and atom radii; predicate decisions; different floating reductions | Exact rational ideal lattice; outward-rounded interval geometry from exact decimal XYZ and constants; exact rational IEEE rounding model for reductions | Candidate only: **every lattice membership and occupancy decision must be certified and identical point by point**, and both tools' scalar outputs must match their independently evaluated IEEE operation sequences bit for bit; any undecidable or differing predicate fails admission | Exact predicate equality proves the same finite scientific calculation. Exact rational rounding accounts for precision differences without choosing a tolerance |
| %Vbur | Same finite occupied fraction, averaged over three frames | Reduction order and sphere constant cancellation | Same certified counts and rational IEEE model | Same mandatory volume gate, including exact reconstruction of each tool's declared percentage expression | Mathematical fraction is shared; direct cross-tool bit equality is not expected |
| Quadrant minima/maxima and maximum adjacent difference | Same finite regional estimator over all three donor-plane grids | f32 direct quadrant ratio vs f64 sum of octant volumes; subtraction | Certified per-point region labels and counts plus rational IEEE model | Same volume gate; exact reconstruction of all plane and aggregate outputs | Matching totals alone cannot establish regional or adjacent-difference equivalence |
| Octant minima/maxima, near/far volume | Same finite regional estimator over all three frames | Normalization, sequential sum and mean rounding | Same certified per-point reference and rational IEEE model | Same volume gate; exact reconstruction of all plane and aggregate outputs | Prevents compensating regional or voxel differences from being hidden by aggregate agreement |
| Pyramidalization P | Same continuous trivalent-donor definition after matching mapping | f32 input vs f64, unit vectors, determinants/dot products and signs; independent algebra | No certified high-precision geometry/conditioning reference established in this revision | **Unresolved; excluded** | Rounding morfeus to f32 and observing equal bits is not a derived bound or a proof of correct rounding |
| Pyramidalization alpha | Same continuous angular definition | Above plus inverse trig, conditioning and degree conversion | No certified inverse-trigonometric/input propagation reference established in this revision | **Unresolved; excluded** | Numerical identity is not guaranteed for different algebra and precision |
| Other/combined descriptors | No admitted combined scope | Includes unresolved constituents or unmatched work | None | **Excluded** | Passing a subset never admits the full descriptor CLI |

## Exact finite-grid reference for the candidate volume gate

1. Read the original decimal coordinates as exact rationals; identify the sole P.
   Independently infer all three neighbors with the declared Cordero table
   (H .31, C .76, N .71, O .66, P 1.07, S 1.05 Å) and factor 1.3.
   Certify bond predicates with outward intervals. Require matching atom mapping.
2. Construct the negative sum of normalized donor-to-neighbor vectors and put the
   center exactly 2.28 Å along its normalized direction. This revision admits
   only nondegenerate, nonplanar geometry; uncertain/degenerate branches fail
   closed. For each neighbor construct +Z from donor to center, projected +X
   toward that neighbor, and +Y = Z cross X. Normalizing exact Y is immaterial
   because exact X and Z are orthonormal. All arithmetic is independently
   enclosed, including input conversion, square roots, products and sums.
3. Use exact radii H 1.20, C 1.70, N 1.55, O 1.52, P/S 1.80 Å, multiplied by
   exact 1.17. Exclude H from occupancy but retain it for neighbor inference.
4. Radius is exactly 7/2 Å, density exactly 1/100 Å³. The ideal side count is
   round(cuberoot(8 R³/density)) = 32, proved by comparing with 31.5³ and 32.5³.
   Index each lattice point by (i,j,k) in 0..31. Its exact coordinates are
   `7*(2*i-31)/62`, likewise j,k. Sphere inclusion is the exact integer test
   `(2*i-31)^2 + (2*j-31)^2 + (2*k-31)^2 <= 31^2`.
   Require both actual grids to map bijectively to exactly these indices. Reject
   zero-plane, missing, duplicate or additional points. Region labels follow the
   signs of these exact coordinates; +Z is far and -Z near.
5. Certify the union-of-spheres occupancy predicate at **every** included lattice
   point using outward binary64 interval arithmetic around the exact-decimal
   geometry. Each primitive operation expands endpoints using nextafter;
   division by an interval containing zero is rejected. Square root endpoints
   are rounded outward. A point is occupied if at least one atom's squared
   distance upper bound is <= its squared-radius lower bound; it is free only
   if all lower bounds are > their corresponding radius upper bounds. Otherwise
   the point is unresolved and the descriptor is not admitted. There is no
   uncertainty threshold or boundary-residual allowance.
6. Compare both tools' occupancy masks to the certified reference, by lattice
   index, for all 56 × 3 frames. All must be identical. Inspect native optimized
   regional outputs against its brute predicate probe as an additional check.
   This is stronger than equal occupied counts, which can hide cancellations.
7. Derive pi independently using Machin's identity
   `pi = 16 atan(1/5) - 4 atan(1/239)` with rational alternating-series bounds
   (80 terms each). Require both bounds to round to the same binary32 and
   binary64 constants. Do not obtain the sphere constant from StericX.
8. Evaluate the inspected per-tool count-to-output operation sequences with
   exact rational arithmetic, rounding each primitive to nearest, ties to even,
   at that tool's precision. Integer counts are exact at these grid sizes.
   StericX uses `(4*pi_f32*R^3)/3`, direct region fractions and sorted sum of
   per-frame values divided by three. Morfeus uses `(4*R^3*pi_f64)/3`, octant
   fractions and octant sums for quadrants, then sum/3 for frame means. Match
   each tool's percentage formula as written, not an algebraic simplification.
   Require exact bit equality between each tool and its own independent model.
   A disagreement, unsupported operation, nonfinite or subnormal result is
   unresolved; it cannot be absorbed into an epsilon. The exact-real target is
   the same rational occupancy fraction times `4*pi*R^3/3`; the different
   rounded results are thus precisely attributable to the declared arithmetic.

This reference resolves the finite-lattice estimator more accurately than either
candidate: exact sphere membership and certified occupancy, with independently
rounded scalar evaluation. It is **not a continuum convergence study**. Passing
would authorize only a comparison explicitly described as this fixed 32-point-
per-edge, density-0.01, three-frame estimator. No accuracy claim about the true
continuous union volume, alternate density, geometries or hardware is authorized.
If the intended scientific quantity requires a continuum error budget, that
portion remains unresolved pending a separate reference convergence study.

## Execution, investigation and publication

Freeze this file before running the new evaluator. Record and verify the hash at
start and end of each evaluation, along with hashes of reference code, observer,
current production sources and binaries. Preserve all 56 rows and every frame,
including failures, exact masks, interval ambiguity and molecule-level residuals.
Retain the original MAE, RMSE, median/max AE, identity-line R², slope and intercept
for every descriptor; relative error is undefined at a zero reference value.
Investigate each predicate/output failure separately; do not trim outliers.

No pass/fail rule is inferred from these summary statistics. A new implementation
of the reference may fix a demonstrable coding error only with an explicit
revision log and preservation of the failed run; changing a scientific rule
requires a new policy. If this candidate gate fails, stop that performance portion
and report unresolved, without further timing or speedup publication.

Even if it passes, the previous timing campaign remains exploratory. A performance
claim requires new alternating paired timings after scientific admission, the
existing fairness/measurement/documentation gates, and explicitly limited scope.
The old approximately 14× result and the internal 5.46× optimization are never
combined. No production implementation may be optimized during this task.
