# Proposed hard-sphere validation policy

Status: **proposed; awaiting explicit approval before release-gate evaluation**.
This policy replaces neither historical evidence nor scientific constants.
The root-cause study has measured inputs and outcomes, but has not scored a
release under this policy. The prospective analysis schedule is fixed in
[ANALYSIS_PROTOCOL.md](resolution_analysis/ANALYSIS_PROTOCOL.md).

## Why separate validation lanes are required

The previous strict lane compared a native already-clipped binary32 sphere grid
with the same coordinates clipped again by DBSTEP's binary64 sphere predicate.
Its requirement of equal denominators therefore did not describe an identical
finite estimator. An arithmetic-only tolerance cannot bound a changed sample
set. The original failed results remain failed in that historical lane.

The proposed correction is to test finite-domain occupancy, sphere selection,
continuous-volume accuracy and default-table conventions separately. It does
not increase the original arithmetic tolerance or change production defaults.

## Required gates

1. **Original failures remain reproducible.** All original requests, geometry,
   accepted coordinates, parameters, outputs and environment hashes must be
   preserved. Each failure must have an operation-level explanation or remain
   explicitly unresolved. Unresolved engine errors block scientific admission.

2. **Common-domain atomic occupancy.** Supply the exact observed native nodes,
   atomic positions and accepted scaled radii to unmodified DBSTEP occupancy and
   the independent union reference. Do not perform another sphere clip in this
   lane. Compare masks, not only total counts. For the original 69 cases/207
   frames, exact mask agreement is required. Any new disagreement must be
   investigated and cannot be hidden by aggregate tolerances.

3. **Floating-point boundary witnesses.** Exact rational squared-distance
   margins determine the mathematical classification of represented inputs.
   The independently derived binary32 error bound is
   `B = gamma_5*s + u*r²`, where `u=2^-24`, `gamma_5=5u/(1-5u)` and
   `s=sum((p-a)²)`. Its finite-normal domain checks are mandatory.
   A disagreement outside this band fails. A disagreement inside it is accepted
   as explained discretization only when a recorded IEEE operation-order trace
   reproduces the actual observed kernel decision and the independent rational
   calculation accounts for the alternative. Unsupported domains and untraced
   disagreements remain unresolved. Being inside a band alone never passes.
   Full derivation: [independent methodology](resolution_analysis/INDEPENDENT_METHOD.md).

4. **Normalization and regions.** Integer occupied/total counts and the declared
   quadrant/octant sign mappings must agree exactly with observed masks.
   Count-to-volume binary32 arithmetic must reconstruct the native output bits.
   The original `32*2^-23*4*pi*R³/3` volume and `32*2^-23*100` percentage
   tolerances additionally remain in force for matched-count, common-domain
   comparisons. Whole-grid DBSTEP `occ_vol` is not a buried-volume oracle.

5. **Independent continuous-volume evidence.** Analytical sphere intersections
   and adaptive-box bounds must use the same accepted balls and center. The box
   reference must enclose analytical cases and give nested bounds as refined.
   Run every fixed resolution; publish volume, percentage, computational cost,
   absolute/relative errors or error intervals and nonmonotonic behavior.
   Each claimed discretization error bound must follow from cell coverage and
   boundary geometry, including rounded nodes/predicates and ratio normalization.
   Check the analytical truths against those independently derived bounds.
   No residual-derived absolute tolerance and no universal default-resolution
   precision claim are permitted. Broad bounds remain broad, explicitly reported
   accuracy limitations; cross-tool agreement alone is insufficient.

6. **Adversarial geometry.** Preserve all requested single-ball, overlapping,
   contained, crossing, tangent, adjacent-float, symmetry, hydrogen, phosphine,
   permutation, rotation/translation and nearly degenerate cases. Permuting atom
   storage while retaining identity/topology must preserve results. Rigid-motion
   differences must be accounted for by accepted-coordinate/frame rounding and
   boundary-point changes; unexplained changes fail. Invalid frames must retain
   structured errors and must not be silently dropped.

7. **Independent-reference defects stay visible.** DBSTEP's separately confirmed
   mixed sphere predicates and percentage clamp are reproduced in an explicit
   reference-defect lane. Its inconsistent public fraction is excluded as an
   exact finite-set oracle for those cases, with raw counts and analytical truth
   retained. DBSTEP source is not patched to force agreement. This does not
   exempt native calculations from the independent reference or analytical tests.

8. **Preservation and engineering.** Complete Rust/Python tests, legacy
   Morfeus/Kraken scientific replay, conformer identities/statistics/populations,
   precision, density and profile tests must pass their existing applicable
   contracts. Original known historical reproduction limitations remain labelled.
   Actual Linux/macOS/Windows wheel installations and supported Python versions
   must pass. No source changes may escape the tested source/hash binding.

9. **Performance disclosure.** Run the established matched baseline/candidate
   benchmark on the fixed workloads, affinity and paired sampling schedule.
   Measure allocations, memory and Python overhead. Investigate the existing
   10% regression trigger without sacrificing a justified scientific correction.
   The historical 10,825 files/s is context, not a fitted acceptance target.

## Admission and claims

Admission requires approval of this policy followed by a separate, source-bound
review of all gates. Approval of the policy is not approval of unseen results.
No remaining unexplained native error may be relabelled as convention or
roundoff. A failed experiment stays in the record even if a later corrected
controller succeeds.

If these gates pass, the release claim is a reproducible hard-sphere numerical
estimator with independently validated finite-point semantics and documented
resolution-dependent accuracy. It is not exact continuum integration, universal
rotational bit-invariance, or equality with every default DBSTEP setting.
If a required gate fails, publication/release recommendation remains withheld.
