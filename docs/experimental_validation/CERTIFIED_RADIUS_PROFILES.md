# Certified radius profiles

`sphere_union_volume_certified_radius_profile(center, atoms, radii, options)` returns
total, quadrant and octant certificates at every requested radius. The center and
Cartesian frame remain fixed. `CertifiedBuriedVolumeFrame::compute_radius_profile`
uses the same exported molecular frame at every radius; it does not re-infer a
radius-dependent coordination center. Each entry retains requested float bits,
result or error, and any exact-duplicate reuse index. Input order is unchanged.

Atoms are validated, canonically ordered and deduplicated once. Atomic relative
coordinate and squared-radius intervals are formed once and reused because these
quantities depend only on the fixed center and atomic geometry. The integration
sphere's interval and all geometric candidate decisions are recomputed for each
distinct radius. Every distinct radius receives its own entire cell budget and
independent XY tree. No clipping or numerical volume bound crosses unequal radii.
Bit-identical duplicate radii, with all common inputs identical, reuse their exact
certificate; their recorded integration time is zero. This reuse is mathematically
identity of the same integral, not an assumption about neighboring radii.

Work counters in the enclosing profile count actual integrations only. An entry's
work counters describe its certificate even when reused. Whole-profile runtime
includes common preparation; entry runtime excludes that shared preparation.
Invalid radius entries (nonpositive/nonfinite/>1e100) remain errors rather than
aborting other points. Invalid common geometry/options abort the whole profile.
An empty requested profile returns zero entries and zero integrations.

## Mathematical consistency

For fixed geometry/frame and R1<R2, B(c,R1) is a subset of B(c,R2). Intersecting
both with the same atomic union and same quadrant/octant preserves inclusion, so
**absolute occupied volume** is nondecreasing separately in every region.
The percentages divide by a denominator proportional to R³ and need not increase.
A centered radius-1 atom gives 100% at smaller integration radii and 12.5% at R=2.
Neither finite-accuracy midpoint estimates nor independently computed lower/upper
bounds need be individually monotone. A contradiction would require a smaller
radius's certified lower volume to exceed the larger radius's upper volume. All
such checks use the unmodified intervals; no monotonic postprocessing hides errors.

An optional future monotonic interval intersection could transfer lower bounds
from smaller radii and upper bounds from larger radii with a valid inclusion proof.
It is not needed for the implemented independent-certificate contract. Reusing a
subdivision tree would require reclassification of each cell and radius-dependent
sphere clipping; reusing its old volume bounds would be invalid in general.

## Independent validation and observed cost

`validate_profiles.py` evaluates **24 profiles / 120 requested radii / 96 distinct
integrations**, including twelve analytical profiles and the twelve frozen molecular
cases. Analytical radii are 0.5/1/1.5/2.5/3.5 Å plus a repeated 1.5; molecular radii
are 2.5/3.5/4.5 Å plus a repeated 3.5. Analytical tolerance is 0.02 Å³ with 262,144
cells, molecular tolerance 0.05 Å³ with 1,048,576 cells. All targets succeed.

All **936 independently integrated analytical total/regional volumes** lie in their
certificates. Each returned certificate exactly matches a separate regional call,
including all bounds, flags and work counters. Molecular total intervals overlap
separate certified-total calculations; this is an independent algorithm-path
consistency test, not an exact molecular truth. Every absolute-volume monotonicity
consistency check passes. Requested-order, errors, exact duplicates, common-input
validation and the nonmonotone percentage example have additional Rust tests.

Shared-workstation diagnostics: median profile runtime 65.43 ms, maximum 442.32 ms;
maximum work at one radius 266,450 cells. Median shared preparation fraction is
0.00000846 (0.000846%). Whole validation observer, including repeated standalone
controls, peaks at 23,184 KiB RSS and uses 6.99 s user CPU. These are exploratory
costs, not controlled equal-accuracy speedup claims. The small measured preparation
fraction gives no basis to claim a material speedup across distinct radii from
preprocessing alone. Duplicate reuse avoids an entire repeated integration.

Raw evidence and hashes: `profiles_independent_v1/` under the campaign evidence
root and `profile_validation.json`. Early test compilation errors are retained.
An initial test incorrectly expected <2 percentage-point uncertainty from a
0.02 Å³ absolute tolerance on a tiny integration ball; the test's requested volume
tolerance was **tightened to 0.001 Å³**, retaining its percentage assertion. The
failed attempt remains in `profiles_v1/`; no original scientific tolerance was
relaxed. Full regression and final profiling are separate gates.
