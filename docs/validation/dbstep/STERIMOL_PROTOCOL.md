# Sterimol cross-validation protocol

This protocol precedes the first StericX–DBSTEP numerical comparisons. The
reference is the unmodified source and environment in `environment.json`.

## Fixed inputs and lanes

Use every retained geometry in the 11-ligand/56-conformer sensitivity manifest,
without selecting by residual. The native source f32 coordinates/radii and
explicit centers are expanded exactly to f64 for the reference. Preserve source
atom order, identities and all original files. Supply a zero-radius `Bq` ghost
at the common center to DBSTEP, index 1; the donor is native index + 2. The
ghost is not a ligand atom. Neither program reoptimizes the geometry.

Add five declared synthetic support witnesses: a forward trivalent P/C
geometry, that geometry with a distant backward carbon, an equilateral planar
P/C projection, a rigidly rotated/translated backward witness, and a permuted
backward witness with the donor index mapped explicitly. They test definitions,
symmetry, rotation and atom mapping; they are not purported optimized molecules.
Freeze their XYZ text, options and hashes before their native/reference runs.

Run two labelled reference lanes for every case:

1. The official `Dbstep.from_coords` API with classic vdW measurement, Bondi
   radii, scale 1, H included, no metal exclusion, default independent rotation,
   and the declared atom1/atom2. Retain its native default radii, including H
   1.09 Å; this tests workflow/convention differences.
2. The unmodified `sterics.get_classic_sterimol` function with explicitly
   supplied native radii and independently DBSTEP-aligned coordinates from
   `calculator.rotate_mol`. The function's radii argument is used directly;
   no constants are patched. This isolates support calculations under the
   common declared radii convention (including native H 1.20 Å). Record the
   reference-aligned coordinates and radius arrays for each lane.

Native coordination L includes +0.40 Å. Remove that documented correction
from the reported result for comparison with DBSTEP raw L, retaining both
original and compared numbers. DBSTEP uses absolute axial extent; native uses
signed positive-axis extent. Label L as sharing the definition only where the
forward envelope dominates the backward envelope. Backward-dominated cases
remain in the all-case residual table and convention-discrepancy list; they
must not be called failed implementations or silently deleted. Compare signed
L and both widths with unmodified Morfeus on the same common inputs to help
locate the earliest divergence.

## Acceptance established before observation

Let epsilon = 2^-23, M = max(1 Å, max_i(norm(atom_i-center)+radius_i)), and
rho_max be the maximum perpendicular distance from an atom center to the
declared axis. The arithmetic length budget is
`64 epsilon M + 2 sqrt(3) * 1e-8 Å`; the second term accommodates DBSTEP's
documented eight-decimal coordinate rounding. This is a conservative numerical
regression gate for these nonsingular fixtures, not a continuum precision claim.

For common-definition L and B5 in the supplied-radii lane, require absolute
error no greater than that budget. For B1 additionally allow
`rho_max * (pi/360 + pi/360)`, the support-function Lipschitz bound for two
independently phased one-degree scans. Native uses 360 unique directions;
DBSTEP uses 361 endpoint-inclusive directions. Do not fit an empirical cutoff
to residuals. Keep the native default scan unchanged.

In the default-radii lane, additionally bound the radius-convention contribution
by `max_i abs(native_radius_i - dbstep_radius_i)`. Report that bound and every
actual radius difference independently of the measured descriptors. The
analytical support functions are 1-Lipschitz in the maximum radius change.
Backward-dominated L is a declared convention difference; its discrepancy is
explained with the two explicit axial envelopes rather than an expanded error
tolerance. Any unexplained excess remains a validation failure requiring
investigation and frozen evidence.

Report MAE, RMSE, maximum and median absolute error, `1-SSE/SST` R², slope,
intercept and every case residual, both for all executed cases and clearly
labelled common-definition subsets. Constant-reference R²/slopes are null.
Retain errors, nonfinite results, missing pairs, atom mapping, stdout/stderr,
input/executable/source hashes, exact commands and all reference settings.
Neither R² nor a large angular bound alone proves high-resolution accuracy.
No production change is justified merely by a smaller DBSTEP residual.
