# Buried-volume cross-validation protocol

Written before StericX–DBSTEP buried-volume comparisons. Use the pinned,
unmodified reference and a fresh source-bound native audit observer. The observer
exports actual native grid points and aligned atoms; it is the system under
test, not an independent oracle. The independent occupied-point and buried-
volume calculations are the official DBSTEP functions.

## Fixed cases and scientific inputs

Base cases are all 56 retained MMFF conformers and the five already frozen
Sterimol witnesses. Use their exact accepted f32 geometries, explicit centers,
donor and three plane-atom indices. Default volume settings are R=3.5 Å,
density=0.01 Å³, radius scale=1.17 and H excluded. Add, without residual-based
selection:

* H included for conformer 0 of ligands 401, 498 and 2067;
* densities 0.08 and 0.00125 Å³ for the forward synthetic witness;
* R=2.5 and 4.5 Å for conformer 0 of ligand 401;
* radius scale=1.0 for conformer 0 of ligand 401.

This gives 69 cases and 207 molecular-frame integrations. Keep all geometries,
indices, radii, grid settings and failures. The extra settings probe distinct
scientific inputs; they are not adjustments selected to improve agreement.

## Matched finite-lattice comparison

Export the native integration points once for each radius/density combination,
and each frame's actual aligned atom positions and squared-radius intermediates.
Restore serialized native numbers to their exact f32 values before expanding to
f64. For DBSTEP, pass the native scaled radius `f32(base_radius * scale)` as a
radius, **not** the square root of the rounded squared-radius intermediate.
This preserves the declared radius input and leaves differences in distance
arithmetic visible. Include exactly the atoms retained by the native H policy.

Call unmodified `sterics.occupied` on those points/atoms/radii and
`sterics.buried_vol` at the declared origin and radius. Record denominator and
occupied counts separately. Report both DBSTEP's occupied-point volume
`N_occ * spacing^3` and its percentage `100 N_occ/N_sphere`.

Native absolute buried volume has a different finite normalization:
`(N_occ/N_sphere) * (4 pi R^3/3)`, with native f32 arithmetic. Therefore compare
native absolute volume only with the DBSTEP fraction times the analytic sphere
volume, clearly labelled as a derived common-normalization quantity. Do not
compare it directly to DBSTEP's whole-molecule `occ_vol`, or silently relabel
grid-count volume as analytic-sphere-normalized volume. Average the same three
frame estimates to compare the native summary.

The independently fixed arithmetic gates are `32 * 2^-23 * analytic_sphere_volume`
for Å³ and `32 * 2^-23 * 100` for percentage points. Native public results must
also agree bitwise with their captured frame/summary outputs where the source
case/settings are identical. Any error beyond the arithmetic gate, denominator
difference or atom mapping discrepancy must be retained and investigated with
the actual point decisions; no empirical residual cutoff is fitted. This gate
validates the stated finite lattice, not the continuous volume.

## Independent lattice and convergence comparison

Separately run official `Dbstep.from_coords` hard-sphere volume on the same
source geometry, explicit center, sphere radius, scaling and H policy at the
fixed spacings 0.20, 0.10 and 0.05 Å. Use the original source coordinate system;
DBSTEP chooses its own lattice. Record its actual radii (including H 1.09 Å
where included), occupied-point and sampled-sphere volumes, and all percentages.
The H-included default-table lane has an explicit radius convention difference
from native H 1.20 Å and cannot be called the same hard-sphere surface.

Compare all common-definition quantities and retain every resolution, even if
convergence is nonmonotonic. Report MAE, RMSE, maximum/median absolute error,
R², slope/intercept and molecule-level residuals. There is no fitted numerical
PASS threshold for two different finite lattices. Classify a difference by the
earliest changed mathematical input (radii, aligned coordinates or lattice)
and preserve unresolved cases. Passing the matched-lattice gate cannot by
itself establish high-resolution continuum accuracy. Existing independent
Morfeus and analytic convergence audits remain separate evidence.

All new native builds, requests, raw responses, reference calls, stdout/stderr,
environment, helpers and source hashes are retained. No original scientific
constants, default settings or reference implementation may be changed to make
this comparison agree.
