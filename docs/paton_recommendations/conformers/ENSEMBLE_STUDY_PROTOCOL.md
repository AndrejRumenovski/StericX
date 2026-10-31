# Ensemble behavior study protocol

Written before collecting the new independent reference outputs. The fixed
cohort is all 11 rows of `data/reactions_raw.csv`, all 56 retained geometries and
all recorded MMFF94 energies, as sealed in `sensitivity_v1/manifest.json`.
Nothing is selected by residual, descriptor range or experimental response.

## Comparisons

1. Compare single-conformer (first recorded state), arithmetic mean,
   MMFF-energy weighted mean, minimum and maximum reductions of each of the
   14 scalar fields with the **same reduction** of fresh Morfeus 0.8.0 outputs.
   Retain every per-conformer and per-ensemble residual and failed reference
   call. This tests direct numerical behavior, not predictive usefulness.
2. For the fixed historical endpoint `vbur_max_delta_qvbur_min` in the retained
   official Ni-hDA/Kraken CSV, compare the five native estimates without fitting
   any correction or selecting a winner for deployment. Only the native
   minimum shares the target reduction definition. Other reductions are
   explicitly different estimators, included to evaluate sensitivity to the
   choice of ensemble summary. This is exploratory historical reproduction;
   original geometries/energies are unavailable. It is not direct descriptor
   validation, a causal test of conformer averaging, or an experimental
   prediction benchmark. No training, hyperparameter selection, feature search
   or held-out prediction is performed. The historical 723/2064 provenance and
   modeling limitations are unchanged.

Report N, missing pairs, MAE, RMSE, maximum and median absolute errors,
coefficient of determination `1-SSE/SST`, and least-squares slope/intercept
(native against reference, descriptive only). Constant-reference R² and slopes
are unavailable. Do not use squared correlation as the R² error metric. Preserve
all molecule-level values, including unfavorable and unavailable comparisons.

## Common scientific inputs and numerical triage

Use the source f32 coordinates and radii accepted by Rust, expanded exactly to
f64, with unchanged atom order. The explicit virtual center and three plane
atom indices from the native report are supplied to both methods; this isolates
descriptor kernels and does not independently validate center inference.
Use the recorded 2.28 Å center distance, 3.5 Å sphere, 0.01 Å³ density, 1.17
radius scale and H exclusion for volume. Morfeus uses its own alignment and
grid arithmetic. Total/near/far volumes are averaged over all three frames,
matching the current symmetric native definition; regional extrema span all
frames. Do not reuse the older first-frame 2.1 Å study outputs as if identical.

Sterimol includes all real atoms, excludes the zero-radius virtual center,
and uses +0.40 Å for L. Native default B1 uses 360 unique directions; Morfeus
uses 3600 endpoint-inclusive directions. Pyramidalization uses the same three
declared donor neighbors. Record all configuration and reference-package hashes.

For triage before observing residuals, define epsilon = 2^-23 and a length
scale M = max(1 Å, max_i distance(atom_i,center)+radius_i). L/B5 differences
above `64 epsilon M` require investigation. The B1 threshold additionally
includes the support-function angular bound
`rho_max * (pi/360 + pi/3599)` for the two independently phased finite scans.
For volume values use `32 epsilon * (4 pi R^3/3)` as the floating-reduction
budget, and `32 epsilon * 100` for percent. Differences beyond that budget
remain explicit grid/rounding investigations, not a fitted tolerance or a
continuum-accuracy PASS. Pyramidalization triage uses 1e-5 dimensionless and
1e-3 degree on this nondegenerate trivalent cohort; these are reporting
thresholds, not universal precision guarantees near singular geometries.

The ensemble numerical budget is the maximum constituent budget (the five
reductions are nonexpansive for fixed nonnegative normalized weights), plus
binary64 reduction roundoff. The independent energy-weight calculation is
the already tested Decimal reference. Failures are preserved rather than
discarded from an ensemble mean. Historical comparisons have no numerical
acceptance threshold because they do not share original scientific inputs.

No claim that ensembles are always superior is permitted. Small spread does
not prove rigidity; one-state ensembles and repeated/symmetry-related retained
states stay visible. MMFF energy weights omit thermochemical contributions and
do not establish experimental equilibrium populations.
