# Scientific audit

Status: campaign in progress. This document distinguishes confirmed defects,
tested behavior, reference conventions and outstanding validation. Existing
scientific-remediation and DBSTEP failure evidence remains unchanged.

## Confirmed correction S1: R² depends on target units

`model::fit::metrics` discarded target variance below 1e-12;
`model::evaluation::score_frozen_predictions` discarded variance below machine
epsilon. Both thresholds have absolute squared target units. They incorrectly
returned `None` for small, nonconstant targets.

Independent witness: actual targets (-s,+s), predictions (-s/2,+s/2).
The mean is zero, RSS=s²/2 and TSS=2s², so R²=1-RSS/TSS=3/4. Using dyadic
s=2^-30 makes every operation exactly representable. Both original functions
failed this test; see [r2_baseline_failure.log](r2_baseline_failure.log).

The correction tests TSS>0 and evaluates the ratio lazily. Constant targets
still return `None`; poor predictions can still have negative R². No fit
regularization, physical constant, prediction or scientific tolerance changes.
Regression coverage checks s=1, 2^-30 and 2^-100, and the fit helper also checks
constant targets and the exact score -3. All 58 model unit tests and six public
training API tests pass, including the unchanged golden training report and
bit-for-bit comparison against the former CLI composition. Raw validation logs
are in `.stericx/ultimate_engine/r2_{model_tests,training_api}.log`.

## Audit coverage and conventions

The remaining entries are being expanded with fresh independent replay evidence:

| Calculation family | Definition and audit focus | Current evidence/status |
|---|---|---|
| Sterimol | Raw L versus external +0.40 Å convention; signed support B1 on 360 directions; radial B5; f64 projection with f32 inputs/results; explicit dummy and axis | Historical remediation read; fresh source-bound baseline captured; continuous-angle method investigation pending |
| Hard-sphere buried volume | Union membership on the exact f32 sphere lattice, radii scaling, H exclusion, three donor planes, independently normalized regions | Fresh analytical convergence and native baseline captured; corrected DBSTEP common-node replay pending |
| Quadrants/octants/hemispheres | Zero belongs to the positive side; each region uses its own finite-grid denominator; near is negative local Z | All bins in baseline replay captured; regional sums need not equal total on odd grids |
| Coordination geometry | Negative sum of unit bonds, trivalent topology, planar sign ambiguity requires explicit center | Existing ambiguity/degeneracy gates preserved; geometric lone pair is not an electronic orbital |
| Pyramidalization | Plane geometry and angular descriptors, proper degeneracy rejection | Existing tests and fresh baseline captured; independent replay pending |
| Ensembles/statistics | Identity/order, minima/maxima/ties, weighted means/distributions, missing-data provenance | Fresh aggregation lane captured; population reference replay pending |
| Populations/thermodynamics | Declared energy units, reference-energy subtraction, temperature, stable exponential normalization | Existing explicit missing-energy policy preserved; independent replay pending |
| Electron density | Explicit total density, Bohr conversion, grid transforms, strict density > threshold, cell-volume integral | 12 fresh independent DBSTEP/analytic Gaussian cases pass; convergence remains discretization-dependent |
| Radius profiles/sensitivity | Fixed center, requested radius order/duplicates and per-radius failures | Focused profile tests pass; fresh reference replay pending |
| Eyring/ee/reaction calculations | Energy/temperature units, signed conventions, stable extreme limits and uncertainty | Fresh kinetics lane captured; independent equation replay pending |
| Predictive models | Feature conventions, AVX2/scalar reduction differences, nested versus fixed-feature validation, intervals, extrapolation | S1 corrected; 64 focused tests pass; complete model reference replay pending |
| Input/radii/storage | Coordinate precision, element/radius policies, parser errors, record serialization | Frozen source/data hashes and exact parser replay retained; public radius fallbacks are conventions, not fitted constants |

Finite-grid agreement is not a certificate of continuous-volume accuracy. The
fresh eight-geometry study shows nonmonotone convergence and an external near
contact that can be missed entirely. The current default remains the validated
finite-grid definition. Density is a separate scalar-field descriptor; it does
not validate the hard-sphere physical model or a universal density threshold.

Only x86_64 hardware is available. Cross-architecture bitwise equivalence is
not established by the current runs. The AVX2 model predictor and its scalar
fallback use different reductions, so identical real-arithmetic formulas alone
cannot establish identical floating-point bits.
