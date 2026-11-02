# Prospective discrepancy analysis

This protocol is written before the new mask/convergence campaign. It fixes the
inputs and analysis rules; it does **not** replace or score the historical release
gate. All original failures and rejected executions remain immutable.

## Inputs and resolution schedule

The original 69 molecular/configuration cases (207 native frames) are selected
by the unchanged `validate_dbstep_volume.load_cases`. Every accepted atom, radius,
center, frame and point is saved. Atomic radii are the actual accepted binary32
base radii multiplied once in binary32 by the declared scale. No radius fitting
is permitted. Hydrogen-included default-table differences remain a separate
scientific surface from common-radius comparisons.

The convergence densities are **0.08, 0.01, 0.00125 and 0.00015625 Å³** for
StericX and Morfeus. DBSTEP spacings are **0.4, 0.2, 0.1 and 0.05 Å**. These
controls do not imply identical point sets. Report actual counts, volumes,
percentages, runtimes, errors and successive differences at every level;
center-point integration need not converge monotonically.

Continuous-volume bounds use the independent adaptive-box reference at depth
**11** for the original ligand-401 R=2.5 witness and depth **8** for analytical
primitive cases. The stopping depth is fixed before results. Publish interval
widths, including widths too broad to distinguish the original 0.024 Å³ residual.
Do not call interval inclusion evidence of finer accuracy than that interval.

Primitive cases include an empty union; a single enclosed, crossing, disjoint,
covering and concentric ball; internal/external tangency and adjacent binary32
positions; separated, overlapping and coincident balls; an asymmetric union;
and exact/adjacent sample-boundary points. Analytical sphere-intersection and
two-ball-union formulas apply only where their geometric assumptions hold.
Use explicit radii and ghost centers in reference APIs. Synthetic geometry is
labelled mathematical test data, not an experimentally characterized molecule.

Molecular adversaries include primary, secondary and tertiary phosphines and
PH3 with explicit hydrogens, symmetric/asymmetric geometries, fixed atom
permutations, rigid signed-axis rotations, a seeded general rotation and
translations, and nearly degenerate donor frames. Preserve expected errors
and every unsuccessful comparison. Seed: **20261102**.

## Point decisions and arithmetic

The observation adapter calls the unmodified Rust kernel to obtain individual
decisions and separately calls its aggregate kernel. Its copied scientific
source prefix must match production byte for byte. It is an observation of the
implementation, not a reference algorithm.

Compare each native atomic occupancy mask with unmodified DBSTEP `occupied`
on **exactly the supplied point domain**, without a second sphere clip, and with
the independently written direct union predicate. Also record the historical
DBSTEP reclip and its original failed criterion. Equal total counts are
insufficient when individual decisions differ.

The a priori binary32 certificate is specified in `INDEPENDENT_METHOD.md`.
Any disagreement outside its proven roundoff band fails. Every ambiguous
disagreement retains an exact rational margin and explicit IEEE operation
trace; an ambiguous classification alone never passes. Input/radius mismatch,
unsupported arithmetic domains and unexplained decisions are failures or
explicitly unresolved. No residual-derived numeric tolerance is introduced.

Validate integer counts and quadrant/octant sign mappings exactly. Reconstruct
the final arithmetic from observed counts, recording binary32 bit patterns;
this diagnoses normalization only and is not an independent occupancy oracle.
For matched counts, retain the original 32-epsilon arithmetic tolerance as an
additional check rather than expanding it to hide changed sample membership.

## Continuous integral and discretization

Point-domain correctness and continuum accuracy are separate questions.
Analytical and adaptive-box references use a different mathematical method.
Record each method's error or interval of possible errors; do not use R² or
cross-tool agreement alone to infer correctness. A theoretical cell-boundary
error estimate must account for actual node spacing, rounding, sphere clipping
and ratio normalization. An ideal-grid shell bound must not be silently applied
to rounded/nonuniform grids.

The independently exposed DBSTEP mixed numerator/denominator predicate issue
is tested separately. Its clamped public value cannot serve as an exact
finite-set oracle. Preserve its public output, unmodified reference sources,
raw counts, consistent diagnostic counts and analytical truth.

## Interpretation and release

Classify all disagreements as confirmed implementation error, convention
mismatch, numerical integration difference, reference comparison error, or
unresolved, with the first differing operation. Change production only when
independent evidence establishes an implementation error. Preserve validated
defaults when the difference is legitimate discretization.

A proposed release policy must be separately documented and approved before
being evaluated. This study does not silently turn the historical failed
equality test into a passing result. Full scientific replay, Python/conformer
validation, actual cross-platform wheel tests and fresh matched performance
measurements remain required regardless of the point-level explanation.
