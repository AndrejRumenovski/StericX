# StericX–DBSTEP buried-volume resolution

Investigation baseline: `f04b156f7b8b31700e1c6ab978e337b00cc734e0`.
Status: scientific and engineering validation complete; matched performance is
pending the quiet window. No release recommendation is issued by this
intermediate report.

## Original failure and root cause

The original strict failure concerns **SIG-NIHDA-401, conformer 0, R=2.5 Å**:
three native orientations, their mean volume/percentage, and the repeated
radius-profile position. It is one geometry/configuration, not six independent
chemical failures. Original requests, accepted atoms, center, radius settings,
all outputs, failed attempts and executable/environment hashes remain frozen in
[the baseline](resolution_baseline/manifest.json).

StericX selected 5,565 sphere nodes using binary32 squared-distance arithmetic.
The comparison then supplied those nodes to DBSTEP and clipped them a second
time using binary64 sphere membership. This removed 32 denominator nodes and
six occupied nodes per orientation. Exact rational margins show that those
nodes lie just outside the exact mathematical sphere; the native operation
rounds them onto its inclusive boundary.

The first divergence is that **second sphere test**. Atom mapping, molecular
frames, radii and atomic occupancy agree on the supplied nodes. The original
lane incorrectly asserted an identical finite estimator while changing its
sample domain. Its arithmetic-only tolerance cannot bound that change.

| Comparison, three-frame mean | StericX Å³ | Reference Å³ | Absolute difference Å³ |
|---|---:|---:|---:|
| Original DBSTEP reclip | 16.4300880432 | 16.4541364734 | 0.0240484302 |
| Corrected common sample domain | 16.4300880432 | 16.4300873655 | 0.0000006777 |

The corresponding common-domain percentage is 25.1033248901% versus
25.1033243486%, a difference of 5.41529e-7 percentage points. No original tolerance
was increased. The old failed comparison is still reproducible and remains
labelled failed in its historical lane.

**Classification:** reference comparison error with an underlying floating-point
integration difference. No StericX production defect was confirmed.

## Implemented correction and independent evidence

The corrected research comparison calls unmodified DBSTEP atomic occupancy on
exactly the supplied native points, without another sphere clip. The independent
reference uses direct union predicates, exact rational geometry, analytical
sphere intersections and adaptive box integration. It imports neither production
engine. A source-bound observer obtains native masks by calling the unchanged
Rust kernel, including separate aggregate calls that exercise its row cache.

Across **69 cases, 207 frames and 3,545,298 point decisions**, all atomic masks
agree exactly in Rust, DBSTEP and the independent reference. Native count
arithmetic reconstructs all aggregate and regional bits. All 138 public
mean-volume/percentage comparisons pass the original arithmetic limits. The
maximum frame volume residual is 7.86485133e-6 Å³.

All **53 scientific Rust/Cargo files**, Python scientific bindings and numerical
defaults are unchanged. Independently built baseline and candidate native
executables have the same SHA-256:

`3b9079338386a2532b17ebe844a28d72344a7bdf91d83526d4b3d0fd00438424`

The new implementation is the corrected comparison, transparent independent
reference, boundary regressions and evidence review. There is no artificial
production-code correction or hidden compatibility mode.

## Additional reference and convention findings

Official unmodified DBSTEP 2.2.0, commit
`4365353e0a2f24d05336e619759e350bc630a3cb`, has a separate volume-only fast-path
inconsistency: numerator and denominator sphere predicates use different
floating-point operation orders. A reproducible one-carbon example counts
510 occupied nodes over 491 denominator nodes, then clamps 103.86965% to 100%.
Its analytical continuous percentage is 98.42890625%.

This is a **confirmed DBSTEP reference implementation issue**, independently
demonstrated by exact rational counts and analytical geometry. It is not the
cause of the original ligand-401 tree-path failure. Reference sources remain
unmodified; the affected public fraction is excluded as an exact finite-set
oracle, with every raw result retained.

Exact external tangency reveals another legitimate finite-grid convention:
Rust and DBSTEP can count a contact node, while Morfeus's strict prefilter removes
the exactly tangent atom. All share zero continuous intersection volume. These
findings are detailed in the [root-cause analysis](ROOT_CAUSE_ANALYSIS.md),
[source review](resolution_analysis/DBSTEP_SOURCE_REVIEW.md) and
[convergence study](resolution_analysis/convergence_v2/CONVERGENCE_REPORT.md).

## Numerical convergence and accuracy limits

The prospective schedule contains **20 geometries × four resolutions × four
estimators = 320 comparisons**. Native/Morfeus nominal densities are
0.08, 0.01, 0.00125 and 0.00015625 Å³; DBSTEP spacings are 0.4, 0.2, 0.1 and
0.05 Å. The methods' actual grids differ and are retained.

All analytical enclosure and interval-nesting checks pass. Every one of 801
disputed native grid-node observations has an exact rational margin and matching
IEEE trace; none remains unexplained. Actual-cell error bounds account for
nonuniform spacing, finite coverage, rounding, occupancy, count normalization and
the reference clamp. They are deliberately conservative, not fitted tolerances.

| Original frame | Certified continuous interval Å³ | Default native Å³ | Finest native Å³ |
|---|---|---:|---:|
| 15 | [16.475986864, 16.609940067] | 16.418327332 | 16.525009155 |
| 13 | [16.475893615, 16.610025123] | 16.441848755 | 16.525650024 |
| 10 | [16.475974859, 16.609943748] | 16.430088043 | 16.525812149 |

Thus default-grid absolute errors are provably at least **0.034–0.058 Å³** for
these frames. The finest estimates lie inside the independent intervals, whose
widths remain about 0.134 Å³. These bounds cannot rank differences much smaller
than their widths. The original 0.00025 Å³ arithmetic limit is not a continuum
accuracy guarantee.

Among the fixed analytical primitives, maximum native absolute error decreases
from 1.26536 Å³ at the coarsest setting to 0.01358 Å³ at the finest. This is a
dataset-specific result. All **31 nonmonotonic refinement steps** are retained;
finer sampling is not automatically closer for every geometry.

## Scientific and engineering validation

The replay includes:

* 348 Rust tests and actual Linux/macOS/Windows Rust CI.
* 128,021 byte-identical observations across 13 legacy scientific lanes.
* Fresh Morfeus/Kraken comparisons for 31,721 volume inputs, 95,163 frames,
  1,141,956 regional bins and the 2,517-case fixed comparison set.
* 93 real CLI cases at each of 1, 2, 4 and 6 threads, with exact equivalent output.
* 7,175 model-equation and 449 supplementary domain checks; 744 kinetic scalars,
  51 aggregation checks and the retained thermodynamic/CLI cases. These test
  equations and contracts, not catalyst-prediction accuracy.
* DBSTEP upstream tests: 413 passed, three upstream skips.
* 61 conditional Sterimol cases, 12 density comparisons and all radius-profile
  outcomes, with the original strict profile failure retained.
* 42 valid and 11 invalid population cases; 11-ligand/56-conformer reference
  comparisons and installed-package conformer goldens.
* 49 molecular adversaries: 48 valid results and one retained collinear-frame
  error; permutations, rigid motions, explicit H and phosphine classes included.

Five fresh local installations pass 71 package tests each: Python 3.11–3.14
wheels plus a clean Python 3.12 source build. Actual remote CI passes all
12 Linux/Windows/macOS Python combinations plus source installation, each with
71 tests. Windows CRLF artifact differences are recorded with raw distinct
hashes and verified identical normalized source/ASTs; the initial byte-check
failures remain preserved.

The package remains **0.4.0.dev0**, backed by Rust engine **0.3.0**. Wheels install
successfully; no PyPI publication or new archived release is claimed.

Final standard-runner discovery finds **156 tests: 143 pass and 13 explicitly
require the reference environment**. Those skips are covered by **28 passing
pinned-reference tests and zero skips**. Global Ruff/format checks pass, and all 382 current source records
match before/after testing. An initially incompatible pytest test module was
converted to discoverable unittest methods without changing its numerical
assertions; its failed discovery attempt remains preserved.

The [final scientific admission](resolution_validation/admission_v2/review.json)
verifies 2,473 evidence files, including the original test source in its immutable
campaign snapshot and the equivalent current test bodies. It admits measurement,
with release explicitly pending the performance review.

## Performance and release decision

Pending the admitted matched performance run. The new baseline is the fresh
pre-correction f04b156f build, not the earlier historical throughput number.
The fixed protocol covers single-molecule latency, conformer processing, six-core
batch throughput, Python overhead, memory and allocation counts. Baseline and
candidate native and diagnostic executable bytes are identical.

The user explicitly [approved the prospective policy](POLICY_APPROVAL.md) before
its evaluation. Historical strict failures are preserved; new common-domain,
boundary, continuum, preservation and installation lanes are evaluated
separately. Final release recommendation requires the completed engineering and
performance receipts and will be recorded here after review.

## Remaining scientific limitations

Hard-sphere radii and inferred coordination centers are modelling choices.
Binary32 input/frame arithmetic and finite point grids have resolution-dependent
boundary behavior. No universal continuum precision or bitwise rotational
invariance is claimed. Regional ratios use their documented separate sample
counts and molecular frames.

DBSTEP's confirmed public fast-path issue remains in the unmodified upstream
reference. Historical Kraken geometry identity and unavailable historical
energies remain unresolved provenance limits; direct descriptor validation does
not reproduce those unavailable historical inputs. Density mode remains a
distinct thresholded-grid quantity with its existing validated grid restrictions.

All original failures, corrected controller attempts, input identities and
unfavorable comparisons remain reproducible. The linked study reports give exact
commands; the final commit ledger records the separate changes.
