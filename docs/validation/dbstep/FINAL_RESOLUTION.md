# StericX–DBSTEP buried-volume resolution

Investigation baseline: `f04b156f7b8b31700e1c6ab978e337b00cc734e0`.
Status: applicable scientific, performance, installation and release gates pass
under the approved policy. **Python 0.4.0rc1 is prepared and unpublished.**
Original scientific failures and unsuccessful infrastructure/reviewer attempts
remain preserved; none has been rewritten as a successful historical result.

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

Per-conformer results, ensemble statistics and sensitivity reports, explicit
energy provenance/populations, density grids, radius profiles, Python adapters,
wheel packaging and the executed research notebook already existed at this
investigation's baseline. They were reused and independently revalidated. This
investigation adds numerical evidence and a corrected validation contract, then
prepares versioned distribution artifacts; it does not duplicate those features.

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

These completed scientific and performance measurements used **0.4.0.dev0**,
backed by Rust engine **0.3.0**. The subsequent **0.4.0rc1** preparation changes
four version declarations only; the engine, binding logic and package tests are
unchanged. Its separately identified artifacts pass five new local installations
and all 13 new platform/source CI jobs, each with 71 tests. One selected Windows
wheel additionally passes all 71 tests on each of Python 3.11–3.14. Linux/macOS
runtime payloads are byte-identical across their four interpreter builds. No
PyPI publication or new archived release is claimed.

Final standard-runner discovery finds **156 tests: 143 pass and 13 explicitly
require the reference environment**. Those skips are covered by **28 passing
pinned-reference tests and zero skips**. Global Ruff/format checks pass, and all 382 source records in that replay
match before/after testing. An initially incompatible pytest test module was
converted to discoverable unittest methods without changing its numerical
assertions; its failed discovery attempt remains preserved.

The [scientific admission](resolution_validation/admission_v2/review.json)
verifies 2,473 evidence files, including the original test source in its immutable
campaign snapshot and the equivalent current test bodies. It admits measurement,
with release explicitly pending the subsequent performance review. The
[completed release review](resolution_release/admission_v2/review.json) now binds
the performance results, exact four version substitutions and tested RC artifacts,
rechecking 2,703 evidence records. Final style/research checks bind 600 source
files before and after execution, with no changes.

## Controlled performance

The user-confirmed quiet window ran **2026-11-02 18:52:46–18:57:08 UTC**.
All 504 native attempts, 28 instrumented allocation runs, 96 Python timing
blocks and two ensemble allocation runs completed. Every raw output and summary
was rechecked; no sample was excluded and no campaign was repeated. The fresh
pre-correction f04b156f baseline and candidate have identical native/diagnostic
executable bytes and identical scientific results.

| Measurement | Fresh baseline | Candidate |
|---|---:|---:|
| Six-thread batch, files/s | 10,656.50 | 10,774.50 |
| Six-thread 10,000-file batch, median ms | 938.394 | 928.118 |
| One-thread small molecule, median ms | 1.889 | 1.981 |
| Six-thread 56-conformer batch, median ms | 7.646 | 7.541 |
| Six-thread batch peak RSS, median MiB | 26.334 | 26.457 |
| Python bridge, six-conformer call, median ms | 3.2727 | 3.2475 |
| Python typed API, six-conformer call, median ms | 4.0311 | 4.0125 |
| Warmed ensemble allocation calls | 871 | 871 |
| Warmed ensemble requested bytes | 2,454,897 | 2,454,897 |

No native or Python method reaches the unchanged 10% investigation trigger.
The largest native side-median slowdown is 4.843%, on the approximately 2 ms
small-molecule process. Maximum median RSS increase is 1.77%. All 14 native
allocation/reallocation/deallocation-count comparisons match exactly; small
process-wide requested-byte and peak differences remain reported. The warmed
ensemble call has 73 reallocations and zero live-byte growth on both sides.

Performance is preserved on this fixed workload set. Since executable bytes are
identical, timing differences are host/sample variation, not an algorithmic
speedup. Python timings exclude input parsing and use the same installed package
bytes; native-to-wheel differences also include build-target effects. Read
[the complete performance report](resolution_performance/PERFORMANCE.md) for all
28 configurations, distributions, scope, host counters and preserved samples.

## Release decision

The user explicitly [approved the prospective policy](POLICY_APPROVAL.md) before
its evaluation. Historical strict failures are preserved; new common-domain,
boundary, continuum, preservation and installation lanes are evaluated
separately. The original gate incorrectly compared two different finite sample
domains. Its historical failure remains; exact common-domain occupancy and
independent boundary/convergence evidence satisfy the approved policy without
weakening any numerical tolerance. No unexplained StericX implementation error
remains in the tested cases.

**Recommendation: admit the scoped Python 0.4.0rc1 release candidate for user
review.** The [candidate bundle](resolution_release/candidate_v1/candidate_artifacts.tar.gz)
contains Linux x86-64, Windows x86-64 and macOS arm64 ABI3 wheels plus a source
archive. Python 3.11–3.14 installations are validated. Exact artifact identities,
installation instructions and reproducible review commands are in
[release preparation](resolution_release/RELEASE_CANDIDATE.md). The local unpacked
bundle is `.stericx/releases/stericx-0.4.0rc1/`. No PyPI upload, GitHub Release,
public release tag or registry-installation claim is made.

The first RC macOS Rust job never acquired a runner: GitHub reported arm64
capacity constraints. Its overall run remains `failure` and its Mac job remains
`cancelled`; Linux/Windows jobs passed. A separate actual run on macOS 14 arm64
executes the identical source and commands, passing **347 applicable Rust tests
plus the 22-test screening repeat**, with zero failures/ignored tests and unchanged
source. The x86-only AVX2 comparison accounts for the one-test difference from
the 348-test x86 suites. The final gate explicitly combines these real platform
executions; it does not fabricate successful metadata for the original run.
[CI evidence](resolution_release/RELEASE_CI.md) retains the original capacity
failure, successful retry and deliberately cancelled duplicate pipelines.

The first final coordinator review rejected an incorrect expectation that the
overall workflow conclusion would also be `cancelled`. Frozen API records and
independent review agree that it is `failure`, while the unassigned Mac job is
`cancelled`. The corrected check preserves both fields and still requires every
test in the separate successful execution. The [failed reviewer attempt](resolution_release/admission_v1/failure.json)
and its source snapshot remain alongside the passing review. No acceptance
tolerance or test requirement changed.

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
commands; the [commit ledger](COMMITS.txt) records the separate changes through
release preparation. The final documentation commit is identified by Git history
and the completion response.

The [requirement-by-requirement completion audit](COMPLETION_AUDIT.md) and fresh
2026-11-12 integrity check confirm the same release decision. All 2,703 evidence
records and the complete candidate archive remain unchanged.
