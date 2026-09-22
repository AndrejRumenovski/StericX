# C3: reject impossible atoms once per integration row

Status: proposal only. C1 and C2 must be accepted and freshly profiled before any
C3 production edit. This proposal uses the corrected pre-C1 baseline at commit
`1ecfbd5be8b354711bdadd8a8d148417e03d8e7a`; its measurements do not describe an
accepted C1/C2 build. No C3 native speedup has been measured.

## Evidence and conditional benefit

The [corrected profile](../CORRECTED_PROFILE.md) attributes 81.773% of one-thread
`descriptors_10000` diagnostic process wall time to `occupied_volumes`. Its native
median is 4,991.033 ms. The [source-bound shadow report](../evidence/corrected_occupancy_shadow/RAW_REPORT.md)
replayed the exact corrected workload arguments for all 10,000 files and the
56 distinct conformers they repeat. It kept the original predicates and scientific
results, adding counters and assertions around a hypothetical row rejection.

| Measured counter | Diverse 56 | Actual 10,000 files |
| --- | ---: | ---: |
| Frame orientations | 168 | 30,000 |
| Integration point visits | 2,588,544 | 462,240,000 |
| Contiguous XY rows | 124,320 | 22,200,000 |
| Original lazy XY preparations | 3,837,215 | 685,507,811 |
| Original Z predicates | 10,884,695 | 1,943,832,492 |
| Removable Z predicates | 5,431,482 | 969,969,725 |
| Added endpoint predicates | 484,253 | 86,487,808 |
| Row/atom rejections | 335,743 | 59,966,781 |
| Net Z predicates removed | 4,947,229 | 883,481,917 |
| Net Z reduction | 45.4512% | 45.4505% |
| Net combined XY/Z comparison reduction | 33.6045% | 33.6009% |

Every frame uses 15,408 accepted integration points in 740 rows. No frame had a
net Z-predicate increase. The complete [summary](../evidence/corrected_occupancy_shadow/summary.json)
also records XY rejections, cached/new Z tests, first hits, misses, per-atom
histograms, nonfinite fallbacks, and a separate eager row×atom upper-work model.
The proposed algorithm retains lazy atom preparation; the eager model is not its
behavior.

A planning hypothesis is **15–25% less one-thread descriptor wall time** (about
1.18–1.33×), conditional on the post-C1/C2 profile still showing this work. This
requires an occupancy improvement of approximately 18.34–30.57% at the recorded
81.773% share. A deliberately simplified equal-cost comparison model gives
`0.81773 × 0.336009 ≈ 27.48%` whole-process reduction before unmodeled work; this is
not a rigorous ceiling or a calibrated prediction. The 15–25% hypothesis discounts
that proxy for row scans, finite checks, branches, cache access and other occupancy
bookkeeping; its upper end remains optimistic. Comparison counts are neither CPU
instructions nor wall time. The actual gain can be smaller or negative.

Acceptance requires reproducible **at least 10% lower target batch wall time**
beyond paired run variation, exact scientific equivalence, and no material memory
cost or regression in other workloads. At the old one-thread share, 10% would
require approximately 12.23% occupancy improvement. Fresh profiling must replace
these estimates before implementation and measurements must decide acceptance.

## Proposed algorithm and exactness proof

Use an outer loop over contiguous points whose X and Y `to_bits()` values match.
Scan the actual row once to obtain minimum/maximum Z and whether every point is
finite; preserve the original inner point order and regional counters. This can
also replace the existing per-point row-change branch. Do not assume a sorted or
regular Z sequence, merge separated equal-XY rows, sort points, or reorder atoms.

Keep the original lazy atom prefix/cache. When an atom first reaches the original
XY-survival point, retain the bit-identical f32 value
`q = dx*dx + dy*dy`. A row-level rejection is allowed only when the row points,
atom position, squared radius and `q` are finite and the atom's Z is strictly
outside the row's actual `[min_z, max_z]`. Test the nearest actual endpoint using
exactly `q + dz*dz <= radius_squared`. An endpoint miss then rejects the atom for
that entire row. Equality, an atom inside the range, or any nonfinite case retains
the original path.

Let `RN` denote the original f32 rounding and
`F(z) = RN(q + RN(RN(z - atom_z) * RN(z - atom_z)))`. All row points lie on the
same side of an outside atom. Rounded subtraction magnitude is monotone with
separation; squaring and addition of nonnegative operands are monotone. Therefore
`F(z) >= F(nearest_endpoint)` for every actual row point. A strict endpoint miss
implies every original predicate misses. Finite subtraction/square overflow to
infinity and underflow/absorbed additions preserve this order. An endpoint hit
proves nothing about other points and must retain the atom.

Preserve the original `(dx² + dy²) + dz²` grouping: no reassociation, FMA, f64
substitution, square-root boundary, or `dz² <= radius² - q` inversion. For example,
`q=1`, `dz=1e-4`, `radius²=1` is an original hit after absorbed f32 addition but
fails the inverted predicate. Nonfinite fallback also preserves NaN comparisons
and `infinity <= infinity`. Identical XY bits preserve signed-zero distinctions.
Removing only predicates proven false leaves the first occupied atom, lazy prefix
extension, point visits, grid/regional denominators and final reductions unchanged.
This holds for arbitrary point order within and across contiguous rows.

## Measured versus static operations

The table above and raw per-frame streams are measured shadow-loop counters. All
452,520 recorded f32 fields across 30,168 frames matched the unchanged occupancy
function bit for bit. Every would-skip predicate was still executed and asserted
false. Original native and probe stdout/stderr matched byte for byte, and both
matched the frozen native fingerprints. Two focused test groups exercised finite
boundary-adjacent values, absorption, subnormal underflow, finite overflow, inside-
range and reversed/interrupted rows, atom order, empty inputs, NaN and infinities.
The [completion binding](../evidence/corrected_occupancy_shadow/completion_binding.json)
connects source, build, tests, execution and binary hashes.

Coordinate transforms, radius operations and allocations were **not instrumented**
by this probe. The [static operation model](../evidence/corrected_occupancy_shadow/static_operation_model.json)
applies source operations to observed shapes: 935,751 aligned-atom transforms,
2,807,253 dot products, and 935,751 radius scales/squares for the 10,000-file
workload. Parsing its 553,228 input atoms entails that many source-level radius-
table lookups. There is no transform or element-radius lookup in the point loop.
The source has 10,000 grid-construction sites, 30,000 aligned-vector collections
and 30,000 initial occupancy-cache `Vec::with_capacity` calls. These are not measured
allocator calls, reallocation counts, bytes, or instruction counts. Actual C3 memory
and allocation effects require their own diagnostic measurement.

A candidate can keep scalar row bounds on the stack and reuse the existing
candidate vector, with no required additional heap allocation. Do not copy the
shadow's row vector, counters, eager model or diagnostic replay into production.
The added row scan is O(points), endpoint checks occur at most once per prepared
XY-surviving row/atom, and the original worst-case point×atom complexity remains.
A performance benefit depends on removing enough repeated cached tests to pay for
that scan and its finite/endpoint branches.

## Required candidate gates and reproduction

Before editing, freeze the accepted C1/C2 source, current profile, proposal, input
identities and diff in a new candidate directory. Preserve checked geometry APIs,
topology/center conventions and all scientific equations. Focused differential
tests must exercise actual endpoint rejection and retention, adjacent f32 boundary
values, signed zero, underflow, overflow, nonfinite fallback, arbitrary point/atom
order, row resets, dense first hits and lazy cache extension. Compare all 15 frame
fields against the unchanged direct predicate; keep the independent reference
kernels and established classifications unchanged.

Then run the complete corrected scientific observations, all private frame/bin
coverage, independent references, CLI/error/invariance cases and engineering
checks. Only after they pass may quiet paired AB/BA native trials run at 1/2/4/6
threads. Retain every sample, raw output and error; assess variation, CPU, RSS,
allocation diagnostics and all other workloads. Worker spans are not additive
wall shares. Reject or revise C3 if the target reduction is below 10%, within noise,
or accompanied by a material regression.

The [evidence package and reproduction instructions](../evidence/corrected_occupancy_shadow/REPRODUCTION.md)
retain exact executed controller/probe bytes, patch, source archive, receipts and
compressed per-frame counters. Original local archives remain unchanged. Rebuilding
or replaying diagnostics is separate from native performance measurement and must
not overlap a benchmark window.
