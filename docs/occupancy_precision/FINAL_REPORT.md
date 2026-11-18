# StericX occupancy and certified continuous volume

**Status: scientific implementation and validation complete; final controlled
performance measurements await a quiet-machine confirmation.** This report does not yet claim completion
of the campaign's performance objectives. All timing claims below are explicitly
classified. The final measurement controller refuses to proceed without the
frozen scientific and quiet-host gates. Its first final attempt failed before
any timing: 0.589 busy CPU cores exceeded the prospectively fixed 0.15 limit
(I/O wait 0.023 was below its 0.05 limit). GNOME Shell and Warp dominated the
background activity. No unrelated applications were stopped. See
[timing admission](timing_admission.json). Phase 12 explicitly calls for a
quiet-window confirmation once all non-timing work is complete.

## Scientific accuracy

The historical finite-grid calculation is preserved. Its f32 point lattice,
inclusive squared-distance predicate, operation order, atomic radii, scaling,
hydrogen policy, center/orientation conventions and regional reductions have
not been changed. Three accepted occupancy changes remove scalar tails with
overlapping AVX2 loads, prepare immutable row/Z/regional metadata, and retain
up to eight exact radius/density grid keys per thread under an 8 MiB buffer cap.
Overlapping SIMD results combine by idempotent OR. Short rows and unsupported
platforms retain scalar paths. No occupancy or molecular geometry is cached.

Continuous volume is a different mathematical quantity:

\[
V=\operatorname{vol}\left(B(c,R)\cap\bigcup_i B(c_i,r_i)\right).
\]

The original midpoint-chord estimate and conservative 3D box reference remain
available unchanged. The new explicit Rust API,
`sphere_union_volume_certified`, integrates certified vertical chord lengths
over adaptive XY rectangles. Conservative chord-range bounds apply everywhere.
Where interval predicates prove fixed overlap topology and endpoint ordering,
a tensor midpoint remainder bounded by hemisphere curvature supplies a tighter
enclosure. Analytical contained-sphere and lens cases avoid unnecessary
subdivision. The algorithm splits the rectangle with the largest uncertainty,
inherits conservative atomic candidates and intersects successive global bounds
to preserve nesting. See [the mathematical derivation](CONTINUOUS_MATHEMATICS.md).

`CertifiedVolumeResult` reports volume, lower/upper bounds, absolute and meaningful
relative error bounds, geometric evaluations, subdivisions, cells, runtime,
status and `target_met`. Budget exhaustion and floating-point resolution failure
are explicit. Contradictory bounds return an error. The certificate concerns the
exact supplied binary64 coordinates and radii under IEEE round-to-nearest with
gradual underflow; it does not bound experimental geometry or radius uncertainty.
`BuriedVolumeCalculator::compute_certified_with_center` applies the established
continuous scaling and hydrogen conventions with an explicit center. Certified
regional integrals are not implemented; finite-grid regional output is unchanged.

Two subsequent isolated improvements reduce the cost of this same certificate:
nonnegative interval products select two provably extremal products instead of
four, and per-integration scratch buffers replace repeated temporary allocations.
Signed/nonfinite intervals retain the general arithmetic fallback. Every bound,
work counter and stopping decision in the 695-case certificate replay is exactly
preserved. Scratch data never survives an integration or backs retained heap cells.

### Independent validation and tolerance outcomes

| Evidence | Result and scope |
|---|---|
| Original analytical enclosure suite | 167 observations, zero invalid enclosures; 18 stress targets remain unmet |
| Additional sphere-chain references | 96 analytical observations including rotations/permutations; rational arithmetic and integer square roots; pairwise inclusion–exclusion used only after proving a triangle-free overlap graph |
| Higher-order overlap references | 108 translated/coaxial cases, exact piecewise polynomial cross-section integration times a rational enclosure of pi; all contain the reference and reach 0.01 Å³ |
| Independent geometric predicates | 6,000 exact-rational rectangle/sphere cases, including extreme exponents and tangencies; zero failures |
| Positive-interval optimization | 100,000 deterministic bit-equivalence cases, including signed zero, overflow, underflow and fallback intervals |
| Original twelve molecular cases | All reach 0.1, 0.05 and 0.01 Å³; all reach 0.001 Å³ at the declared expanded budget |
| Broader molecular set | 220 cases across 56 real conformers, multiple radii/centers, hydrogen and radius-scaling policies; 212 reach 0.05 at 1,048,575 cells; eight R=5 cases fail that cap and succeed at the expanded 4,194,303-cell cap |
| Budget ladders | Nested certified intervals checked; no monotonicity claim for point estimates |

The old box method missed all twelve 0.05 Å³ targets at 262,143 boxes, with
reported error bounds 1.176–3.065 Å³. The new method reaches those same targets
using 30,609–141,719 XY cells, within the original numerical count cap. Cells and
boxes have different costs; their counts are not a runtime ratio. The original
twelve 0.01 cases use at most 504,455 cells; 0.001 uses 780,709–3,193,705.
The eight expanded R=5 cases require 1,079,129–1,332,145 cells. Acceptance
tolerances were never relaxed. Original failed attempts remain in the evidence.

Independent old 3D enclosures overlap all 220 broad molecular results. Eight
preselected rational octree checks also overlap. These are useful independent
consistency checks, not exact molecular truth values. The fine 0.001 molecular
reference used in accuracy plots comes from the same certified algorithm and is
labeled accordingly. Analytical references use independent equations.

The range-only chord ablation still misses all twelve 0.05 targets at 262,143
cells (errors 0.215–0.627 Å³). Smooth-region curvature bounds, rather than merely
changing dimensionality, are essential. All-atoms traversal preserves bounds but
requires 1.16–3.98 times as many geometric queries. The old reference's boundary
boxes cannot resolve partial volume tightly enough within its budget.

## Performance

Production baseline: `ffc698f479cc5b0d75ae4f6856edc81ad5ca9c2a`.
Final production implementation: `aeb3ca74` (full source hashes in the frozen
`final_profile_v2/build.json`). Portable release compiler settings and inputs are
identical between the final A/B builds. Historical six-thread throughput of
23,783.8 files/s is context only; it is not this campaign's measured denominator.

Controlled one/two/four/six/twelve-thread throughput, occupancy/regional/profile
latency, conformer throughput, CPU time, memory and accuracy-adjusted runtime
remain pending. The final protocol retains one warmup, eight alternating pairs,
full distributions and unsuccessful measurements. Fresh-process measurements use
warmed OS file caches; they do not measure cold-disk startup. No speed claim is
drawn from correctness preflights, simulator runtimes or concurrent validation.

Exploratory candidate results and unfavorable samples are retained in
`o1_measurements.json` through `o5_measurements.json` and `c3_measurements.json` /
`c4_measurements.json`. They justified further validation, not final controlled
speed claims. X-slab filtering regressed the representative path by about 12%.
General BVH/cell-list paths and internal parallelism did not amortize well on
small real molecules. Large sparse synthetic gains remain documented with their
scope; those methods are not silently selected in production.

Default prepared-grid buffer retention grows from 393,216 to 495,808 bytes
(26.09%). The multi-key cache can retain 8 MiB per thread plus fixed bookkeeping;
this is not a whole-process memory limit. Native allocation and memory diagnostics
are separate from native speed measurements. Scratch reuse reduces median allocation calls by 82.80% at 0.05 and 83.85% at
0.01 across the original twelve molecules. Peak live bytes at 0.01 change from
29,371,032 to 29,372,128; the adaptive heap dominates. Final accuracy/cost Pareto
results remain pending controlled timing.

## Profiling

Fresh baseline occupancy accounts for 115,933,006 of 130,723,468 modeled Callgrind
instructions (88.69%). The prior campaign's 81.55% is a different historical
measurement. Algorithmic counters and full-mask comparisons cover 42 initial and
101 spatial-research configurations, including 56 real conformers and synthetic
sizes through 1,024 atoms. Per-row candidate counts, tail decisions, duplicate
hits and preprocessing/execution costs are retained.

The first accepted-source reprofile reduces total instructions to 78,307,346
(40.10% fewer) and exclusive occupancy to 62,986,231 (45.67% fewer). Occupancy
remains 80.43% of modeled instructions. The final-v2 reprofile records 78,320,019 total instructions and the identical
62,986,231 exclusive occupancy instructions. [Before](baseline_flamegraph.svg) and [after](final_v2_flamegraph.svg)
flamegraphs show instruction counts, not CPU time. Assembly inspection finds no
FMA in the occupancy kernels. Static SIMD opcode counts cannot establish dynamic
hardware SIMD utilization.

Cachegrind's L1 data misses and branch counts fall, but its last-level data misses
increase and its branch misprediction fraction rises. These unfavorable findings
are preserved in [the profile](OCCUPANCY_PROFILE.md). `perf stat` is denied at
`perf_event_paranoid=4`. Native cycles, IPC, branch/cache events and memory
bandwidth are unavailable. The evidence identifies repeated candidate/predicate
and mask work; it does not prove a hardware arithmetic-port or bandwidth limit.

Continuous profiling identified interval arithmetic, chord evaluation and
allocation as substantial modeled costs. Arithmetic specialization precedes
scratch reuse in separate commits and paired experiments. Future work should
target singular/tangent boundary regions and proven topology changes; extremely
tight targets can still exhaust the cell budget or floating-point resolution.

## Competitive results

The scientific admission compares identical frozen inputs and common-grid point
decisions before timing. The 56-conformer gate validates 2,588,544 independent
exact-decimal/IEEE point decisions; the separately declared twelve-conformer,
three-radius profile validates 1,664,064. Thread-count preflights check ordered
outputs at 1/2/4/6/12 workers. Morfeus uses its public buried-volume calculation;
DBSTEP's occupied-point kernel receives the admitted common points.

DBSTEP's default grid has different semantics and is not substituted into this
comparison. The radius policy fixes R=2.5/3.5/4.5 Å and
`density=0.01*(R/3.5)^3`, keeping a 32-side lattice; it is distinct from a
constant-density profile. Matched finite-grid output is not proof of equivalent
continuous-volume accuracy. No equivalent certified continuous API was admitted
for the external tools, so no continuous speed-superiority claim is made.
Controlled warm/fresh-process competitive timing remains pending. Results will
apply to these workloads and versions, not universal tool rankings.

## Regression and reproducibility

All 128,021 historical outputs in thirteen lanes are exactly preserved on both
the first accepted-source replay and the fresh final-v2 replay. Full independent Morfeus/DBSTEP/Kraken reference,
regional-bin, model and thermodynamic replays preserve their scientific results,
including known historical disagreements. The final-v2 transition changes only
the separate certified implementation; fresh exact-output and engineering
receipts bind that transition. [Final scientific admission](final_admission.json)
records the completed gates and their hashes. This does not relabel historical discrepancies as
new passes or claim that every external reference program ran twice.

The final native and Rust 1.85 suites each pass 372 tests. Python bindings pass
71 tests and research checks pass 156 without skips. macOS/Windows target checks
pass, and QEMU AArch64 passes 369 applicable tests (three x86-specific tests are
not compiled on ARM). Emulated/compile-only checks are not physical-platform benchmarks.
The upstream DBSTEP suite retains 413 passes and three skips. The optimized
certificate replay preserves all 695 scientific responses exactly, excluding
runtime and allocator diagnostics only.

Accepted commits: `29bc3313` (AVX2 tails), `34e32d33` (prepared grid), `d732d1cb`
(bounded grid LRU), `15534fb8` (certified chords), `2ae0b3fb` (positive intervals),
`aeb3ca74` (scratch reuse). Freeze and validation evidence have separate commits.
No releases or publications were made.

[Reproduction commands](REPRODUCIBILITY.md), [baseline](BASELINE.md),
[candidate ledger](CANDIDATES.md), [failed experiments](EXPERIMENTS.md) and
[work log](WORK_LOG.md) identify the retained controllers and receipts. Raw
measurements, source archives, binary hashes and frozen datasets are under
`.stericx/occupancy_precision/`. A bare Git clone does not include those large
local historical datasets. Input hashes, compiler/hardware details and workload
identities are tracked in `baseline_manifest.json` and `workload_identity.json`.
No historical evidence or failed experiment was removed.
