# Occupancy profile before optimization

Baseline source: `ffc698f4`; fresh builds and input identities are in
`baseline_manifest.json` and `.stericx/occupancy_precision/baseline/build.json`.
This report is being completed as measurements finish. Native timing/scaling
tables remain pending; the instruction results below are fresh measurements.

## Instruction and assembly findings

The public six-conformer observer performs one warmup, five iterations and one
output evaluation. Callgrind records **130,723,468 modeled instructions**, of
which **115,933,006 (88.69%)** are exclusive to AVX2 occupancy. Widths in the
[before flamegraph](baseline_flamegraph.svg) are instructions, not elapsed time.
Inline instructions stay attributed to their containing function. The sum of
self costs reconstructs the Callgrind total exactly.

Fresh Cachegrind records 130,748,513 instructions, 27,225,275 data references,
409,087 L1 data misses, 9,628 last-level data misses, and 927,617 branch
mispredictions in its configured simulation. These are **not hardware events**.
Its split 32 KiB/8-way L1 and 16 MiB/16-way last-level model omits actual L2,
speculative execution and hardware prefetch behavior. See
[event data](baseline_profile_events.json) and the
[Cachegrind documentation](https://valgrind.org/docs/manual/cg-manual.html).

The actual `perf stat` attempt was denied at `perf_event_paranoid=4`; cycles,
IPC, real branch/cache behavior and memory bandwidth are unavailable. No CPU
throughput or bandwidth ceiling is inferred from the simulator.

The disassembled occupancy function contains separate packed multiply/add,
ordered comparisons and POPCNT, with no FMA. It also contains permutations,
blends and variable shifts for gathering Z values and constructing positive-Z
masks. Static opcode counts are not dynamic SIMD utilization.

Callgrind's source attribution identifies 19.30 million instructions on XY
distance/filter source lines, 4.63 million on Z gathering/positive masks,
2.03 million on row detection, and 15.65 million on scalar-tail loop/distance/
mask lines. Optimized debug-line attribution is approximate: instruction motion
and intrinsic/inlined library attribution prevent treating these as disjoint
complete subroutine timing measurements. Raw assembly, annotated source and
profiler output remain under `.stericx/occupancy_precision/baseline_*`.

## Algorithmic counters

An isolated observer appends instrumentation to an exact source copy. A separate
replay of the AVX2 algorithm counts XY atom tests, retained candidates, full
eight-point vectors, scalar tail decisions, repeated covered hits, row lengths
and occupied points. It compares every vector mask against scalar operations
and every union bit against the independent original direct-distance loop.
All **42 configurations** (36 molecular orientations and six synthetic sizes
1–1,024 atoms) have zero mask mismatches. Counts are in
[baseline_counters.json](baseline_counters.json). Zero-iteration counter runs
have no timing interpretation.

The default grid has 15,408 points in 740 XY rows. Triphenylphosphine's first
frame executes 14,060 XY tests, retains 2,765 candidates and performs 59,660 Z
decisions, including 9,148 scalar-tail decisions. Its 4,517 occupied points
receive 5,224 redundant hits from overlapping atoms. For methylphosphine the
corresponding counts are 1,480 XY tests, 518 candidates, 13,452 Z decisions and
1,948 scalar-tail decisions. Preprocessing must therefore remain inexpensive
for small molecules.

## Ranked experiments

1. Remove scalar tails for rows of at least eight points using a final overlapping
   vector load. This preserves each predicate and can replace expensive scalar
   mask assembly with one vector operation sequence.
2. Reuse exact squared-X prefixes and filter candidates per X slab before XY
   rows. The implication uses nonnegative squared terms and preserves original
   f32 grouping; no approximate geometric cutoff is acceptable.
3. Cache immutable row/Z/region metadata along with the grid. Measure preparation,
   retained bytes and cold overhead separately. Geometry-specific occupancy
   cannot be cached under a grid-only key.
4. Compare interval/bitset unions and larger spatial structures for dense grids
   and large molecules. The previous 15.6%-slower endpoint-search experiment is
   retained historical evidence, not a predicted gain.
5. Evaluate internal parallelism only where a large individual integration
   amortizes scheduling. Batch parallelism already exists; nested pools must
   not be introduced implicitly.

Current evidence establishes substantial repeated predicate and mask work. It
does not establish whether hardware arithmetic ports, branch speculation or
memory bandwidth are the limiting microarchitectural resource. Candidate
timings and re-profiling will test the actionable hypotheses.
