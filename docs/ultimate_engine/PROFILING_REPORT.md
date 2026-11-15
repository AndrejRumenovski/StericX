# Engine profiling

Fresh baseline profiling is complete on Linux x86-64. Final candidate profiling
and controlled A/B acceptance remain pending. These are exploratory measurements,
not a final speedup claim. The baseline source is `c03630cb` and its frozen hashes
are in [baseline_manifest.json](baseline_manifest.json).

## Measurement scope and host

AMD Ryzen 5 5600G: six physical cores, twelve logical processors; Rust 1.97.0,
LLVM 22.1.6. Portable release builds use thin LTO and one codegen unit. Baseline
CLI used `RUSTFLAGS=-C debuginfo=1`; Cargo's release strip setting removed line
tables, so the separate observer uses `debug=1` in its release profile. Both use
the same frozen engine source. Assembly and symbol-build instructions are
retained with the builds. Native-target builds are not silently substituted.

175 CLI attempts cover seven scientific/I/O workloads at 1/2/4/6/12 threads:
one warmup, three warmed native launches and one instrumented launch each.
Every run has identical scientific fingerprints across threads/instrumentation.
71 kernel runs isolate computation, regions, frame construction, parsing,
statistics, serialization, radius profiles, density and in-memory batches.
Every mode's native and allocation-instrumented outputs agree. Output writing
and setup are outside the internal kernel timer; complete CLI measurements
include process startup, parsing and regular-file output. All attempts remain.

Full distributions, CPU time and RSS are in [baseline_measurements.json](baseline_measurements.json)
and [BASELINE_TABLES.md](BASELINE_TABLES.md). Raw commands, host pressure, binary
hashes, outputs, allocation snapshots and wait4 records are under
`.stericx/ultimate_engine/baseline_{cli,kernels,valgrind,python,python_layers}`.
The existing native wait4 launcher avoids attributing Python's heap to the child.
Background desktop processes remained running; no isolated-host claim is made.

## CPU and memory findings

The instrumented single-thread 10,000-file batch attributes 81.85% of its
5.107-second wall interval to occupancy, 9.08% to integration grids, 2.93% to
Sterimol support scans, 2.60% to XYZ parsing and 1.44% to file parsing/I/O.
Scopes on one thread are exclusive; overlapping worker elapsed times must not
be added to process wall time. Instrumentation overhead is assessed against
the paired workload's native distribution, not assumed zero.

That diagnostic requests 4,115,596,082 bytes through 3,056,774 allocations and
104,600 reallocations. Peak live requested bytes are 14,250,613. These are
allocation traffic and live heap, respectively; neither equals native RSS or
measures time spent in malloc. Separate DHAT on seven six-conformer evaluations
(warmup, five measured iterations, output evaluation) records 17,142,207 total
requested bytes, of which 16,515,072 are attributed to grid construction.
Allocation snapshots show zero net live growth for warmed kernel calls.

On the same ensemble workload, Cachegrind reports 297,778,376 modeled
instruction events, 39,936,664 conditional branches, 2,437,765 conditional
branch misses, 531,234 L1 data-cache misses and 12,891 last-level data misses.
These are **simulated**, not hardware measurements. Its explicit model is
32 KiB/8-way I1 and D1, 16 MiB/16-way last level, 64-byte lines; it omits the
host's separate L2, prefetching and speculative execution.

The [baseline instruction flamegraph](baseline_flamegraph.svg) comes from a
separate Callgrind run with 40-frame caller contexts, totaling 297,765,939
modeled instructions. The renderer sums self costs, excludes call-edge costs,
and verifies equality to Callgrind's total. Its widths are instruction counts,
not sampled CPU time; inline work remains attributed to the containing function.
The context-free run and full source attribution are retained separately.
See [Callgrind's documented semantics](https://valgrind.org/docs/manual/cl-manual.html).

Hardware cycles, retired instructions, IPC, native cache/branch misses and
memory bandwidth are unavailable: `perf stat` was actually denied at
`perf_event_paranoid=4`. No system permissions were changed. Samply and
cargo-flamegraph are absent; the retained Valgrind 3.26.0 installation supplies
call graphs, modeled events and DHAT instead. Memcheck reports zero errors and
zero lost bytes; 544 runtime bytes remain reachable at exit in this baseline
workload. This is a dynamic check of that workload, not a proof for all inputs.

## SIMD

[Disassembly counts](baseline_simd.json) show SSE packed instructions in frame
construction and the Sterimol scan, while the row-candidate occupancy inner
loop uses scalar subtraction/multiplication/addition and early-exit branches.
Packed instructions elsewhere in the containing function do not prove that
the occupancy predicates themselves are vectorized. There are no AVX packed
multiplications or fused operations in the inspected portable descriptor paths.
The existing explicit AVX2 implementation is in regression prediction, with
runtime feature detection and a portable fallback; it does not accelerate
buried-volume occupancy. Fallback and AVX reduction orders differ and must
not be advertised as bit-identical for arbitrary model inputs.

## Scaling and Python

For the 10,000-file workload, exploratory medians are 1,977.2, 3,745.4, 6,885.3,
10,264.5 and 12,360.9 files/s at 1/2/4/6/12 threads. Relative to one thread,
parallel efficiencies are approximately 100%, 94.7%, 87.1%, 86.5% and 52.1%.
The last value uses twelve logical threads sharing six physical cores. More
workers improve throughput here, with diminishing benefit. The full table
retains single-file and single-ensemble workloads that do not parallelize
internally; scheduling, I/O and compute are not conflated into one scaling claim.

Fresh local binding build: the six-conformer native descriptor/statistics call
has median 2.933 ms, the JSON bridge 3.274 ms and typed public API 3.973 ms.
Both transport layers exactly reproduce the native result. Portable native and
wheel profiles differ in metadata/build configuration, so their difference is
not a pure PyO3 transition cost. In the same Python process, independent lanes
measure JSON decode at 422 µs, typed-object construction at 176 µs, NumPy export
at 59.1 µs and options encoding at 5.91 µs. These isolated timings must not be
summed into a guaranteed end-to-end prediction.

The 56-file Python comparison preserves every typed conformer result. Median
per-file calls take 39.55 ms; batch calls take 39.54/29.72/17.10/13.84/13.24 ms
at 1/2/4/6/12 threads. These include input parsing and Python objects, unlike
the in-memory Rust kernel measurement. NumPy output values are checked exactly
against the result records, and conversion is timed without recalculation.

## Accuracy versus cost

Eight declared geometries × four resolutions are evaluated through the actual
public Rust API. Independent rational spherical-cap formulas with a rigorous
pi interval give analytical references; overlapping contained spheres use
inclusion–exclusion. All five measured calls per setting are retained. The
analytical values interpret the actual accepted f32 coordinates/radii, not
unavailable decimal originals. Regions are also captured.

The overlapping pair has absolute errors 0.8276, 0.05328, 0.08033 and 0.01486 Å³
as density decreases from 0.08 to 0.00015625 Å³. The nonmonotonic step remains.
A near-tangent intersection has continuous volume about 3.31e-10 Å³, but one
resolution reports 0.0013103 Å³. Finite-grid identity and continuous accuracy are
therefore separate acceptance requirements. No target accuracy is inferred
solely from density or a successive-grid difference.

## Ranked opportunities, before production edits

1. **Occupancy:** highest CPU cost; investigate rowwise union masks, explicit
   SIMD and interval searches with the exact existing f32 predicate. Retain a
   portable fallback. The historical row-endpoint candidate was slower despite
   fewer predicates; its failure is not overwritten.
2. **Grid reuse/layout:** about 9% CPU and most allocation traffic. Investigate
   bounded reusable grids and row metadata while preserving every grid node and
   boundary decision. Account for retained per-thread memory.
3. **Density parsing:** Ne cube parsing takes 8.509 ms and 140,636 allocations,
   versus 177.7 µs for density integration. The source unconditionally allocates
   a replacement string for each numeric token. Test a borrowed ordinary-token
   path while retaining identical Fortran D/d parsing and error behavior.
4. **Python transport:** measured JSON and typed-object construction are
   meaningful for small molecules. Investigate avoiding repeated conversions
   and use existing batches without moving chemistry to Python.
5. **Sterimol and general parsing:** 13.70 µs per measured Sterimol call and
   60.17 µs per six-conformer parse. Investigate reusable scan directions and
   allocation reduction after the larger occupancy opportunities, with exact
   angle arithmetic and unchanged input/error conventions.

All opportunities require isolated changes, scientific validation and native
A/B evidence. Assembly appearance or modeled operation counts alone cannot
admit an optimization. Actual final speedups, final profiles, compiler
vectorization remarks and the final evidence seal will be appended after work.
