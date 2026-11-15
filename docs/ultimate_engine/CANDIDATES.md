# Prospective candidate ledger

Each entry predates its corresponding production edit. Baseline source/archive,
measurements and raw profiler results are immutable under `.stericx/ultimate_engine`.

## U1: borrowed density numeric tokens

Measured target: Ne cube parsing, 8.509 ms and 140,636 allocations/call. Source
`density::number` creates a replacement String for every token even without D/d.
Borrow ordinary tokens; retain the old replacement exactly for Fortran exponents.
Expected benefit: remove most numeric-token allocations; native parsing should
improve. No descriptor arithmetic, constants, thresholds or node ordering changes.
Admission: exact parsed metadata/values and errors, existing density tests and
independent DBSTEP comparisons, memory check, alternating native timing blocks
and allocation deltas. A slower candidate is rejected; no tolerance changes.

## Later occupancy/grid hypotheses (not yet implementations)

Fresh profiling supports bounded grid reuse and SIMD occupancy over contiguous
Z rows. IEEE grouping must remain `(dx*dx + dy*dy) + dz*dz <= r*r`; algebraic
subtraction of the squared radius is not an exact replacement near absorption.
No early atom rejection may remove an originally occupied point. Candidate
grids, masks, regional counts and every public result/error must be compared
against frozen scalar behavior and independent common-domain occupancy.

## U1 experiment outcome

Retain for full final replay: 23 complete parsed-grid/error responses are byte
identical, including three official cubes and D/d, subnormal, underflow,
nonfinite and malformed numeric tokens. Eight alternating native AB/BA pairs
on CPU 2: medians 7.697106 ms baseline and 5.744153 ms candidate (1.3400×,
25.37% less time). Every sample is in [u1_measurements.json](u1_measurements.json).
This is an exploratory desktop result, not the final controlled claim.

Allocations per Ne cube fall from 140,636 to 12; requested bytes from 5,742,068
to 4,195,248. Reallocations remain 18 and peak live memory is effectively
unchanged: this removes temporary allocation traffic, not the retained density
array. All eight focused Rust tests and all 12 pinned DBSTEP density/reference
cases pass. Memcheck reports zero errors and zero lost bytes; 544 runtime bytes
remain reachable. Raw outputs, compressed full grids, binary hashes, reference
receipt and logs are under `.stericx/ultimate_engine/u1_validation` and
`u1_density_reference`. `cargo fmt` subsequently joined one Rust expression;
the experiment snapshot preserves the pre-format source. Final replay uses the
formatted committed source.

Reproduction: build the documented observer against each frozen source, then
run `python docs/ultimate_engine/check_density_candidate.py` in a fresh output
location (the controller deliberately refuses to overwrite prior evidence).

## S1 scientific witness: scale-dependent R²

Both fit metrics and revealed-target evaluation use an absolute TSS threshold.
For targets ±s and predictions ±s/2, exact arithmetic gives R² = 3/4 for any
nonzero s. Dyadic s = 2^-30 triggers `None` in the baseline despite representable,
nonconstant input. A regression test is run against the original metric code
before correction; it also checks constant targets and a legitimate negative
score. Preserve `.stericx/ultimate_engine/r2_baseline_failure.log`.

## U2: bounded per-thread integration-grid reuse (pre-edit)

Measured target: grid construction consumes 9.08% of the diagnostic 10,000-file
run; grid allocations account for 96.34% of DHAT requested bytes in the ensemble
workload. The grid depends only on the exact f32 radius and density bits.
Keep one grid per thread, keyed by those bits, shared locally through Rc. Retain
at most 8 MiB of grid capacity per thread; larger grids use the unchanged builder
without retention. No global lock, approximate key, node recomputation or public
grid definition change. This trades bounded retained memory for fewer allocations.
Preserve the uncached builder for independent observers. Check exact grid/region
results, key switching, thread equivalence and finite memory retention. Benchmark
warm kernels and cold/full CLI workloads; reject if the tradeoff is not useful.

### U2 experiment outcome

Retain for final replay. Eight alternating pairs show speedups of
1.1390× for volume,
1.1319× for regions and
1.1306× for six-conformer ensembles.
The changing-radius scan is essentially unchanged (1.0027×);
it repeatedly replaces the one-entry cache. Sterimol is an unchanged control.
All native results are exact against the frozen baseline. Full samples are in
[u2_measurements.json](u2_measurements.json). These are exploratory measurements.

Requested bytes per warmed six-conformer call fall from 2,415,541 to
55,927; allocation calls from 620 to 614.
The default retained grid capacity is 393,216 bytes per participating thread,
with an 8 MiB maximum cached capacity per thread. The tradeoff is persistent
bounded capacity instead of repeated temporary allocations. No warm live-byte
increase or allocation failures; Memcheck reports zero errors/lost bytes.

All 352 Rust tests pass. The fresh independently bound DBSTEP point-domain
replay has zero disputed/unexplained atomic points, zero unexplained sphere
points and zero arithmetic errors; raw evidence is at
`.stericx/ultimate_engine/u2_point_domains`. Source integrity checks pass.
Final full CLI/cold-start A/B and complete scientific replay remain mandatory.

## U3: AVX2 row occupancy masks (pre-edit)

Measured target: atomic occupancy was 81.85% of the pre-cache diagnostic CPU
time. The source and disassembly show scalar branch-heavy point/candidate loops.
Try gathering at most 64 contiguous-row Z coordinates onto the stack, testing
eight points per AVX2 vector, OR-ing atomic hit masks and counting integer region
occupancies. Reuse `(dx*dx + dy*dy)` exactly; compute `xy + dz*dz <= r²` with
separate multiply/add and ordered comparison. No FMA, radius subtraction or
approximate geometric rejection. Stop only when all row points are occupied.
Retain the original scalar kernel for unsupported CPUs. Dispatch requires AVX2
and POPCNT support; unaligned loads must stay inside the slice. The SIMD path
must be proved active in assembly and dynamically observed, with scalar/SIMD
adversarial masks, complete finite-grid/regional comparison, independent DBSTEP
replay, memory checking and native A/B distributions before acceptance.

Compiler/intrinsics references: [LLVM vectorizer diagnostics](https://llvm.org/docs/Vectorizers.html)
and [Rust x86_64 intrinsics](https://doc.rust-lang.org/stable/core/arch/x86_64/index.html).

### U3 experiment outcome

Retain for full replay. Eight alternating pairs versus U2: volume 2.5418×,
regions 2.5512×, ensemble 2.2838×, radius profile 2.1141×. Combined versus frozen
baseline: volume 2.9295×, ensemble 2.5332×. The unchanged Sterimol control differs
by about 1%; do not call this a Sterimol improvement. All result JSON is exact.
See [u3_vs_u2.json](u3_vs_u2.json) and [u3_vs_baseline.json](u3_vs_baseline.json).
These are exploratory desktop kernel measurements, not final files/s claims.

All 353 Rust tests pass; the added tail/boundary test compares 5,200 vector masks
(lengths 0–64, signed zeros, subnormals, absorption, adjacent floats and nonfinite
raw predicates). The existing independent direct-distance tests now also force
the original scalar fallback, including random non-grid ordering and all real
fixture orientations. The fallback test assertion was added after the timing
snapshot; it changes only test code. Twenty-three focused geometry tests pass.

The 69-case/207-frame DBSTEP replay again has zero disputed atomic points and
zero unexplained/arithmetic differences; historical second-clip failures remain
preserved. [u3_point_domains_summary.json](u3_point_domains_summary.json) records
this scope. Memcheck reports zero errors/lost bytes, 544 runtime bytes reachable.

Assembly contains AVX packed subtraction, multiplication, addition, ordered
comparison and POPCNT, without FMA in the occupancy function. Callgrind includes
executed `occupied_volumes_avx2` instructions, proving dispatch is active under
the profiler. See [u3_simd.json](u3_simd.json); raw assembly and all four profiler
runs are retained under `.stericx/ultimate_engine/u3*`. Full final scientific,
thread-scaling, cold-start and controlled CLI A/B replay remain pending.

## N1: explicit continuous angular Sterimol

The pre-implementation definition/proof is in NUMERICAL_ACCURACY.md. Retain
this explicit accuracy option: 39 independent 80-digit cases, 39 matched
Morfeus/DBSTEP cases, analytical and rigid-transformation tests and Memcheck
pass. Preserve the two failed convention-mismatched external lanes. This does
not replace the default scan; the cubic algorithm is slower for some larger
molecules and its numerical error is not certified. Full final replay pending.

## Compatibility review and S2 outcome

Clippy caught use of adjacent-float APIs newer than the manifest's declared
Rust 1.85 minimum in the uncommitted volume prototype. They were replaced with
explicit IEEE bit stepping, without changing bounds. The acceleration's new
cache lookup was rewritten without a let-chain, and boundary-test adjacent
floats now use bit construction. Existing baseline let-chains elsewhere were
also discovered. The real Rust 1.85 baseline failure is preserved. Six existing
let-chains were rewritten equivalently and intrinsic unsafe blocks made
compatible; no minimum version was raised. All 364 tests pass with actual
Rust 1.85 and 1.97. All 362 applicable ARM tests pass with Rust 1.85-generated
AArch64/musl binaries under QEMU; this is functional evidence, not hardware timing.

## N2: continuous-volume quadrature and conservative reference

Retain the two explicit methods with the documented limitations. Independent
analytic/rational enclosures, 10,000 integer box checks, exact permutation tests,
native/instrumented comparisons and Memcheck pass. General-box strict targets
are often unmet and reported as such; no accuracy guarantee is inferred from
midpoint resolution or convergence differences. The original finite-grid
method remains the default. Initial invalid stress input and all unsuccessful
target attempts are preserved. See NUMERICAL_ACCURACY.md for the proof and scope.

## Final replay status

The prospective entries and exploratory measurements above retain their original
chronology. Final source `f9c952e2` reproduces all 128,021 original responses
exactly and passes the full engineering and explicit-method replay described in
VALIDATION_REPORT.md. Final controlled timing remains a separate admission step.
No fourth production optimization is introduced during the final replay.

The historical C3 row-endpoint filter remains rejected (about 15.6% slower).
Uncertified midpoint convergence extrapolation, sampled Sterimol as an exact
continuous reference, and a global replacement of the finite lattice were not
adopted. General adaptive-box integration remains an explicit reference option:
its rigorous enclosure is useful even when its requested tolerance is not met,
but the measured molecular budgets do not justify making it the default.
