> Historical evidence for audited `b515c4f`, superseded as the optimization baseline by corrected `1ecfbd5`. See [BASELINE_STATUS.md](BASELINE_STATUS.md). Do not use these timings or hypotheses as current-build measurements.

# Candidate admission against the freshly profiled audited build

These are **pre-implementation proposals**, not accepted optimizations or measured
speedups. [Baseline selection remains unresolved](BASELINE_STATUS.md). The numbers
below come from the fresh September 21 runs of `b515c4f`, not the September 16
profile. No production source has been changed in this task.

## 1. Reuse the screening model's Student-t multiplier

**Measured work:** five one-thread diagnostic launches attribute a median 50.757 ms
and 15.651% of process wall to `student_t_two_sided_quantile`, with 1,000 calls per
screen. Each call executes 200 bisection iterations, in addition to its bracketing
loop: 200,000 bisection iterations per 1,000-candidate workload. At six threads,
the same serial component represents 33.60% of diagnostic process wall.

**Theoretical ceiling:** eliminating this scope gives `1/(1-f)`, or 1.1855× at one
thread and 1.5059× at six. These are diagnostic bounds, not native predictions;
profiling overhead is larger at six threads. The native screen medians are
323.953 ms and 138.913 ms respectively.

**Expected realistic savings:** approximately 40–50 ms per complete screen on this
machine, contingent on native paired measurements. This suggests roughly 12–15%
less one-thread wall time. Retaining one evaluation reduces the quantile work from
1,000 calls to one (200 bisection iterations), with negligible per-candidate cache
access. These estimates have not been tested.

**Exactness argument:** `alpha=0.05`, `observations` and `parameters` are immutable
for a single loaded model. The current function is deterministic and has no
scientific side effects. A context bound to that same immutable model can compute
and reuse the identical `Option<f64>`. Initialize lazily after the existing
leverage guard; distinguish an uninitialized cache from a computed `None`.
Preserve `(multiplier * residual_standard_error) * sqrt(1 + max(leverage, 0))`,
the half-width finiteness test, endpoint arithmetic, and every existing error,
feature and candidate-validation step. The quantile algorithm and its documented
extreme-tail limitation must remain unchanged.

**Risk and complexity:** low arithmetic risk; small API/context change. A cache
must not be reused with another model, alter model serialization, or suppress
validation. Reusing leverage is a separate candidate and should not be bundled.

**Gate:** exact endpoint bits and complete screening outputs, including missing
geometry, invalid degrees of freedom, unavailable leverage, exclusions, bootstrap
uncertainty, applicability, ranking and exported decks. Independent statistical
residuals and classifications must stay fixed. Then require a reproducible native
paired improvement of at least 10% on this screen workload, exceeding noise,
without a meaningful regression on the other fixed workloads.

## 2. Eliminate additional buried-volume point–atom predicates

**Measured work:** `occupied_volumes` takes a median 4,161.406 ms and 81.753% of the
10,000-descriptor diagnostic process wall across 30,000 calls. Its ideal
elimination ceiling is 5.480×; halving this scope would give about 1.691× overall.
Native descriptor-batch median is 5,044.513 ms. The integration grid remains
unchanged. A separate, nonproduction shadow probe has now measured the current
lazy XY loop on all 10,000 input files, covering 30,000 frames:

| Operation | Count |
| --- | ---: |
| Current Z predicates | 1,943,834,272 |
| Predicates a row-endpoint rejection would remove | 969,971,505 |
| Added endpoint predicates | 86,487,986 |
| Net predicate reduction | 883,483,519 (45.4506%) |
| Additional XY preparations under the eager model | 6,947,929 |

The 56-conformer replay independently gives a 45.4513% net predicate reduction.
Every proposed skipped predicate was still executed and asserted to miss; all
15 per-frame scientific f32 values and both workloads' native stdout/stderr
matched exactly. Boundary, nonfinite, finite-overflow, arbitrary-order, atom-order
and empty-input fixtures also passed. Evidence, the isolated instrumentation
patch and hashes are under the new optimization root's `occupancy_probe/`.

This is a measured operation reduction, **not measured runtime improvement**.
If it translated proportionally to the entire 81.753% occupancy scope, the
descriptor workload would be about 1.591× faster. That is an optimistic model:
it ignores unchanged work and added scans, finite checks, branches and storage.
A realistic pretrial hypothesis is 10–25% less descriptor wall time, with a
minimum acceptance threshold of a repeatable 10% native improvement. A result
outside that estimate must be reported, including a failed speed gate.

Two exact approaches warrant that investigation:

1. For an atom outside a finite row's z range, evaluate the original occupancy
   predicate at the closest endpoint. If that endpoint misses, every point in
   the row misses. This can eliminate a candidate from all subsequent point tests.
2. For finite nondecreasing-z rows, find each atom's occupied index intervals by
   binary search of the original f32 predicate on either side of its z position.
   Combine intervals in a bitset and obtain the exact integer regional counts.
   Route unordered, short or nonfinite cases through the unchanged loop.

The predicate must remain `(dx*dx + dy*dy) + dz*dz <= radius_squared`, with its
original rounded intermediate values. Nonnegative additions and correctly rounded
subtraction/squaring are monotone on each side of the atom, providing the rejection
proof. Do not replace it with a square root, a subtracted-radius threshold, an
epsilon, different precision or reassociation. Duplicate points, signed zeros,
boundary plateaus and complete-row early hits need focused differential coverage.

**Unresolved performance risk:** default rows have at most 32 z coordinates, so
binary-search branches and bitset bookkeeping may cost more than the saved
arithmetic. Dense first-atom hits already perform well. Neither approach is yet
admitted for production implementation. Fresh operation counts must establish
the realistic savings before an isolated exactness-first experiment.

## 3. Process independent descriptor files concurrently

**Measured work:** `descriptors_command` processes files in a serial loop. Fresh
native `descriptors_10000` medians are 5,044.513 ms at one configured thread and
5,008.436 ms at six, using approximately 100% of one CPU in both cases. Inclusive
per-file descriptor work accounts for 99.2547% of one-thread diagnostic wall in
the batch and 87.8106% for the 56-file workload.

**Theoretical ceiling:** with parallel fraction `f` and `p` workers,
`speedup <= 1 / ((1-f) + f/p)`. The resulting 2/4/6-worker bounds are
1.985×/3.913×/5.784× for the batch and 1.783×/2.929×/3.728× for 56 files. These
conditional bounds ignore scheduling, frequency and shared-resource costs.
Single-file XYZ and SDF inputs receive no benefit from this candidate.

**Realistic hypothesis:** approximately 1.7–1.9×, 3.0–3.7× and 4.3–5.2× batch
speedup at 2/4/6 workers, subject to native paired measurements. This eliminates
no scientific operations; it distributes independent files across available
cores. CPU work may increase slightly through scheduling. Do not conflate
wall-time reduction with an algorithmic reduction in distance predicates.

**Exactness argument:** the full `descriptors_for_file` path already runs in
parallel database workers. Its molecule, grid and conformer state is local.
Collect outcomes with an indexed parallel iterator, then emit results and errors
serially in the original input order. Preserve duplicates, validation before
work, all failures and partial-success status, serializers and every within-file
floating-point reduction. Do not print errors from workers, short-circuit a
collection, sort inputs or change ensemble aggregation. Keep single-file and
one-thread serial paths to avoid unnecessary pool/buffer overhead.

**Memory and behavior risks:** six default integration grids retain 2.25 MiB
versus 0.375 MiB serially. Concurrent whole-file parsing, worker stacks, allocator
retention, result slots and error strings add further memory. The current native
batch RSS is 18.13 MiB and diagnostic peak requested heap is 12.43 MiB; concurrent
scratch does not simply add to the later JSON serialization peak. Measure real
RSS at each thread count, including larger SDF/finer-density cases. Fixed-input
error/output order must remain exact. Resource exhaustion is a particular risk
for unconstrained worker counts and large grids.

**Admission gate:** freeze shuffled mixed-success/failure and duplicate-input
batches, all-failure batches, all three output formats, and the multi-input donor
index precheck at every requested thread count. These cases are additional to
the successful full batch's existing scientific fingerprint. Require exact
scientific/behavioral equivalence, all independent references and engineering
checks before native timing. Require at least 20% reproducible batch improvement
at two or more workers, a much larger expected gain at six, and no meaningful
single-file or one-thread regression. Memory growth must remain reasonable for
the supplied inputs; reject or constrain the candidate if it does not.

This candidate is proposed, not implemented or accepted. Its large measured
parallel fraction makes it a priority once the baseline prerequisite is resolved.

## Lower-priority opportunities

Grid preparation now represents 9.134% of descriptor-batch wall (1.1005× ideal
ceiling), compared with 3.42% in the older profile. Its one-thread screening worker
time is 44.102 ms, or 13.733% of process wall. Worker spans overlap the loader wait
and must not be added to it. An exact, bounded grid cache could become worthwhile,
but synchronization, retained memory, configuration keys and validation/error
ordering require a separate admission record.

Parsing and Sterimol dominate ensemble packing, but repeated-file caching would
benefit this deliberately repeated corpus disproportionately. No representative
distinct-input performance evidence currently justifies that architectural change.

All accepted candidates must subsequently pass the complete scientific oracle,
independent-reference recheck and engineering checks before native A/B timing.
None of these proposals repairs or reclassifies the audit's negative findings.
