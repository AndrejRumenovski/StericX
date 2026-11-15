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
