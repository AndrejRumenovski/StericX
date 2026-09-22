# C1b: carry the newer CLI/CI commit into the optimization sequence

While optimization was paused, commit
`3ea936d553940dc51a39550cc2cc77ed6be54018` added standalone scientific test
fixtures and moved three large clap argument groups into separate derived
parsers. That existing work is retained. C1b is the current commit plus the
already implemented C1 scheduling change; it is a compatibility bridge, not a
new scientific algorithm or an independently claimed optimization.

The old accepted C1 executable remains immutable. Its measurements cannot be
silently assigned to the changed executable. All future cumulative performance
comparisons still use the scientifically corrected `1ecfbd5` baseline; incremental
C2 comparisons must use the actual C1b executable to avoid attributing CLI/layout
effects to quantile reuse.

## Static and behavioral evidence

The independent source review found only `src/cli.rs` and `src/main.rs` changed
between the 84 frozen C1/C1b Rust/Python entries. All 67 moved argument fields
(19 fit, 28 search, 20 screen), attributes, defaults, comments and order match
after visibility/whitespace normalization. Dispatch matches after the three
destructuring syntax changes. The 15 added test payloads match the already sealed
audit witnesses byte for byte; the tests strengthen error assertions.

Generated clap parser/argument-group code nevertheless changes, so static
similarity does not establish runtime compatibility. An additive 24-case native
comparison covers help, missing arguments, conflicts, aliases, invalid values,
the full 1,543-record search workload, search aliases and filtered constraint-only
search. All exit codes, stdout and stderr match exactly, without normalization.
The ordinary 93-case CLI oracle has no search cases, so these are meaningful
additional checks rather than an assumption about its coverage.

The first helper attempt compared a baseline executable named `native` with one
named `stericx`. Eleven usage messages differed in that executable name. Both
raw attempts are retained. The successful replay executes the same baseline
bytes at their original `stericx` path; it does not normalize away output text.
This was a harness invocation mismatch, not evidence of a scientific change.

## Captures and gate status

Paths below are under
`.stericx/profiling/scientifically_validated_optimization/`:

- Source/fixture review: `candidate_c1b_static_review_v1/manifest.json`, SHA-256
  `10b8bab1fb0af8ae8e2aa8fc4e01d82f917939e4cd556c467c2d1f2e3e279dd6`.
- Additional CLI comparisons: `candidate_c1b_cli_group_checks_v2/manifest.json`;
  failed invocation retained at the corresponding `v1` root.
- Frozen build: `candidate_c1b_science_v1/validated_build_corrected_vC1b/manifest.json`.
  Native SHA-256:
  `e7e7f00a1b4802e219144535d5c6567eb5e97bd269d1073e931b880288c5023d`.
- Complete exact observation comparison:
  `candidate_c1b_science_v1/exact_comparison/comparison.json` passes all 128,021
  responses with tolerance zero.
- Diagnostic/export build: `candidate_c1b_diagnostics_v1/builds.json`;
  all eight full million-value native/diagnostic prediction streams match the
  original `4afd2f0b…1d97` baseline digest. The exporter-native executable is
  byte-identical to this C1b build.

Full CLI comparisons, independent reference replay and engineering checks pass.
The reference replay retains the same six timestamp-derived model hash differences
as C1; the additive resolution proves identical model content after the previously
declared creation-time exclusion. No scientific value or tolerance changes.

The engineering controller initially failed its final Ruff formatting check on a
Python example in the occupancy evidence documentation. The original document and
publication index are preserved under that package's `publication_history_v1/`.
Formatting the example, updating its publication index and rerunning the package
verifier and both Ruff checks completed the gate. All 308 Rust tests, 124 Python
tests, Clippy, rustdoc and Rust formatting had already passed on unchanged source.
Four additional reference-integrity tests and the study verification also pass.

The full pre-timing admission is
`candidate_c1b_performance_gate_v1/manifest.json`. The 40-configuration native
comparison completed all 720 launches. The independent review verifies all raw
results and recommends acceptance with disclosed costs:
`candidate_c1b_independent_review_v2/review.json`, SHA-256
`a3693a581c003b63ba897162fd53e60732e8d3dfba3ec841020d55129c5c24d3`.

The bridge is accepted in `accepted_c1b_v1/acceptance.json`, SHA-256
`a6a69381f79b3e82310c4b89b894ca4e54b331eaa554abd412155246eb8582cb`.
Paired median speedups for 10,000 files are 1.937×/3.732×/5.422× at 2/4/6
threads; for 56 files they are 1.817×/3.063×/4.006×. Every measured pair in these
six configurations improves. One-thread large batches remain approximately
neutral. These compare this executable directly with the corrected baseline.

Costs remain visible: roughly 1% parser regressions in several settings, about
0.11 ms slower single-small launches in two settings, more total CPU work in
parallel batches, and batch peak RSS of 26.57 MiB. Screening and four-thread
database-build wall penalties are retained; near-neutral CPU does not establish
that wall differences are merely noise. Search is too dispersed for a precise
wall claim and has a small two-thread CPU cost. These are accepted bounded
tradeoffs for the large batch improvement, not claims that all workloads improve.

A fresh profile of this accepted executable is required before C2 changes source.
