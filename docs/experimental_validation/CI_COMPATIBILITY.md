# GitHub CI compatibility follow-up

The first push at `7d9e56da` passed all fifteen platform/package jobs but failed
the Linux quality job under Rust 1.100. The campaign's local Rust 1.97 check had
passed. New Clippy diagnostics identified two indexed loops in regional integration;
running Rust 1.100 locally also identified one in the regional tests.

The follow-up uses array iterators in those three loops. Iteration order, arithmetic,
assertions, numerical tolerances and defaults are unchanged. Rust 1.100 formatting,
Clippy with warnings denied, documentation with warnings denied, and all 385 Rust
tests pass locally.

A fresh source-bound replay also reproduces the archived geometric, regional,
radius-profile, published-conformer, generated-conformer and statistical results
exactly, excluding only the controller's named timing/allocation fields. This replay
uses the campaign's Rust 1.97 toolchain. Commands, source hashes, row counts and
results are retained in [ci_compatibility_replay.json](ci_compatibility_replay.json).
Reproduce it with `reproduce_results.py --output` pointing to a new directory under
`.stericx/experimental_validation`, as described in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

The original campaign manifest and reports remain unchanged historical evidence
for the source at campaign completion. This follow-up records the subsequent source
changes separately. CI verifies the pinned original manifest's 83 published files
and syntax of its 23 archived Python helpers. Exact-name Ruff exclusions preserve
those recorded bytes; new helpers inherit repository rules. The archive verifier's
35 tests cover tampering, path traversal, missing files and exclusion scope.
