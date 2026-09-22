# Final accepted optimization: C1/C1b and C2

The final decision [ACCEPTANCE.json](ACCEPTANCE.json) retains ordered per-file
descriptor work and per-screen Student-t multiplier reuse. The accepted build
manifest SHA-256 is
`82b63f16c90e0767c5cb5b60daed61e2ce52f123d4f1a9f3156e791c5a03fd81`;
the measured and delivered native executable SHA-256 is
`b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`.
C3 was rejected for slower performance; its complete evidence remains in the
[separate rejected package](../rejected_c3/README.md).

This package preserves the final re-admission of the unchanged C2 executable.
[PERFORMANCE_REVIEW.json](PERFORMANCE_REVIEW.json) contains every configuration,
including regressions and outliers. The [final scientific report](FINAL_SCIENCE.md)
keeps optimization equality, independent implementation checks, methodological
limitations and predictive validity separate. `scientific_global_pass` remains
false in the independent-reference gate; finite exactness is not a universal
scientific guarantee.

## Retained scope

- All **720 final native launches**: 40 configurations, 320 measured pairs and
  80 warmup processes; ten workloads at 1/2/4/6 threads. No sample was excluded.
- Fresh post-performance replay: all 13 scientific lanes and **128,021 responses**,
  including all 31,721 conformers in each inferred/topology corpus and private
  stream, complete frames/bins, and the eight retained inferred-connectivity errors.
- All 25 independent-reference commands, 7,624 model comparisons, full geometry
  and thermodynamic reviews, original six-timestamp failure and the unchanged
  field-specific metadata resolution. No scientific tolerance was widened.
- Four original 93-case CLI captures, the additional 24 cases, and the final
  **31-case boundary replay**: 18 successes and 13 exact validation errors.
  Boundary raw fingerprints are recomputed under the frozen named process-metric
  policy, including complete direct-volume artifacts and DB/predict outputs.
  The separate eight million-value f32 prediction streams remain retained from
  the same accepted build; command previews do not replace those full vectors.
- Final Rust/Python/Clippy/rustdoc/format/lint checks, integrity tests, Study 011
  verification, exact rebuilt-release identity, helper reviews and documentation
  checks. Publication tool failures and their later passing checks are preserved.

## Verify and restore

[manifest.json](manifest.json) pins the exact plan, restoration index and helper
bytes. [STORAGE_REPORT.json](STORAGE_REPORT.json) lists provider reuse, local
compressed size and the ten largest new contributors. Immutable
[C1](../accepted_c1/README.md), [C1b](../accepted_c1b/README.md),
[C2](../accepted_c2/README.md) and [rejected C3](../rejected_c3/README.md) supply
repeated bytes. The rejected package is only a byte provider; it cannot supply
final acceptance.

From the repository root, with those packages and shared audit files present:

```sh
python3 docs/scientific_accuracy_audit/publication/manage.py restore
python3 docs/performance/evidence/final_optimization/publish_shared.py.txt verify \
  --package docs/performance/evidence/final_optimization --shared-root "$PWD"
final_restore_parent=$(mktemp -d)
python3 docs/performance/evidence/final_optimization/publish_shared.py.txt restore \
  --package docs/performance/evidence/final_optimization --shared-root "$PWD" \
  --output "$final_restore_parent/restored"
python3 -m unittest discover \
  -s "$final_restore_parent/restored/.stericx/profiling/scientifically_validated_optimization/publication_final_v1/controllers" \
  -p 'test_*.py'
```

[VERIFICATION.json](VERIFICATION.json) and [RESTORATION.json](RESTORATION.json)
record the completed full verification and fresh restoration. The 53 passing
tests comprise 31 storage assertions and 22 scope checks using a small actual
receipt fixture. The storage-test setup additionally copies the new final-scope
helper; original storage assertions are unchanged. The fixture tests reject
missing cases/commands, wrong builds, altered raw boundary outputs, incorrect
executed native/context, widened model tolerances and a different release binary.
An independent review is preserved alongside its own test logs.

Historical receipts retain the original absolute workspace prefix
`/run/media/andrej-rumenovski/New Volume/Code/StericX`. Restoration destinations
are mapped only by safe repository-relative paths; sealed JSON strings are not
rewritten. Native executables, disposable targets, virtual environments and
external links are not restored. Source snapshots, build/input/helper identities
and reproduction instructions remain. Raw Python is stored as `.py.txt` or exact
compressed payloads, preserving its executed bytes.

## Reproduce assembly or run a new experiment

The captured specification and controllers restore under
`.stericx/profiling/scientifically_validated_optimization/publication_final_v1/`.
The controller preserves shared storage format 2 and adds only final coverage
checks to the accepted-candidate contract.

```sh
B=.stericx/profiling/scientifically_validated_optimization/publication_final_v1
python3 "$B/controllers/publish_shared.py" plan \
  --repo "$PWD" --spec "$B/spec.json" --output "$B/new_plan.json"
python3 "$B/controllers/publish_shared.py" build \
  --plan "$B/new_plan.json" --output "$B/new_package" \
  --journal "$B/new_journal.jsonl" --pause-file "$B/PAUSE"
```

Use unused destinations. The pause file stops at a file/chunk boundary; resume
requires the same plan and a verified journal. Only actually executed final root
controllers are frozen; no mutable future-controller directory is swept in.

Portable verification reconstructs the completed evidence without a scientific
or performance rerun. A fresh experiment follows the
[remediation recipe](../../../scientific_remediation/REPRODUCE.md), frozen final
observer/reference helpers and [paired native protocol](../../../../scripts/benchmark_corrected_pair.md).
New paths, toolchains and timestamps produce new receipts; outputs, locked
scientific classifications and all measured samples must be evaluated again.
The additive final report and archive receipts were written after the immutable
core finished; its earlier pre-publication report remains preserved in the core.
