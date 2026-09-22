# Accepted C2: one Student-t multiplier per screen

This package retains the exact evidence for accepted C2, which reuses a lazily
computed Student-t multiplier within one screen operation. The accepted decision,
including measured tradeoffs, is [ACCEPTANCE.json](ACCEPTANCE.json). C2 preserves
the corrected baseline's scientific outputs and documented limits. C3 is outside
this package's admission scope.

The accepted build manifest SHA-256 is
`82b63f16c90e0767c5cb5b60daed61e2ce52f123d4f1a9f3156e791c5a03fd81`;
the native executable SHA-256 is
`b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`.
Acceptance SHA-256 is
`071f29b5e52cc4940b9581eb3de90430b3c1d7d4d2d1c4e7d4c9d30750f8cabe`.

## Retained evidence

- All **1,320 native launches**: 720 launches over the original baseline's full
  40 configurations, plus 600 incremental launches against accepted C1b over
  12 configurations. These include 608 measured pairs, 104 warmup processes,
  exact pair order, raw output, native resource usage and unfavorable samples.
- All 13 scientific lanes and 128,021 responses, public and private-bin
  conformer observations, complete independent references, original failed
  reference receipt and its explicit resolution, source/build and engineering logs.
- All 93 CLI cases at each of 1, 2, 4 and 6 threads, 24 additional CLI cases,
  and all eight million-prediction bitstreams from both baseline and C2.
- The pre-edit source and 46 uncached endpoint vectors, focused tests, independent
  source review, diagnostic cache counts and allocation records, including failed
  attempts. Counts are diagnostic observations, not native timing claims.
- The exact historical C1b test-path alias, its sealed pre-edit proof, original
  failed lookup, focused 30-test review and separate full prior-gate verification.
  The alias applies only to the old engineering record; current C2 verification
  remains strict. Original receipts are never rewritten.
- The complete accepted C2 reprofile: 40 configurations, 200 measured diagnostic
  samples and 40 diagnostic warmups, with output/resources/stage/allocation data.
  Its native denominators reference the retained native campaign.

[manifest.json](manifest.json) and [STORAGE_REPORT.json](STORAGE_REPORT.json)
record file counts, exact logical bytes, provider reuse, local compressed storage
and the largest contributors. Identical scientific output and prediction streams
reuse the immutable [C1](../accepted_c1/README.md) and
[C1b](../accepted_c1b/README.md) provider packages.

Native executables and disposable build/runtime caches are omitted with their
identities retained. External directory links are recorded without following or
restoring them. The first publication build refused a concurrent edit to an unexecuted C3
controller admitted by an overbroad file glob. Its plan/log/journal and
[failed-attempt receipt](FAILED_PUBLICATION_ATTEMPT.json) remain retained; the
final plan freezes only C2-executed controllers. No scientific or timing
results changed. No original output bytes, error strings or manifest paths are
normalized. Python sources are archived as `.py.txt` or exact compressed payloads.

## Verify and restore

Use a repository checkout containing C1, C1b, C2 and the original shared audit
and occupancy evidence. Restore the original audit payload if needed:

```sh
python3 docs/scientific_accuracy_audit/publication/manage.py restore
python3 docs/performance/evidence/accepted_c2/publish_shared.py.txt verify \
  --package docs/performance/evidence/accepted_c2 --shared-root "$PWD"
python3 docs/performance/evidence/accepted_c2/extended_coverage.py.txt \
  --package docs/performance/evidence/accepted_c2 --shared-root "$PWD" \
  --output "$(mktemp -d)/extended_coverage.json"
```

The full verifier checks every planned file, all local/shared file and chunk
hashes, provider metadata before and after reading, storage counts and the actual
1,320-launch/13-lane/four-thread CLI coverage contract. The additional checker
verifies the 24 CLI captures against their original baseline, full prediction
streams, actual cache-count records, historical alias-review bindings and all
240 reprofile samples. It supplements the full byte verifier. These commands
perform no scientific calculations or native benchmarks.

[RESTORATION.json](RESTORATION.json) records a complete verified restoration;
[EXTENDED_COVERAGE.json](EXTENDED_COVERAGE.json) records the additional checks.
Original and restored controllers must pass all 31 storage/provenance tests;
14 focused coverage tests exercise missing, duplicate and changed records.

Restore only into a fresh directory:

```sh
c2_restore_parent=$(mktemp -d)
python3 docs/performance/evidence/accepted_c2/publish_shared.py.txt restore \
  --package docs/performance/evidence/accepted_c2 --shared-root "$PWD" \
  --output "$c2_restore_parent/restored"
python3 -m unittest discover \
  -s "$c2_restore_parent/restored/.stericx/profiling/scientifically_validated_optimization/publication_c2_v1/controllers" \
  -p 'test_*.py'
python3 -m unittest discover \
  -s "$c2_restore_parent/restored/.stericx/profiling/scientifically_validated_optimization/publication_c2_v1" \
  -p 'test_*coverage.py'
```

The original workspace prefix is
`/run/media/andrej-rumenovski/New Volume/Code/StericX`. The index maps only
repository-relative paths into the new destination; provider paths resolve
under `--shared-root`. Sealed JSON retains its exact absolute strings, which
the restorer neither rewrites nor follows. Traversal, escaped symlinks, changed
provider metadata and mismatched hashes fail. Original evidence is never
replaced because the destination must not already exist.

## Exact assembly and fresh experiments

The restored specification and frozen controllers live under
`.stericx/profiling/scientifically_validated_optimization/publication_c2_v1/`.
The exact plan is published as `plan.json.gz`; assembly logs and readable helper
copies are retained beside this README. Assembly commands were:

```sh
python3 .stericx/profiling/scientifically_validated_optimization/publication_c2_v1/controllers/publish_shared.py plan \
  --repo "$PWD" \
  --spec .stericx/profiling/scientifically_validated_optimization/publication_c2_v1/spec_v2.json \
  --output .stericx/profiling/scientifically_validated_optimization/publication_c2_v1/plan_v2.json
python3 .stericx/profiling/scientifically_validated_optimization/publication_c2_v1/controllers/publish_shared.py build \
  --plan .stericx/profiling/scientifically_validated_optimization/publication_c2_v1/plan_v2.json \
  --output docs/performance/evidence/accepted_c2 \
  --journal .stericx/profiling/scientifically_validated_optimization/publication_c2_v1/build_journal_v2.jsonl \
  --pause-file .stericx/profiling/scientifically_validated_optimization/publication_c2_v1/PAUSE
```

Use new destinations for another assembly. A pause file stops an incomplete
build at a file/chunk boundary; `--resume` requires the unchanged plan and
verified journal/blob identities. Provider dependencies must match their pins
and form an acyclic graph. Later candidates may reuse C2 payloads through a new
package; storage reuse does not establish scientific or performance admission.

For a fresh experiment, follow the
[scientific reproduction recipe](../../../scientific_remediation/REPRODUCE.md)
with the captured C2 source and input hashes, then the
[paired benchmark protocol](../../../../scripts/benchmark_corrected_pair.md)
after scientific admission. Build in new scratch locations. New toolchains,
paths and timestamps can produce new metadata and binary hashes; retain the
mapping and receipts, and verify scientific values/bits under the same contracts.
Portable archive verification does not assert a new executable or timing replay.
