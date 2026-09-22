# Rejected C3: exact row-endpoint rejection

C3 preserved the corrected scientific outputs but was **rejected for slower
native performance**. [REJECTION.json](REJECTION.json) records the decision;
[RESTORATION_SOURCE.json](RESTORATION_SOURCE.json) binds the exact restoration to
accepted C2. All candidate evidence remains intact. C3 is absent from the final
production source and contributes no accepted speedup.

The rejected build manifest SHA-256 is
`bace96d478096c9db5db9985e33bc557a733c754119842a6d4b93dd39819b987`;
its native executable SHA-256 is
`91c7d3b4e1d53481be87b67c3574377b580d211154d4871174817710754a0856`.
Accepted C2's build is `82b63f16…`, native `b2335016…`. The complete identities
and independent interpretation are retained in
[PERFORMANCE_REVIEW.json](PERFORMANCE_REVIEW.json).

## Retained scope

- All **972 native launches**: 720 across the full 40 configurations and 252
  across the predeclared 14 incremental configurations; 432 measured pairs and
  108 warmup processes. Every unfavorable sample, process output and resource
  record remains. Both target batch medians were about 15.6% slower than C2.
- All 13 scientific lanes and 128,021 responses, 25 independent-reference stages,
  93 CLI cases at each of 1/2/4/6 threads, 24 supplemental CLI cases, and all eight
  million-prediction f32 streams. Scientific passes retain their scoped meaning.
- The separate 31-case command-boundary baseline/C3 oracle: 18 successes and
  13 exact validation errors, including complete direct buried-volume packed
  output and 115 conformer rows. Its independent review is retained.
- Pre-edit proposal/source, rejected patch, focused tests, source review,
  build/engineering logs, and actual 30,168-frame/452,520-field occupancy probe.
  Fewer counted predicates did not predict a native speedup; instrumentation
  counts are not timings or machine-instruction counts.
- Failed occupancy v1 fixture/test build, boundary v1 harness JSON expectation,
  original reference timestamp-mismatch gate and explicit metadata resolution.
  These failures are preserved without classifying them as scientific changes.

This package has no accepted C3 reprofile. Final fresh checks of the restored
accepted C2 executable belong to the separate final-optimization publication.
The current science/methodology/predictive limitations remain those documented
in the scientific remediation and final scientific report.

## Verify and restore

[manifest.json](manifest.json) pins the plan, exact restoration index and helper
bytes. [STORAGE_REPORT.json](STORAGE_REPORT.json) records local compressed bytes,
provider reuse and the ten largest new compressed contributors. Repeated payloads
reuse immutable [C1](../accepted_c1/README.md),
[C1b](../accepted_c1b/README.md) and [C2](../accepted_c2/README.md) providers.
No executable, disposable target, virtual environment or external directory link
is restored. Raw Python is archived as `.py.txt` or exact compressed payloads.

From the repository root, with those providers and shared audit files present:

```sh
python3 docs/scientific_accuracy_audit/publication/manage.py restore
python3 docs/performance/evidence/rejected_c3/publish_shared.py.txt verify \
  --package docs/performance/evidence/rejected_c3 --shared-root "$PWD"
c3_restore_parent=$(mktemp -d)
python3 docs/performance/evidence/rejected_c3/publish_shared.py.txt restore \
  --package docs/performance/evidence/rejected_c3 --shared-root "$PWD" \
  --output "$c3_restore_parent/restored"
python3 -m unittest discover \
  -s "$c3_restore_parent/restored/.stericx/profiling/scientifically_validated_optimization/publication_c3_v1/controllers" \
  -p 'test_*.py'
```

The full verifier checks every planned byte, all providers and chunks, exact
full40/incremental14 native coverage, all 13 observer lanes and four CLI captures,
and the actual rejected decision, independent reviews and restored C2 identity.
It requires 972 launches/432 pairs and rejects omitted or unrelated campaign
bindings. The 46 tests include unchanged storage/path tests and explicit rejection,
missing-scope, wrong-build and historical-source alias tests.
[RESTORATION.json](RESTORATION.json) records the actual fresh restoration.

Historical JSON strings keep the original workspace prefix
`/run/media/andrej-rumenovski/New Volume/Code/StericX`. Only repository-relative
restoration destinations are used. The old live BV source record is resolved to
its exact frozen C3 bytes `1df92d49…`; restored C2 resolves to `2215d273…`.
Those two aliases are constrained by the restoration proof, not arbitrary path
rewrites. Sealed receipts remain byte-identical. Portable verification does not
perform a new scientific calculation or native timing run.

## Assembly and fresh experiments

The captured specification and controllers restore under
`.stericx/profiling/scientifically_validated_optimization/publication_c3_v1/`.
The reviewed rejected-candidate controller retains storage format 2 and unchanged
provider handling; it requires `accepted=false` with the real rejection and
restoration receipts. It does not fabricate an acceptance-shaped proxy.

```sh
B=.stericx/profiling/scientifically_validated_optimization/publication_c3_v1
python3 "$B/controllers/publish_shared.py" plan \
  --repo "$PWD" --spec "$B/spec.json" --output "$B/new_plan.json"
python3 "$B/controllers/publish_shared.py" build \
  --plan "$B/new_plan.json" --output "$B/new_package" \
  --journal "$B/new_journal.jsonl" --pause-file "$B/PAUSE"
```

Use unused destinations. A pause file stops assembly at a file/chunk boundary;
resume requires the same plan and verified journal. Only actually executed root
controllers were frozen; no mutable future-controller directory was included.
The original rejection helper is preserved unchanged, including a harmless long
descriptive string; archived Python is not reformatted after execution.

For a fresh scientific experiment use the captured C3 source and the
[remediation reproduction recipe](../../../scientific_remediation/REPRODUCE.md),
then the [paired benchmark protocol](../../../../scripts/benchmark_corrected_pair.md)
if scientific admission is repeated. New toolchains/paths/timestamps produce new
receipts; scientific output and timing results must be evaluated anew. This
package preserves the completed negative performance result.
