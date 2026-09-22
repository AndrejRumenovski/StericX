# Accepted C1b and complete C1/C1b reprofiles

This is a lossless evidence package for the accepted C1b executable: ordered
descriptor file parallelism on retained commit `3ea936d`. The exact decision,
including parser/startup, CPU and memory costs, is [ACCEPTANCE.json](ACCEPTANCE.json).
It is separate from the earlier [accepted C1 package](../accepted_c1/README.md),
which remains unchanged. These receipts do not admit C2 or C3, transfer earlier
timings to a different executable, or declare all scientific limitations resolved.

The accepted C1b build manifest SHA-256 is
`def1aa57289814fd24e787dd95092cd0115a67091f250e4be9c49f0f4a5e663d`;
the native executable SHA-256 is
`e7e7f00a1b4802e219144535d5c6567eb5e97bd269d1073e931b880288c5023d`.
The frozen admission decision SHA-256 is
`a6a69381f79b3e82310c4b89b894ca4e54b331eaa554abd412155246eb8582cb`.

## Retained scope

- All 720 native launches: 40 workload/thread configurations, one warmup pair
  and eight measured pairs each. This retains all 320 measured pairs and 80
  warmup processes, raw streams, resource usage, pair order and unfavorable samples.
- All 13 scientific lanes and 128,021 responses, including complete public and
  private-bin conformer campaigns, full independent references and their reviewed
  limits. All source/build, helper and engineering evidence is retained.
- All 93 candidate CLI cases at each of 1, 2, 4 and 6 threads, complete native
  prediction vectors, 24 additional CLI bridge cases, and the retained failed
  21-case bridge attempt. Original failures, the documentation-only Ruff repair,
  reference metadata resolution, and both performance review attempts remain visible.
- Both accepted C1 and accepted C1b diagnostic reprofiles: each has all 40
  configurations, 200 measured diagnostic samples and 40 diagnostic warmups, with
  raw per-process output, resources and stage/allocation records. The original
  C1 package omitted its later reprofile; this additive publication supplies it
  without altering that package. Native denominator adapters refer to existing
  accepted paired samples; the reprofiles launch no additional native timings.

The archive represents **8,144 original files / 5,275,609,854 bytes**. It adds
23,278,364 compressed blob bytes, with a largest blob of 2,912,658 bytes.
There are 2,425 complete files supplied by the pinned C1 package, one inherited
chunk, and 3,871 new local chunks. Equal scientific observation streams and
prediction vectors are shared, not republished. [manifest.json](manifest.json)
and the compressed restoration index record every path, hash and dependency.

Native executables and disposable caches are omitted explicitly, with binary
identities retained in the receipts. External input directory links are recorded
but never followed or restored. The index preserves exact file bytes, including
unfavorable outcomes, error strings, timestamps and historical absolute paths.
All Python controllers are archived as `.py.txt` or inside exact-byte blobs.

## Verify and restore

Use a complete repository checkout containing both evidence packages and the
shared original audit/occupancy dependency files. Restore the original audit
payload first if it has not already been unpacked:

```sh
python3 docs/scientific_accuracy_audit/publication/manage.py restore
python3 docs/performance/evidence/accepted_c1b/publish_shared.py.txt verify \
  --package docs/performance/evidence/accepted_c1b --shared-root "$PWD"
python3 docs/performance/evidence/accepted_c1b/extended_coverage.py.txt \
  --package docs/performance/evidence/accepted_c1b --shared-root "$PWD" \
  --output "$(mktemp -d)/extended_coverage.json"
```

The main verifier checks complete original-plan membership, all local and inherited
file/chunk identities, pinned provider metadata before and after reading, source/
controller bindings, storage counts, and the actual 720-launch/13-lane/four-thread
CLI coverage contract. The additional checker verifies the retained failed and
successful CLI bridge cases and both complete 240-process diagnostic sample
matrices, including stream/stage/resource bindings. It does not substitute for
the main full-byte check. Neither command runs scientific calculations or timings.

The successful full verification and restoration receipt is
[RESTORATION.json](RESTORATION.json); the additional coverage result is
[EXTENDED_COVERAGE.json](EXTENDED_COVERAGE.json). Both original and restored copies
of the controller passed all 31 storage/provenance tests. Six additional tests
exercise missing/duplicate diagnostic samples, altered stage hashes and missing
or inconsistent CLI captures.

Restore only into a new directory:

```sh
c1b_restore_parent=$(mktemp -d)
python3 docs/performance/evidence/accepted_c1b/publish_shared.py.txt restore \
  --package docs/performance/evidence/accepted_c1b --shared-root "$PWD" \
  --output "$c1b_restore_parent/restored"
```

The original workspace prefix is
`/run/media/andrej-rumenovski/New Volume/Code/StericX`. Restoration maps that prefix
to the new destination using the index's repository-relative paths. Provider paths
are resolved beneath `--shared-root`. Sealed JSON content retains its original
absolute strings; the restorer never rewrites those strings or opens their paths.
Absolute restoration paths, traversal, escaped symlinks, mismatched hashes and
changed provider metadata are rejected. Recorded external links are not created.
The output directory must not already exist, so old `.stericx` evidence cannot be
silently replaced.

For restored tests, from the repository root:

```sh
python3 -m unittest discover \
  -s "$c1b_restore_parent/restored/.stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/controllers" \
  -p 'test_*.py'
c1b_test_dir=$(mktemp -d)
cp docs/performance/evidence/accepted_c1b/extended_coverage.py.txt "$c1b_test_dir/extended_coverage.py"
cp docs/performance/evidence/accepted_c1b/test_extended_coverage.py.txt "$c1b_test_dir/test_extended_coverage.py"
python3 "$c1b_test_dir/test_extended_coverage.py"
```

## Exact assembly and future providers

The original specification and controller/test snapshot restore under
`.stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/`.
The original plan is also retained directly as `plan.json.gz`; readable assembly
logs are published as `publication_plan.stdout.txt` and `publication_build.stdout.txt`.
The executed frozen controller was `controllers/publish_shared.py`; its readable
copy and all three required helper sources are published beside this README.
The actual assembly commands were:

```sh
python3 .stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/controllers/publish_shared.py plan \
  --repo "$PWD" \
  --spec .stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/spec.json \
  --output .stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/plan_v1.json
python3 .stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/controllers/publish_shared.py build \
  --plan .stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/plan_v1.json \
  --output docs/performance/evidence/accepted_c1b \
  --journal .stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/build_journal.jsonl \
  --pause-file .stericx/profiling/scientifically_validated_optimization/publication_c1b_v1/PAUSE
```

For another assembly, retain these originals and use new plan/output/journal
paths. A pause file stops an incomplete build between files/chunks; `--resume`
requires the same plan and rechecks retained journal/blob identities. New paths
and timestamps legitimately produce new metadata, while every restored original
file must keep its recorded bytes. This package's C1 provider registry pins the
unchanged C1 manifest, index and plan; no repository-wide search or implicit
provider substitution is allowed.

C2/C3 packages can register this package and C1 as immutable providers, then
reuse complete file identities or ordered content-addressed chunks. Dependencies
must be present with matching pins, must form an acyclic graph, and must remain
unchanged during verification. A storage dependency alone does not admit a new
scientific or performance candidate. [STORAGE_REPORT.json](STORAGE_REPORT.json)
records local/reused accounting and the largest new compressed contributors.

## Reexecute the experiment separately

Use the [scientific remediation reproduction recipe](../../../scientific_remediation/REPRODUCE.md)
for fresh builds, observations and independent references, followed by the
[paired benchmark protocol](../../../../scripts/benchmark_corrected_pair.md)
after full scientific admission. The restored native manifests preserve exact
workload arguments, environments and pair order; raw reprofile manifests preserve
their complete diagnostic commands and native denominator bindings.

Build in a fresh scratch tree using the captured C1b source, pinned dependencies
and original input hashes. Do not execute historical controllers against retained
output destinations. Relocation or another toolchain can change file paths,
model/build timestamps and binary hashes. Record the explicit mapping and new
receipts; verify unchanged scientific fields and bits under the existing contracts.
Portable archive verification does not claim a fresh executable or timing replay.
