# Accepted C1: lossless evidence publication

This package retains the evidence for **ordered descriptor file parallelism on
the frozen corrected baseline**. The exact decision is readable in
[ACCEPTANCE.json](ACCEPTANCE.json). It accepts large batch improvements with the
measured parser, screen CPU, search CPU and memory costs stated there. It does
not transfer those measurements to a later source revision or claim that all
scientific limitations have disappeared.

The accepted candidate build manifest has SHA-256
`13b99d5157a57fbcc2372e8f25d09841160db0a75396fc507bf8edd6f3309221`.
Its native executable has SHA-256
`2f65e0aba46a502bb653c61ac3090fa0acbb0244dedcf4aba18ce06e458be358`.
The subsequent HEAD `3ea936d` bridge and any C2/C3 candidate require separate
receipts. This package contains no newly measured timings or scientific results.

## Coverage and storage

The restoration index retains 6,479 regular files, totaling 4,566,739,926 bytes
before deduplication and compression. It includes:

- All 970 native process launches: 720 in the original 40-configuration campaign
  and 250 in the five-configuration supplement. All 440 measured pairs and 90
  warmup launches retain their metrics, stdout, stderr and native resource usage.
- All 13 scientific observation lanes, totaling 128,021 rows, with original
  public observations and private frames/bins; all independent reference results,
  including full volume/bin and fixed historical variant coverage.
- All 93 CLI cases at each of 1, 2, 4 and 6 threads, their raw artifacts, exact
  comparisons, and engineering logs.
- Complete native prediction vectors, build and source snapshots, exact executed
  helpers, the original failed build, and the failed reference gate followed by
  its reviewed timestamp-only resolution. Both reference decisions remain
  readable as [FAILED_REFERENCE_GATE.json](FAILED_REFERENCE_GATE.json) and
  [RESOLVED_REFERENCE_GATE.json](RESOLVED_REFERENCE_GATE.json).

`manifest.json` records exact compressed sizes and hashes. The gzip blobs are
addressed by the SHA-256 of their uncompressed bytes. Files are reconstructed by
concatenating the ordered chunk list; equal entire files and equal chunks share
storage. Every blob is smaller than 50,000,000 bytes. Compression uses gzip level
6 and a zero timestamp. No record, number, newline or absolute path inside a
sealed file is rewritten.

There are 356 references to byte-identical files already published in the
original audit or the corrected occupancy shadow package. The index records
their repository-relative provider, encoding, and original and encoded hashes.
Native executables and disposable build/runtime caches are omitted explicitly;
binary identities remain in the original receipts and omission inventory.
External directory links are recorded literally but never followed or restored.
The published files contain no native executables. Exact Python source is kept
as `.py.txt` or inside blobs so publication does not alter the executed source
to satisfy later lint rules.

## Verify and restore exact artifacts

Run from a complete repository checkout with Python 3. Restore the shared audit
payload first using its [publication instructions](../../../scientific_accuracy_audit/REPRODUCE.md).
The corrected occupancy shadow dependency files must also be present in
`docs/performance/evidence/corrected_occupancy_shadow/`.

```sh
python3 docs/scientific_accuracy_audit/publication/manage.py restore
python3 docs/performance/evidence/accepted_c1/verify_publication_v2.py.txt verify \
  --package docs/performance/evidence/accepted_c1 --shared-root "$PWD"
```

Verification checks compressed and uncompressed chunk hashes, complete restored
file hashes, the native launch matrix and raw stream bindings, all observation
lane bindings, CLI raw artifact bindings and the acceptance/failure/resolution
receipts. The additive v2 verifier also binds the complete inventory and metadata
to the exact original plan, recomputes storage counts, and checks that the package
metadata stays unchanged throughout verification/restoration. It does not execute Rust, run scientific references or measure native
performance. Its retained result is [VERIFICATION.json](VERIFICATION.json).

Restoration must target a **new** directory. This example creates an unused
child of a temporary directory; it does not change existing `.stericx` captures:

```sh
c1_restore_parent=$(mktemp -d)
python3 docs/performance/evidence/accepted_c1/verify_publication_v2.py.txt restore \
  --package docs/performance/evidence/accepted_c1 --shared-root "$PWD" \
  --output "$c1_restore_parent/restored"
```

The origin workspace prefix is
`/run/media/andrej-rumenovski/New Volume/Code/StericX`.
For artifact access, map that prefix to the new restoration root and preserve
every following repository-relative component. For example, the historical
`.../StericX/.stericx/profiling/scientifically_validated_optimization/accepted_c1_v1/acceptance.json`
is restored beneath that same relative path. The verifier performs this mapping
through the index; it never opens an embedded absolute path. Shared providers
are resolved beneath `--shared-root`. Absolute paths, parent traversal and
symlink escapes are rejected for restored/provider paths.

Sealed JSON still contains the original absolute strings after restoration.
Do not edit those strings or point old controllers at existing captures.
`RESTORATION_RECEIPT.json` names the fresh destination and lists links that were
not created and binaries/caches that were omitted. Restored file permissions and
modification times are retained where supported; scientific identity is checked
against exact bytes, not filesystem timestamps.

The synthetic verifier/restorer tests are archived as
[test_publish.py.txt](test_publish.py.txt). Run their restored originals with:

```sh
python3 "$c1_restore_parent/restored/.stericx/profiling/scientifically_validated_optimization/publication_controller_v1/test_publish.py"
```

These tests cover path traversal, parent symlink escape, corrupted/missing chunks,
compressed and uncompressed hash mismatch, bounded decompression, shared-provider
identity, preservation of embedded absolute strings, refusal of an existing
destination, exact chunk concatenation and safe pause behavior.

The separate [v2 tests](test_verify_publication_v2.py.txt) add omission, changed
inventory metadata, forged storage counts, altered plan bytes/scope, metadata
replacement, modified controller and verified-index reuse cases. All 10 original
tests passed again from the fully restored 4.57 GB tree; all 9 additional tests
passed. The original v1 controller, compressed plan, index and blobs remain
unchanged. [RESTORATION_V1.json](RESTORATION_V1.json) retains the successful full
restoration receipt; v2 verification independently checks the same full package.

Run the additive tests from the repository root using fresh copies:

```sh
c1_test_root=$(mktemp -d)
cp docs/performance/evidence/accepted_c1/verify_publication_v2.py.txt "$c1_test_root/verify_publication_v2.py"
cp docs/performance/evidence/accepted_c1/test_verify_publication_v2.py.txt "$c1_test_root/test_verify_publication_v2.py"
python3 "$c1_test_root/test_verify_publication_v2.py"
```

## Recreate the publication

The exact builder is [publish.py.txt](publish.py.txt). The compressed plan binds
the original source inventory, publication specification, shared manifests and
builder hash. Its original commands were:

```sh
python3 .stericx/profiling/scientifically_validated_optimization/publication_controller_v1/publish.py plan \
  --repo "$PWD" \
  --spec .stericx/profiling/scientifically_validated_optimization/publication_controller_v1/spec.json \
  --output .stericx/profiling/scientifically_validated_optimization/publication_controller_v1/plan_v1.json
python3 .stericx/profiling/scientifically_validated_optimization/publication_controller_v1/publish.py build \
  --plan .stericx/profiling/scientifically_validated_optimization/publication_controller_v1/plan_v1.json \
  --output docs/performance/evidence/accepted_c1 \
  --journal .stericx/profiling/scientifically_validated_optimization/publication_controller_v1/build_journal.jsonl \
  --pause-file .stericx/profiling/scientifically_validated_optimization/publication_controller_v1/PAUSE
```

For a new assembly use new output/plan/journal paths. Relocate the specification's
controller file paths if needed, then generate a **new plan** against the restored
root; do not rewrite the retained plan. The shared provider repository files must
be available beneath that root. New plan timestamps and path metadata will differ.
With the same uncompressed bytes and gzip implementation, content-addressed blob
contents and whole restored file identities are reproducible. Creating `PAUSE`
stops before the next file/chunk with exit 75; completed evidence is retained and
`build --resume` can continue the same incomplete plan after the pause ends.

[STORAGE_REPORT.json](STORAGE_REPORT.json) lists the ten largest compressed
contributors, attributing each unique blob to its first indexed owner. For
future C1b/C2/C3 publications, pin this package's manifest and restoration-index
hashes as a shared provider. Match complete files by SHA-256 and size, then reuse
their existing chunk sequence; remaining equal chunks can also reference this
provider. Do not republish equal observation or prediction streams. A later
publication must explicitly support these confined provider references and
verify both packages; this original v1 builder is preserved unchanged.

## Fresh scientific and native reproduction

Artifact restoration is separate from executing the experiment again. Follow
the [corrected scientific reproduction recipe](../../../scientific_remediation/REPRODUCE.md)
with fresh build, observation and reference destinations. The restored
`candidate_c1_science_v2/complete.json` contains the original argument arrays,
environment overrides and engineering command logs. Its source/build manifest
identifies the accepted C1 source; the original baseline and workload manifests
identify the required baseline source, input files and toolchain. The archived
controllers provide the exact historical orchestration, including independent
reference and model/thermodynamic review programs.

The candidate source tree is restored under
`candidate_c1_science_v2/validated_build_corrected_vC1/source` beneath the common
optimization prefix. The original descriptor implementation and both dependency
files are preserved under `candidate_c1_ordered_files/before/`; the retained
`candidate.patch` and pre-edit/implementation manifests bind the change. Use these
in a new scratch tree and verify its complete source inventory against the frozen
baseline/candidate manifests before building. Do not build in the restored archive.

Use a clean scratch source tree, restore the shared inputs by their recorded
hashes, and build new binaries with pinned Cargo and Python dependencies. Preserve
full observer/CLI/prediction vectors and independent reference results before
admitting a new build for timing. Run paired native measurements only in a quiet
window, retaining warmups, both execution orders, every sample and every failure.
The original benchmark manifests preserve all actual workload arguments and
pair orders. [benchmark_corrected_pair.md](../../../../scripts/benchmark_corrected_pair.md)
documents its strict scientific admission requirements.

A fresh path, toolchain, build timestamp or model metadata timestamp can change
receipt or binary hashes. Input relocation can also change emitted filenames.
Record an explicit mapping and new receipts; compare unchanged scientific
fields and IEEE bits under the established contracts. Such a replay cannot
claim the historical raw-byte identity unless its raw hashes actually match.
