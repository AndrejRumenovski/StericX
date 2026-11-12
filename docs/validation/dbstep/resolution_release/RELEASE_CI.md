# Release-candidate platform validation

Package **0.4.0rc1**, binding crate **0.4.0-rc.1**, unchanged Rust engine **0.3.0**.
The tested scientific/package source is
`03557ed5279abf8487012ddb39fab019f0538416`. Its four changes from the validated
development package are version strings. No scientific implementation, constant,
default or numerical test was changed.

All required platform tests have actual successful executions. The first Rust
workflow retains its infrastructure failure; separate source-bound evidence
records the subsequent macOS execution. No synthetic combined workflow result
replaces either real run. Final release admission additionally requires the
coordinator's scientific, performance and artifact review.

## Actual executions

| Scope | Source or workflow run | Outcome |
|---|---|---|
| Linux Rust/quality | [41505795122](https://github.com/AndrejRumenovski/StericX/actions/runs/41505795122) | 348 tests plus the repeated 22-test screening guard pass; formatting, Clippy and documentation checks pass |
| Windows Rust | Same original run | 348 tests plus the repeated 22-test screening guard pass |
| Original macOS Rust attempt | Same original run | Hosted runner never acquired; infrastructure cancellation preserved |
| macOS 14 arm64 Rust | [41508791423](https://github.com/AndrejRumenovski/StericX/actions/runs/41508791423) | Actual 347 tests plus the repeated 22-test screening guard pass, zero failures/ignored tests; tracked source unchanged |
| Wheels and source archive | [41505795041](https://github.com/AndrejRumenovski/StericX/actions/runs/41505795041) | All 12 Linux/Windows/macOS × CPython 3.11–3.14 wheel jobs and the source-installation job pass 71 tests each |
| One identical Windows wheel across interpreters | [41507681984](https://github.com/AndrejRumenovski/StericX/actions/runs/41507681984) | The same wheel hash installs and passes all 71 tests on each of CPython 3.11–3.14 |

The arm64 suite has one fewer test because
`model::regress::tests::avx2_matches_scalar_fallback` applies only to x86. The
review checks the complete multiset of actual macOS test names against Windows,
with precisely that architecture-specific difference. This is the same coverage
as the earlier successful arm64 validation.

## Preserved infrastructure failure and actual replacement execution

The original `macos-latest` job waited 15 minutes and one second. GitHub reports
“The job was not acquired by Runner of type hosted even after multiple attempts”
and macOS arm64 capacity constraints. Its API record has `runner_id=0`, an empty
runner name, no steps and a cancelled conclusion. There was no test execution
to reinterpret as a scientific result.

[The original review](ci_v1/ci_review.json) deliberately remains `passed=false`,
with one failure: the Rust workflow/job conclusion was not successful. Its
original run metadata, annotations and successful Linux/Windows logs are
preserved. The full-run log download also failed because the never-started job
has no log; that error is retained alongside separately downloaded successful
job logs.

The additive macOS workflow uses `macos-14`, the arm64 runner platform already
used by the successful wheel tests. It checks out the same frozen source and
runs exactly the original cross-platform commands:

```sh
cargo test --all-targets --all-features
cargo test --test screening_regression
```

Explicit Bash `pipefail` preserves command failures while retaining complete
logs. The workflow also asserts arm64, records Rust/Cargo versions and checks
`git diff --exit-code` after testing. Its workflow commit is
`a371345a2f489744b3112f8cfc46197e086c67ab`; its checked-out scientific source remains
`03557ed5279abf8487012ddb39fab019f0538416`.
[The separate fallback review](macos_fallback_v1/macos_fallback_review.json)
requires both original platform successes, the exact capacity failure, actual
replacement success, source identity, commands and complete named-test coverage.
It does not rewrite the original failed workflow.

## Artifact identity and installation scope

[The artifact review](ci_v1/summary.json) verifies 13 downloaded distributions,
96 wheel RECORD entries, 36 Python source entries and 58 source-archive
scientific/binding/test files. It binds 214 source files to the tested commit.
[Version checks](ci_v1/rc_versions.json) independently verify wheel METADATA,
source PKG-INFO, Python `__version__`, project/crate/lock versions, and 13 actual
successful installation log lines. All raw artifact hashes remain available.

Within each of Linux and macOS, native extension bytes and Python runtime
payloads are identical across the four interpreter builds. Their SBOM and RECORD
entries differ and are preserved. Windows extension bytes differ between builds:
`.text`, `.data`, `.pdata` and `.reloc` match, while headers and `.rdata` differ.
No whole-binary equality claim is made for those Windows builds.

To validate a concrete distributable Windows wheel, the additional workflow
downloads the Python-3.11-built artifact from run 41505795041 and requires its
unchanged SHA-256:

`900fa8616446459cfd544edb5633c167b20bae0d0d0e01d19c78c9b0a6e69864`

Every interpreter independently checks the hash, installed runtime and package
metadata versions, `pip check`, and all 71 unchanged installed-package tests.
[The same-artifact review](same_artifact_v1/same_artifact_review.json) preserves
the four installation receipts, package inventories and test logs. Its workflow
commit is `6c2eddfadc3329b6ddf2e40dc2eb981ef865e4f2`; tests use the same frozen
scientific source. All four original Windows builds are also retained.

The tested architectures are Linux x86-64, Windows AMD64 and macOS arm64. No
additional architecture or interpreter support is inferred from the ABI tag.

## Reproduction and retained attempts

The committed workflows specify all runner, source, artifact and command inputs.
Raw downloads remain under
`.stericx/dbstep_resolution/release_candidate_v1/{ci_raw,same_artifact_raw,macos_fallback_raw}`.
For example:

```sh
gh run view 41505795041 --repo AndrejRumenovski/StericX --log
gh run download 41505795041 --repo AndrejRumenovski/StericX --dir WHEEL_EVIDENCE
gh run download 41507681984 --repo AndrejRumenovski/StericX --dir WINDOWS_RECEIPTS
gh run download 41508791423 --repo AndrejRumenovski/StericX --dir MACOS_RECEIPTS
```

The unchanged original artifact reviewer is
`scripts/review_dbstep_resolution_ci.py`; its expected nonzero result preserves
the original runner failure. The supplemental version, Windows and macOS review
programs are archived verbatim with their receipts and input hashes. Their
manifests preserve each committed evidence file.

Pushing the two additive workflows to dedicated validation branches also
automatically queued duplicate full pipelines. Only those duplicate runs
(41507681761, 41507681961, 41508791307 and 41508791361) were cancelled with explicit
coordinator authorization; their run metadata, commands and reasons remain
archived. Required platform tests were executed and reviewed separately. No
remote main change, release tag, GitHub Release or PyPI publication occurred.
