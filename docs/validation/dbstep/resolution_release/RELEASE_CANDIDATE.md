# StericX Python 0.4.0rc1

**Prepared and unpublished.** The approved scientific policy, controlled
performance review, actual platform tests and artifact checks pass. The
[final admission](admission_v2/review.json) rechecks 2,703 evidence records.
The Rust scientific engine remains 0.3.0; the binding crate is 0.4.0-rc.1.

The package source is commit `03557ed5279abf8487012ddb39fab019f0538416`.
Its only changes from the validated measurement source are four version
declarations. Scientific code, constants, numerical defaults and package test
assertions are unchanged. Later commits add evidence, validation workflows and
release documentation.

## Installable artifacts

The [candidate archive](candidate_v1/candidate_artifacts.tar.gz) contains exactly
these four installable files, plus a manifest, checksums and instructions:

| Artifact | Validated target |
|---|---|
| `stericx-0.4.0rc1-cp311-abi3-manylinux_2_28_x86_64.whl` | Linux x86-64, Python 3.11–3.14; manylinux 2.28 compatibility tag |
| `stericx-0.4.0rc1-cp311-abi3-win_amd64.whl` | Windows x86-64, Python 3.11–3.14 |
| `stericx-0.4.0rc1-cp311-abi3-macosx_11_0_arm64.whl` | macOS arm64, Python 3.11–3.14; executed on macOS 14 |
| `stericx-0.4.0rc1.tar.gz` | Clean source installation tested on Python 3.12 with Rust |

Extract the archive into a new directory, verify [SHA256SUMS](candidate_v1/SHA256SUMS),
then install the wheel matching your platform:

```sh
python -m pip install ./stericx-0.4.0rc1-cp311-abi3-manylinux_2_28_x86_64.whl
python -c "import stericx; print(stericx.__version__)"
```

That command shows the Linux filename; select the Windows or macOS filename
from the table on those systems. The unpacked local bundle is also available at
`.stericx/releases/stericx-0.4.0rc1/`. Its
[manifest](candidate_v1/manifest.json) binds each selected artifact to the actual
CI output and release review. The archive is 2,989,553 bytes, SHA-256
`262791009a166331a7aa1a62982b52626035ca09679d11463af5ce440cc37868`;
all seven archive members were verified after writing.

## Gates and retained failures

* Scientific: approved common-domain, boundary, analytical/convergence,
  invariance, conformer and independent-reference gates pass. All 128,021 legacy
  observations remain unchanged. The [scientific report](../FINAL_RESOLUTION.md)
  retains original DBSTEP failures and the resolution-dependent accuracy limits.
* Performance: all 504 native attempts, 28 allocation runs, 96 Python blocks and
  two ensemble allocation runs are preserved. No unchanged 10% trigger fires.
  Measurements used the validated dev0 package; the RC metadata transition and
  rebuilt artifacts are separately checked, not described as separately timed.
* Local packaging: five clean installations each pass 71 tests, base imports
  without optional dependencies, dependency checks and exact version assertions.
  [Local evidence](local_v1/RESULTS.md) also checks wheel RECORD and 64 source members.
* Remote packaging: all 13 wheel/source jobs pass 71 tests each. Linux/macOS
  runtime payloads match across four interpreter builds. Windows build variants
  differ in import/debug metadata, so the exact selected wheel is additionally
  installed and tested on all four supported interpreters, with 71 tests each.
* Rust platforms: Linux/Windows pass; the initial Mac job never obtained a
  runner and remains failed infrastructure evidence. The identical suite and
  source actually pass on macOS 14 arm64: 347 tests and a 22-test screening repeat.
  [CI audit](RELEASE_CI.md) records the separate real runs and cancelled duplicates.
* Final research/style checks: 156 discovered tests, 143 passes and 13 explicit
  reference-environment skips; the mandatory pinned-reference suite separately
  passes 28 tests with no skips. Ruff/format pass. The final checks bind 600 source
  files before and after execution without changes.

The [first coordinator attempt](admission_v1/failure.json) is retained. It checked
the overall workflow conclusion as `cancelled`, confusing it with the Mac job
conclusion. The frozen metadata agree on overall `failure` and Mac `cancelled`.
The corrected reviewer requires those exact values plus the complete successful
Mac execution. No test count, numerical tolerance or scientific result changed.

## Reproduce review and bundling

Use new output destinations; original evidence is immutable. The exact final
review entry point is:

```sh
.venv/bin/python scripts/review_dbstep_resolution_release.py \
  --science .stericx/dbstep_resolution/admission_v2/review.json \
  --performance docs/validation/dbstep/resolution_performance/final_v1/review.json \
  --local .stericx/dbstep_resolution/release_candidate_v1/local/receipt.json \
  --ci docs/validation/dbstep/resolution_release/ci_v1/ci_review.json \
  --ci-versions docs/validation/dbstep/resolution_release/ci_v1/rc_versions.json \
  --windows-wheel docs/validation/dbstep/resolution_release/same_artifact_v1/same_artifact_review.json \
  --macos-completion docs/validation/dbstep/resolution_release/macos_fallback_v1/macos_fallback_review.json \
  --final-checks .stericx/dbstep_resolution/release_candidate_v1/final_checks_v3/receipt.json \
  --output docs/validation/dbstep/resolution_release/admission_repeat
.venv/bin/python scripts/package_stericx_release_candidate.py \
  --review docs/validation/dbstep/resolution_release/admission_v2/review.json \
  --versions docs/validation/dbstep/resolution_release/ci_v1/rc_versions.json \
  --output .stericx/releases/stericx-0.4.0rc1-repeat \
  --archive docs/validation/dbstep/resolution_release/candidate_repeat
```

No PyPI upload, GitHub Release or public release tag was created. Publication
requires separate user approval, a fresh package-name availability check and
post-publication installation verification. No registry installation is claimed.
