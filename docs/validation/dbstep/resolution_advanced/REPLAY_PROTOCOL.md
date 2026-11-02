# Advanced functionality replay protocol

This replay uses the existing scientific controllers, input geometries, descriptor
settings and acceptance criteria without modification. Original strict hard-sphere
DBSTEP failures remain failures in this lane. Completion of a command is distinct
from scientific admission, and a known failing lane does not prevent independent
lanes from running. The revised point-domain and continuum policy is assessed
separately after approval; it does not rewrite these results.

The new orchestration script is `scripts/replay_dbstep_resolution_advanced.py`.
Each phase refuses an existing output directory, captures full stdout/stderr,
exit codes, commands, source and executable hashes, resource limits, and before /
after source inventories. Runtime observations collected while other validation
runs are active are not performance benchmarks.

## Commands

Run from the repository root, using the already pinned DBSTEP / Morfeus environment
and fresh destinations:

```sh
.venv/bin/python scripts/replay_dbstep_resolution_advanced.py references \
  --build .stericx/dbstep_resolution/baseline/build \
  --output .stericx/dbstep_resolution/advanced_v1/references
.venv/bin/python scripts/replay_dbstep_resolution_advanced.py distribution \
  --output .stericx/dbstep_resolution/advanced_v1/distribution
.venv/bin/python scripts/replay_dbstep_resolution_advanced.py installs \
  --distribution .stericx/dbstep_resolution/advanced_v1/distribution \
  --output .stericx/dbstep_resolution/advanced_v1/installs
```

The controller records exact expanded commands. All subprocesses use CPUs 0,1,
Rayon / Cargo at two threads, and BLAS / OpenMP / NumExpr at one thread.

The references phase freshly compiles four Rust examples with the locked offline
dependency graph. It executes DBSTEP Sterimol, hard-sphere volume, density volume,
radius profiles, independent Decimal populations, and the Morfeus conformer study.
The study recomputes independent references against preserved conformer reports;
it does not itself regenerate Rust ensemble outputs. Current compiled Python
bindings are separately tested against all 56 frozen conformer rows by installed
wheel tests. Historical Kraken comparisons retain their provenance limitations.

## Installation and CI

The distribution phase rebuilds a manylinux 2.28 / CPython 3.11 ABI3 wheel using
maturin and Zig, and a source archive. The installation phase creates five new
virtual environments: wheel installation on CPython 3.11, 3.12, 3.13 and 3.14, plus
source installation on 3.12. Each checks dependency-free base import, optional
NumPy / pandas / RDKit adapters, the copied installed-package test suite, and
`pip check`. Artifacts, package versions and tests are hashed.

The existing `.github/workflows/python-wheels.yml` will run on the dedicated
`dbstep-resolution-validation-20261102` branch after the parent coordinator freezes
and supplies its committed source SHA. Its twelve wheel jobs cover Linux, Windows
and macOS with CPython 3.11–3.14; a thirteenth job independently builds, installs and
tests the source archive. The remote branch is for validation only. No publication,
release tag or update to remote main is part of this replay. A successful local
Linux test cannot substitute for those actual remote platform runs.

All original evidence is preserved. Results, including failures and incomplete CI
jobs, will be reported explicitly; no revised tolerances are introduced here.

## CI evidence collection and review

The approved source snapshot was `5c521fe6bf864d7abeebdaacf9a1bcc02b5c5814`.
The dedicated-branch push receipt records the exact command. The existing push
triggers started wheel run `41484774517` and Rust run `41484774606`.
Use `gh run view RUN_ID --json databaseId,workflowName,headSha,headBranch,status,conclusion,url,createdAt,updatedAt,jobs`
and `gh run view RUN_ID --log` to retrieve the saved metadata and logs. Download
wheel/source artifacts with `gh run download 41484774517 --dir DESTINATION`;
the artifact API response records their remote identities and retention dates.

`scripts/review_dbstep_resolution_ci.py --input RAW_CI_DIR --output NEW_REVIEW_DIR`
checks the source SHA, 3 Rust / 13 wheel-source job results, 71 installed tests per
wheel/source job, artifact counts, wheel ABI tags and RECORD digests, committed
Python source identity and source-archive Rust / binding / test identity. Input
`source_bindings.json` consists of SHA-256 over exact `git show COMMIT:path` bytes.
Raw source hashes always remain recorded. Python text differing solely by exact
CRLF-to-LF serialization can be source-equivalent only when the normalized bytes
match the Git source digest and ASTs match with source-location attributes omitted.
The first byte-exact review's eight Windows line-ending failures remain in
`ci_v1`; the explicit serialization review is preserved separately as `ci_v2`.
No numerical gate or artifact was modified in this review.
