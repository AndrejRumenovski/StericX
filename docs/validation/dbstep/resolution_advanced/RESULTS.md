# Advanced descriptor and distribution replay

The fresh local replay preserves all prior scientific results. The original strict
DBSTEP hard-sphere lane still fails on the same radius-2.5 Å case. This replay does
not apply the proposed replacement policy and does not convert that failure into
a pass.

| Check | Fresh outcome |
| --- | --- |
| DBSTEP Sterimol | 61 / 61 cases completed; conditional numerical gate passes; documented axial-convention differences preserved |
| DBSTEP hard-sphere volume | 69 cases / 207 frames; original 3 frame and 2 summary failures preserved exactly |
| Electron-density volume | 12 / 12 declared-grid cases pass; threshold's physical choice remains outside validation scope |
| Radius profiles | 20 positions: 16 reference comparisons and 4 expected errors preserved; strict hard-sphere radius-2.5 position remains failing |
| Boltzmann populations | 42 valid and 11 invalid cases pass independent Decimal checks; largest absolute weight residual 6.05 × 10⁻¹⁷ |
| Conformer study | 11 ligands / 56 conformers; independent references and ensemble reductions have no failures or investigations |
| Linux wheel installation | Python 3.11.15, 3.12.13, 3.13.14 and 3.14.4: all 71 tests pass in each new environment |
| Linux source installation | Fresh CPython 3.12.13 source compilation and installation: all 71 tests pass |

The five comparable scientific summary files—Sterimol, hard-sphere volume, density,
profiles and the conformer study—are byte-for-byte identical to the preserved
historical summaries. Installed package tests independently preserve every one of
the 56 frozen native conformer descriptor results, including binary32 bits. They
also check failed conformer retention, ordering, thermodynamic provenance, free
energy temperature mismatches, structured errors, tables, RDKit and density APIs.

Every local installation imports the dependency-free base package before optional
NumPy / pandas / RDKit dependencies are added, uses isolated `python -I` test
execution against a copied suite, and passes `pip check`. The manylinux 2.28 ABI3
wheel was rebuilt using valid unchanged Cargo incremental artifacts; source-archive
installation compiled from a new temporary source tree. No scientific or binding
source changed during these phases. The orchestration receipt records intervening
analysis-only Git commits separately from the unchanged scientific source hashes.

The build environment lacked `pip`, so its requested `python -m pip list` inventory
command failed. That failure and its log are preserved. A successful standard
library `importlib.metadata` inventory supplements it without altering the build
environment. Wheel and source build commands, all installation commands and all
installed-package tests succeeded. The other nonzero command is the deliberately
preserved original strict hard-sphere comparison.

## Reproducible local evidence

[Protocol](REPLAY_PROTOCOL.md) gives the exact entry points. The controller records
expanded commands, resource limits, package versions, hashes, output and errors.
[Machine-readable summary](local_v1/summary.json) separates every gate and records
historical comparisons. [Raw evidence manifest](local_v1/raw_evidence_manifest.json)
contains every file hash; [the archive](local_v1/raw_evidence.tar.gz) preserves all
143 logs, receipts, geometries, requests, raw responses and comparison tables.
Each archive entry has been checked against the manifest.

Fresh local artifacts remain in `.stericx/dbstep_resolution/advanced_v1/distribution/`:

- Wheel: `stericx-0.4.0.dev0-cp311-abi3-manylinux_2_28_x86_64.whl`, SHA-256
  `690d160a28b0e3914ccbfb835ef47b17c10b6175c58470f57cb246d999adc6f7`.
- Source archive: `stericx-0.4.0.dev0.tar.gz`, SHA-256
  `d2d733a3a25d82a5b2abc59089440e6b98ae6f5cb8a497bb6ef67d2019976df6`.

These validation timings were collected alongside other workloads and provide no
performance acceptance evidence.

## Fresh cross-platform CI

Snapshot `5c521fe6bf864d7abeebdaacf9a1bcc02b5c5814` was pushed to the dedicated
`dbstep-resolution-validation-20261102` branch with explicit coordinator approval.
No remote main update, release tag or publication occurred.

- [Wheel and source matrix](https://github.com/AndrejRumenovski/StericX/actions/runs/41484774517):
  **all 13 jobs passed**: 12 wheel jobs across Linux / Windows / macOS and
  CPython 3.11–3.14, plus an independent source-installation job. Every installed
  wheel and the source installation passed all 71 tests.
- [Rust and quality checks](https://github.com/AndrejRumenovski/StericX/actions/runs/41484774606):
  **all 3 jobs passed**: Linux quality / tests and Windows / macOS tests.

[Reviewed CI evidence](ci_v2/ci_review.json) binds 200 source / package / test /
example / Cargo / workflow files to the exact tested commit. The 12 downloaded
wheels and source archive are hashed individually. All 96 wheel entries were
checked against RECORD, 36 Python source entries against the committed source,
and 58 source-archive scientific / binding / test files against the commit.
Full run metadata, installed-package logs and source inventories are preserved
in [ci_v2](ci_v2/manifest.json); local copies of downloadable artifacts remain at
`.stericx/dbstep_resolution/advanced_v1/ci_v1/artifacts/`.

The [initial artifact review](ci_v1/ci_review.json) is preserved with eight raw
source-byte mismatches: `__init__.py` and `_types.py` in each Windows wheel use
CRLF line endings. The raw hashes remain distinct. Exact CRLF-to-LF substitution
produces the committed bytes in all eight cases, and Python ASTs are identical
before / after that substitution. The final review explicitly records this text
serialization equivalence; all other source bytes must match exactly. No wheel,
scientific result, reference output or numerical tolerance was changed.
