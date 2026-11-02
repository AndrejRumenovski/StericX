# Python distribution validation

A release-mode `stericx-0.4.0.dev0-cp311-abi3-manylinux_2_28_x86_64.whl`
built with PyO3 0.28.2 and maturin installs successfully in clean Linux x86_64
environments. The separate source archive builds and installs without access to
the repository's original source paths. Nothing was uploaded to PyPI.

| Artifact / interpreter | Base install and import | Complete adapter suite | Dependency check |
|---|---|---|---|
| Wheel / CPython 3.11.15 | Pass, without optional packages | 71 passed | Pass |
| Wheel / CPython 3.12.13 | Pass, without optional packages | 71 passed | Pass |
| Wheel / CPython 3.13.14 | Pass, without optional packages | 71 passed | Pass |
| Wheel / CPython 3.14.4 | Pass, without optional packages | 71 passed | Pass |
| Source archive / CPython 3.12.13 | Build/install/import pass | 71 passed | Pass |
| CI Linux x86_64 / Python 3.11–3.14 | Pass, four jobs | 71 passed per job | Pass |
| CI Windows x86_64 / Python 3.11–3.14 | Pass, four jobs | 71 passed per job | Pass |
| CI macOS arm64 / Python 3.11–3.14 | Pass, four jobs | 71 passed per job | Pass |
| CI source archive / Python 3.12 | Build/install pass | 71 passed | Package suite passes |

The 355 passing installation tests include NumPy, pandas and RDKit adapters,
exact native descriptor values on all 56 frozen conformers, structured failures,
thermodynamic provenance and validated density/profile cases. Tests execute from
a copied suite using isolated Python mode; the module path is recorded to guard
against importing the repository source instead of the installed package.
Each environment starts empty, installs the base artifact using **pip**, verifies
that optional packages are absent, then installs adapters and runs the complete
suite. Exact Python and dependency versions, commands and stdout/stderr are in
the [installation receipts](distribution_receipt.json).

The first Python 3.11 adapter download failed with a network `IncompleteRead`.
That environment and log remain in `install-3.11/`; a separate fresh environment
`install-3.11-v2/` passed. No scientific test failure was discarded. The first
source archive included unrelated local README/license files because Cargo
include patterns were unanchored. Its inventory is retained; the corrected
archive has 131 members, no local environments/build directories/user media,
and passed clean-source installation.

## Build and reproduce

From the repository, with Rust and maturin installed:

```bash
maturin build --release --locked --manifest-path python/Cargo.toml --out dist
maturin sdist --manifest-path python/Cargo.toml --out dist
python -m pip install dist/stericx-0.4.0.dev0-*.whl
```

The admitted Linux portability build additionally used `--zig --compatibility
manylinux_2_28`. Maturin checked the platform tag; no audit bypass was used.
Zig 0.17.0 came from the isolated `ziglang` package. The baseline native compiler
and mathematical constants were not changed. Source installations build for the
local host and need Rust; wheels do not need Rust or optional chemistry packages.

`python/pyproject.toml` declares CPython 3.11–3.14 (`>=3.11,<3.15`), ABI3 from
3.11, optional table/RDKit/notebook dependencies and the development version.
The root research `pyproject.toml` is unchanged. All original registry packages
in the core Cargo lock retain their versions/checksums in the binding lock.
The source package includes the shared engine, not a Python chemistry rewrite.

```bash
python scripts/validate_python_install.py --python 3.12 \
  --artifact /absolute/path/to/stericx-wheel.whl \
  --environment /new/isolated/environment --output /new/evidence/directory
```

This local controller targets Linux environments. The separate
[wheel workflow](../../../.github/workflows/python-wheels.yml) configures actual
Windows, macOS and Linux runner jobs for all four Python versions, plus a source
archive build/install job. It uploads CI artifacts and does not publish to PyPI.
The actual [validated run](https://github.com/AndrejRumenovski/StericX/actions/runs/41384461108)
passes all 13 jobs on commit `aa78ecb09d337488255937f6eb85e9027f0984db`.
All 12 wheels carry the expected `cp311-abi3` tag. Raw logs, runner/interpreter
versions, installed package versions, wheel-member hashes and downloaded artifact
identities are retained in [the CI receipt](ci/artifact_receipt.json). The separate
[Rust CI](https://github.com/AndrejRumenovski/StericX/actions/runs/41384461152) also
passes on Linux, Windows and macOS. These are actual execution results, not
inferences from a configured matrix.

The initial CI run passed Windows/macOS and source installation, but Linux selected
the container's Python 3.9 through a bare `--interpreter python` and built an
incompatible cp39 wheel. The corrected workflow selects `pythonX.Y` explicitly
inside manylinux and asserts the ABI3 wheel tag. This follows the
[official maturin-action interpreter guidance](https://github.com/PyO3/maturin-action/blob/main/README.md).
The initial failed run and complete logs remain in `ci/`; no chemistry or
numerical acceptance rule was changed.

PyPI's `stericx` JSON endpoint returned 404 at the recorded check time. This is
neither a reservation nor permission to claim `pip install stericx` currently
works from PyPI. Local wheel/source installation is the tested adoption path.

Artifacts remain under `.stericx/paton/python-wheels/`; their exact SHA-256 hashes,
wheel-member hashes, source snapshot hashes and install evidence are in
[distribution_receipt.json](distribution_receipt.json). The original host wheel
(`manylinux_2_35`) and development wheels remain separate from the tested portable
release-mode wheel. The tested development distribution is ready for publication review. No PyPI
publication or release candidate is implied; [release readiness](../RELEASE_READINESS.md)
separately records the retained scientific gates and limitations.
