# Local release-candidate installation validation

**Passed:** Python package `0.4.0rc1`, native binding crate `0.4.0-rc.1`, unchanged
Rust scientific engine `0.3.0`. No publication was performed.

The [receipt](receipt.json) binds 90 tracked scientific, package, test, installation
controller and approved-policy files. Before and after bytes agree exactly. Against
validated commit `497de873`, the only permitted changes are four exact version
substitutions in `python/Cargo.toml`, `python/Cargo.lock`, `python/pyproject.toml`
and `python/stericx/__init__.py`. All other covered files match that commit byte
for byte. Existing scientific assertions and the 71-test package suite are
unchanged.

| Fresh installation | Package tests | Optional-free base import | pip check | Installed version assertions |
|---|---:|---|---|---|
| Linux wheel, Python 3.11 | 71 passed | Passed | Passed | Passed |
| Linux wheel, Python 3.12 | 71 passed | Passed | Passed | Passed |
| Linux wheel, Python 3.13 | 71 passed | Passed | Passed | Passed |
| Linux wheel, Python 3.14 | 71 passed | Passed | Passed | Passed |
| Source distribution, Python 3.12 | 71 passed | Passed | Passed | Passed |

Every installation first imports and calculates a descriptor without NumPy,
pandas or RDKit installed. It then installs the optional testing adapters, runs
the unchanged suite, checks dependencies, and verifies package `__version__`,
installed distribution metadata and native binding `__version__`. There are no
skipped package tests. All 16 orchestrated commands exit successfully; subordinate
installation commands and complete outputs are retained.

The clean source installation builds its own wheel through the standard build
backend. The release wheel uses maturin with `--release --locked --zig
--compatibility manylinux_2_28`, the existing build environment, and an isolated
Cargo target. Build environment packages are inventoried using
`importlib.metadata`; no pip executable is required in that build environment.
The installation environments are separately created with seeded pip.

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `stericx-0.4.0rc1-cp311-abi3-manylinux_2_28_x86_64.whl` | 583117 | `830e1aeb825875c2e94e6978f8025fc1c922ea08f735427ae9771cdac8cb7e1e` |
| `stericx-0.4.0rc1.tar.gz` | 1506346 | `4e922e997c02b527b42945a6a955ecb8a75224d4dec3183d4bcea674a8cc7b91` |

The separate [archive review](artifact_review.json) checks all eight wheel RECORD
entries, exact installed Python source bytes and 64 source-distribution members,
including scientific Rust sources, binding code, package code and test inputs.
Maturin's two documented packaging transformations are checked structurally:
`package.readme = "README.md"` is inserted in the binding manifest, and the
source archive's root project metadata points to `python/Cargo.toml`. Neither
changes calculation code.

Reproduce from the repository root, using new destinations:

```sh
.venv/bin/python scripts/validate_dbstep_rc_local.py \
  --output .stericx/dbstep_resolution/release_candidate_repeat/local \
  --archive docs/validation/dbstep/resolution_release/local_repeat
```

The retained `controller.py` is a byte-identical provenance copy of that script;
execute the repository script so its repository-root detection is correct. The
separate `artifact_review.py` is a read-only verifier for the original artifact
paths. `archive_manifest.json` and `artifact_review_archive.json` bind the compact
copies to the raw campaign. Full wheel/source artifacts and isolated environments
remain under `.stericx/dbstep_resolution/release_candidate_v1/local/`.

These five Linux installations do not replace actual Windows/macOS installation
CI or the final scientific/release review. Those are separate release gates.
