# StericX Python development package

A thin Python interface to the Rust steric-descriptor engine. This development
version is not a published PyPI release. Build a local wheel with
`maturin build --release --manifest-path python/Cargo.toml` from the repository,
or `python -m pip install ./python` for a source installation with Rust installed.

Python 3.11–3.14 is tested on Linux x86-64, Windows x86-64 and macOS arm64.
Each of the 12 actual CI wheel-install jobs passes all 71 package tests, including
NumPy, pandas and RDKit adapters. Source-archive installation also passes.
The repository distribution report preserves exact artifacts, versions and logs.

Source coordinates use Å and zero-based indices. Native molecular descriptors
retain f32 kernels; statistics use f64. Density mode is a separate validated
node-quadrature calculation with explicit units, threshold and provenance.
Missing conformers and unavailable populations remain visible. See the repository
scientific methods and independent-validation reports for the admitted scope.
