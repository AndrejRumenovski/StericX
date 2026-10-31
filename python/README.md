# StericX Python development package

A thin Python interface to the Rust steric-descriptor engine. This development
version is not a published PyPI release. Build a local wheel with
`maturin build --release --manifest-path python/Cargo.toml` from the repository,
or `python -m pip install ./python` for a source installation with Rust installed.

Python 3.11+ is the initial compatibility target. Actual tested platforms and
versions will be listed in the distribution validation report; a CI matrix is
not evidence that every platform has passed.

Source coordinates use Å and zero-based indices. Native molecular descriptors
retain f32 kernels; statistics use f64. Density mode is a separate validated
node-quadrature calculation with explicit units, threshold and provenance.
Missing conformers and unavailable populations remain visible. See the repository
scientific methods and independent-validation reports for the admitted scope.
