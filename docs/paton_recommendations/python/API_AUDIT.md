# Python binding audit and core API

The baseline had no PyO3 bindings or importable `stericx` package. Its root
`pyproject.toml` remains the existing research-tools environment. The new
`python/` mixed maturin project depends directly on the root Rust crate through
a separate Cargo lockfile. It does not implement a second descriptor engine.

| Capability | Implementation and verification |
|---|---|
| Molecule loading | Native XYZ/V2000 reader; `Molecule` requires one record; `Ensemble` preserves all records and failures |
| Descriptors | Existing native Sterimol, volume/regions and pyramidalization; typed units/results |
| Ensembles | Rust statistics, explicit missing-data policy, optional checked populations and sensitivity CSV |
| NumPy/pandas | Optional imports; ordered conformer rows, numeric missing values plus explicit error columns |
| RDKit | Existing coordinates serialized with 17 significant digits into checked native parsing; authoritative adjacency and original conformer IDs; input object unchanged |
| Batch | Ordered Rayon collection, explicit optional worker count, failed files retained |
| Errors | Input, calculation, population and density exception classes; per-descriptor typed failure records |
| Precision | Native f32 values widened exactly before JSON transport; no decimal-rounding drift on Python decode |
| Runtime | GIL released during native descriptor, batch and density calculations |

RDKit integration intentionally avoids V2000 coordinate rounding, embedding,
optimization, conformer pruning or energy inference. The core stores adjacency,
not bond orders, charge or stereochemical identity; a matching adjacency check
is not proof of identical chemical species. XYZ inputs cannot establish mapped
regional ensemble correspondence. Density/grid and radius-profile public
wrappers are a separate advanced API admission step.

Core verification: **64 tests pass**, including all 56 frozen native MMFF
conformers with exact widened f32 descriptor values, analytical single-sphere
Sterimol, failed-middle conformers, missing populations, analytical degeneracy,
transactional metadata assignment, structured failures, tables, RDKit precision
and identity, and ordered batch failures. These transport checks complement the
independent scientific validations; they do not replace them.

The development package version is `0.4.0.dev0` (Cargo `0.4.0-dev.0`), following
the existing engine's `0.3.0` version. This is neither a release candidate nor a
PyPI publication. PyPI's name endpoint returned 404 at the recorded check time;
that does not reserve the name. The ABI3 target starts at CPython 3.11; actual
OS/version install tests, source distribution and release gates are recorded
separately. Optional NumPy, pandas and RDKit dependencies are not required to
import the base descriptor package.
