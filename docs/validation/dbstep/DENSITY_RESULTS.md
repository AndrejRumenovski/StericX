# Validated density-grid buried volume

The restricted density mode passes all **12 independent comparisons**:
three official DBSTEP total-SCF-density cubes at R=1.0, 1.5 and 2.0 Å, and
three analytical Gaussian grids. Sphere-node and occupied-node counts agree
exactly with unmodified DBSTEP support functions and its public API.
Percentages and explicitly unit-corrected grid volumes pass the unchanged
`64*epsilon_f64*max(1,abs(value))` allowance. No criterion was fitted to results.

| Cube | Radius (Å) | Sphere nodes | Occupied nodes | Buried volume (Å³) | Buried (%) |
|---|---:|---:|---:|---:|---:|
| Ne_medium.cube | 1 | 6100 | 6100 | 4.18487112 | 100 |
| Ne_medium.cube | 1.5 | 20603 | 19726 | 13.5329127 | 95.7433383 |
| Ne_medium.cube | 2 | 48797 | 19726 | 13.5329127 | 40.4246163 |
| benzene_coarse.cube | 1 | 751 | 751 | 4.12171793 | 100 |
| benzene_coarse.cube | 1.5 | 2559 | 2219 | 12.1785514 | 86.71356 |
| benzene_coarse.cube | 2 | 6081 | 3792 | 20.811657 | 62.3581648 |
| CH2CMe3_100.cube | 1 | 5283 | 5283 | 4.18708611 | 100 |
| CH2CMe3_100.cube | 1.5 | 17832 | 15771 | 12.4994388 | 88.4421265 |
| CH2CMe3_100.cube | 2 | 42287 | 28599 | 22.6663781 | 67.6307139 |
| gaussian_0.25.cube | 0.95251898 | 1551 | 251 | 0.581161914 | 16.1831077 |
| gaussian_0.125.cube | 0.95251898 | 12533 | 2103 | 0.608657124 | 16.7797016 |
| gaussian_0.0625.cube | 0.95251898 | 100137 | 17071 | 0.617593067 | 17.0476447 |

Volumes use original grid cells and node quadrature. They are not the existing
hard-sphere analytic-sphere-normalized volumes. Native CODATA 2022 conversion
is 0.529177210544 Å/Bohr; DBSTEP uses 0.529177249. Both are recorded, and the
volume conversion factor is explicit. No constants were changed to obtain
agreement. The original quantum methods for the official fixtures are not
invented; their comments identify total SCF density and their exact files and
official MIT license are retained.

The analytical Gaussian superlevel set has radius one Bohr and known volume.
The boundary-layer bound is derived before measurement from half the voxel
diagonal, not from observed errors.

| Spacing (Bohr) | Absolute volume error (Å³) | A priori boundary-layer bound (Å³) |
|---|---:|---:|
| 0.25 | 0.0395527525 | 0.818930942 |
| 0.125 | 0.0120575424 | 0.40474087 |
| 0.0625 | 0.00312159917 | 0.20177986 |

The errors decrease on this fixed refinement sequence; arbitrary grids need
not converge monotonically. Six Rust tests additionally cover cube ordering,
original vectors/units, malformed and incompatible fields, strict threshold
equality, empty/full occupation, translation, determinism, unsupported grids,
missing sample nodes and incomplete domains. Clippy passes with warnings denied.

The first controller run completed calculations but failed while serializing
NumPy boolean check values. Its source, partial output and logs remain in
`controllers/validate_density_volume_v1.py` and `density_v1/`. The corrected
controller converts those checks to Python bool; the native implementation
and scientific criteria were unchanged.

Reproduce with a freshly built `cargo build --example density_report`, then:
```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .stericx/paton/reference/venv/bin/python scripts/validate_density_volume.py \
  --output docs/validation/dbstep/density_replay
```

See [method and admission contract](../../methods/electron_density_buried_volume.md),
[summary](density_v2/summary.json), [reference outputs](density_v2/reference_raw.json),
[receipt](density_v2/receipt.json), and [complete artifact hashes](density_v2/manifest.json).
All request settings, source/binary hashes, exact cube files and raw native
responses are preserved. Public support is limited to declared total density
in scalar Gaussian cubes, positive axis-aligned isotropic grids, fully covered
spheres and explicit thresholds. Other grid types fail clearly. Physical
threshold universality, quantum-method accuracy and catalyst prediction are
not validated by this numerical agreement.
