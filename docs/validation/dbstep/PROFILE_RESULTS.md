# Fixed-center radius profiles

The public Rust APIs expose hard-sphere and validated density profiles through
`Conformer::buried_volume_scan` and `DensityGrid::buried_volume_scan`. The scan
variable is explicitly the sampling-sphere radius in Å. Center coordinates,
source geometry, radius requests and actual engine radii are recorded. Molecular
hard-sphere calculations preserve native f32 arithmetic; density calculations
use the separate validated f64 node quadrature.

Twenty requested positions were tested in four profiles: one hard-sphere ligand
and three official density cubes. Sixteen valid positions, including duplicates,
match their individual calculations exactly and were checked with fresh DBSTEP
calls. Four expected failures remain explicit records (negative radius or
out-of-domain density sphere). Nonfinite/empty radius lists fail before producing
a nonserializable or ambiguous profile. Caller order is preserved.

The hard-sphere R=2.5 Å case retains its known sphere-boundary discrepancy with
DBSTEP. Profile preservation passes; **exact common-grid agreement at every
hard-sphere radius does not**. See [discrepancies](DISCREPANCIES.md). This API
extends the existing method without changing it to fit a reference. All density
profile positions pass the original independent arithmetic criteria.

Two Rust tests cover fixed centers, order, duplicate requests, malformed source
geometry, invalid radii and absent domain coverage. Clippy passes. Complete
[summary](profiles_v1/summary.json), [reference calls](profiles_v1/reference_raw.json),
raw native outputs, exact source geometry and [hashes](profiles_v1/manifest.json)
are retained.

```bash
cargo build --example volume_profile
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .stericx/paton/reference/venv/bin/python scripts/validate_volume_profiles.py \
  --output docs/validation/dbstep/profiles_replay
```

No center-displacement or shell-volume interpretation is implied. A profile is
not an ensemble average and does not introduce thermodynamic populations.
