# Native point-decision observation adapter

This is the system under investigation, **not an independent reference**. The
existing remediation adapter is copied unchanged except for an include and two
dispatch arms in `src/main.rs`. Its historical geometry and raw-occupancy
operations remain available. The generated build contains a byte-identical copy
of production `src/geometry/buried_volume.rs` followed by the unchanged historical
audit append. `src/resolution_append.rs` only invokes private native functions;
it does not reimplement an occupancy predicate.

Input is newline-delimited JSON. `grid_occupancy` takes `config`, already-aligned
`atoms` with `position: [x,y,z]` and `radius_squared`, and `include_points`.
It invokes native `validate_config`, `integration_grid` and `occupied_volumes`.
`raw_point_decisions` additionally takes `points` and sends them unchanged to the
same native aggregation; it does not apply integration-sphere membership again.
Input numeric values are converted to f32 exactly as in the historical adapter.
Squared radii should be captured native intermediates, with original radii/scales
retained externally. This adapter does not infer chemistry or geometry frames.

Both operations return native aggregate volumes and binary32 bits, accepted
aligned atoms/bits, grid length, native point quadrant/octant IDs and counts,
occupancy mask, and optionally points/bits. Each mask decision comes from a
separate actual `occupied_volumes` call on the single native point; zero versus
the sphere volume identifies occupancy. Positive normal sphere volume and a
nonempty point set are required so underflow cannot fabricate a mask. Returned
regional totals preserve the native sign-zero conventions. Full aggregate and
pointwise masks are captured separately so aggregation or row-cache defects are
detectable; matching totals alone is insufficient.

Build reproducibly with a new target and evidence directory:

```sh
.venv/bin/python docs/validation/dbstep/resolution_observer/build_observer.py \
  --baseline .stericx/dbstep_resolution/baseline/build \
  --oracle .stericx/profiling/scientifically_validated_optimization/scientific_oracle_v3 \
  --output .stericx/dbstep_resolution/mask_build \
  --evidence .stericx/dbstep_resolution/mask_build_evidence
```

The controller verifies the complete production source set against the frozen
baseline, retains the exact adapter diff/hash and build commands, then checks
the compiled observation source's unchanged production prefix. Build logs,
snapshot source, dependency locks and binary hashes are preserved outside this
adapter directory. Do not modify bound sources or adapter files during a build.
