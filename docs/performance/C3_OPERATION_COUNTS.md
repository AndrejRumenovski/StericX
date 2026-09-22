# Actual C3 occupancy operation counts

The frozen C3 candidate executes the row-rejection work predicted by the earlier
corrected-baseline shadow. This diagnostic establishes source-operation counts
and exact outputs on the two fixed workloads. It establishes no native speedup,
memory saving, or candidate acceptance; the complete scientific gates and paired
native measurements remain separate.

The probe runs the actual C3 helper and loop copied with counters into a separate
source tree. The admitted buried-volume source is retained as an exact byte prefix
and still produces every CLI output. Each instrumented frame is compared with
the unmodified C3 function and the retained original shadow: all 30,168 frames
and 452,520 f32 fields match bit for bit. Full native/probe stdout and stderr are
identical and match the original same-thread scientific fingerprints.

The 10,000-file throughput workload repeats 56 conformers under distinct copied
paths. It is not 10,000 independent chemistries. Full input paths, conformer and
orientation order, unique frame identities and explicit neighbor order are
checked against the old records; aggregate matching alone cannot pass.

| Source operation | 56 conformers | 10,000 files |
| --- | ---: | ---: |
| Frame orientations | 168 | 30,000 |
| Integration point visits | 2,588,544 | 462,240,000 |
| Rows | 124,320 | 22,200,000 |
| Original/current XY preparations | 3,837,215 | 685,507,811 |
| Original cached Z predicates | 10,257,800 | 1,831,868,275 |
| Actual C3 cached Z predicates | 5,162,061 | 921,865,331 |
| Original newly prepared Z predicates | 626,895 | 111,964,217 |
| Actual C3 newly prepared Z predicates | 291,152 | 51,997,436 |
| Actual endpoint predicates | 484,253 | 86,487,808 |
| Actual row/atom rejections | 335,743 | 59,966,781 |
| Actual Z predicates including endpoints | 5,937,466 | 1,060,350,575 |

The actual 10,000-file loop removes 883,481,917 Z predicates after charging
86,487,808 endpoint predicates: **45.4505% fewer Z predicates**, or **33.6009% fewer
combined XY/Z preparations and predicates** under the explicitly counted model.
Every frame's actual cached and new Z counts equal the old counts minus the
old shadow's would-eliminate counts. XY preparation/rejection, hit/miss and
integration-point counts also match per frame.

## Row-discovery and finite-check work

These are actual C3 counts. The full-row scan adds a traversal, but not every
listed check is net new: row discovery replaces the old per-point row-change
test. A tuple-XY check is counted once, and a vector-finiteness call is counted
once; neither is presented as a hardware instruction count.

| Actual C3 source operation | 56 conformers | 10,000 files |
| --- | ---: | ---: |
| `row_scan_points` | 2,588,544 | 462,240,000 |
| `row_bound_checks` | 2,588,544 | 462,240,000 |
| `row_xy_comparisons` | 2,588,376 | 462,210,000 |
| `row_first_point_finite_checks` | 124,320 | 22,200,000 |
| `row_remaining_z_finite_checks` | 2,464,224 | 440,040,000 |
| `row_helper_calls` | 626,895 | 111,964,217 |
| `helper_row_finite_checks` | 626,895 | 111,964,217 |
| `helper_atom_finite_checks` | 626,895 | 111,964,217 |
| `helper_radius_finite_checks` | 626,895 | 111,964,217 |
| `helper_xy_finite_checks` | 626,895 | 111,964,217 |
| `inside_z_range_lazy` | 142,642 | 25,476,409 |
| `nonfinite_fallback_lazy` | 0 | 0 |

Helper calls include inside-range and nonfinite exits; endpoint distance
predicates count only the branch that actually evaluates the endpoint. The
ordinary workload has no nonfinite fallback; focused tests exercise that path.
Added scans, branches, finite checks and cache effects prevent converting the
predicate reduction directly into a wall-time prediction.

## Provenance and limits

The exact C3 buried-volume prefix SHA-256 is
`1df92d49859c3d1a8a17feeaa9d3ac49409ca275feea3b8d5fdd53f6635bb7d9`. The probe manifest is
`candidate_c3_occupancy_counts_v2/manifest.json`, SHA-256
`485188e4d2aff482453effb3e36e9c2c43856d2ba6960a436daa6485c3dc8887`, under the optimization evidence root.
Its inventory retains source, instrumentation diff, test/build commands, binaries'
identities, inputs, raw frame counters and complete output bytes. The compact
[count data](current_c3_operation_counts.json) retains every aggregate without
rounding and binds those receipts.

The probe's two focused Rust tests, one retained boundary/nonfinite test and nine
Python comparison-contract tests passed. Attempt v1 compiled the probe CLI but
could not compile library tests because the build snapshot omitted two test-only
C2 fixtures; no workloads ran in that attempt. V1 remains unchanged. V2 copies
exact archived fixtures from C3's pinned pre-edit manifest, with no expected
value changes, and retains that explicit resolution.

The probe deliberately performs extra calculations and writes a sidecar. Its
elapsed time and allocations are excluded from performance or memory claims.
Actual C3 allocation effects require the separate diagnostic binary and matching
capture phases; executable/session-path byte differences are not algorithmic
memory savings. The separate frozen native/diagnostic/symbol builds completed,
and all eight one-million-value prediction exports match the original baseline.
