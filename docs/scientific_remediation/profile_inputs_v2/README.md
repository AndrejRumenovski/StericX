# Corrected profiling inputs, preparation v2

Preparation completed successfully against frozen corrected build v4, without
running benchmarks. The [verification receipt](receipt.json) binds the completed
manifest, exact input bundle and helper, eight passing tests, and successful lint
and formatting checks. All 41,718 manifest files were reverified. Preparation v1
and historical workloads remain unchanged.

The new tree is
`.stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2`.
Its manifest is compatible with the unchanged profiling harness. It contains:

- Two single-geometry jobs, the original 56 conformers, a 10,000-file descriptor
  job and an 86-conformer original bonded SDF ensemble (Kraken molecule 1170).
- A 1,000-row parse job, fresh packed records repeated to one million predictions,
  a 1,000-row ensemble screen and a 56-geometry/11-ligand database build.
- Search against a newly calculated 1,543-ligand database containing every one
  of the 31,618 conformers admitted by the predeclared whole-ensemble topology
  rule. All 23 excluded molecules and reasons are enumerated before native
  execution. There were zero skipped selected geometries. Search query and
  database use the same coordination-axis convention.

The model is freshly fitted on the original ten training records and one blind
record, with explicit `supplied_weight_mean` and 353.15 K response metadata.
The historical 298.15 K supplied populations and original decimal XYZ bytes are
preserved. The script verifies exact sealed source scientific tokens, geometry
bytes/order, source-temperature evidence and the unresolved ligand-2064 note;
equivalent metadata bundles may reside at different paths, with their actual
manifest hashes retained.

Full source/build/evidence proof is checked before and after preparation.
Six executing helpers, the native binary and all 67 copied Ni-hDA geometries
have 74 explicit source/copy pairs checked before native use and at completion.
All four native preparation commands—parse, fit, model validation and full
database build—succeeded. Any failed selected geometry would abort preparation;
the selection is never narrowed after execution.

Cyclic repetition adds throughput volume, not chemical diversity. Database
means are unweighted. Population transferability, empirical calibration and
unresolved source limitations are not established by preparation. The full
unfiltered topology oracle is separate. Scientific admission and optimization
acceptance remain separate gates; no comparison to the old invalid workload
inputs or speedup claim is made here.
