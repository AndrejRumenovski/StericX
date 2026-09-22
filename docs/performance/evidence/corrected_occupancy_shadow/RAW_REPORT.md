# Corrected occupancy shadow: diagnostic results

No production candidate was implemented. This report counts source-level operations;
it contains no native performance measurement or accepted speedup.

Both exact corrected-workload argument vectors were replayed, including all 10,000
throughput files (56 distinct content hashes repeated) and all 56 distinct conformers.
The baseline is accurate_baseline_v1 at corrected commit 1ecfbd5, with its admitted
scientific scope and retained numerical/convention limits. It is not the old
uncorrected audit baseline. The private copy's entire buried_volume.rs remains an
exact prefix. Four isolated files have additive diagnostic hooks; all other frozen
source bytes are unchanged. Original source and original predicates still produce
all CLI scientific output. The shadow never skips a predicate.

## Measured loop counts

All frames contain 15,408 actual integration points in 740 contiguous XY rows.
Row lengths range 2–30. These are actual accepted sphere points, not the 32³ grid's
candidate cube sites. Three frame scans reuse one generated grid per conformer.
Counts below describe one original computation; duplicated shadow/recheck work is
excluded. Each would-skip original predicate was executed and asserted false.

| Counter | Diverse 56 | Actual 10,000 files |
| --- | ---: | ---: |
| Conformers/files actually replayed | 56 | 10,000 |
| Frame orientations | 168 | 30,000 |
| Original Z predicates | 10,884,695 | 1,943,832,492 |
| Shadow-removable Z predicates | 5,431,482 | 969,969,725 |
| Added endpoint predicates | 484,253 | 86,487,808 |
| Rows | 124,320 | 22,200,000 |
| Point visits | 2,588,544 | 462,240,000 |
| Lazy XY preparations | 3,837,215 | 685,507,811 |
| XY-only rejections | 3,210,320 | 573,543,594 |
| Cached Z predicates | 10,257,800 | 1,831,868,275 |
| New-atom Z predicates | 626,895 | 111,964,217 |
| First hits in cache | 688,227 | 122,908,551 |
| First hits on new atom | 160,417 | 28,648,248 |
| Unoccupied point misses | 1,739,900 | 310,683,201 |
| Row/atom endpoint rejections | 335,743 | 59,966,781 |
| Inside-range retained row/atoms | 142,642 | 25,476,409 |
| Nonfinite fallback row/atoms | 0 | 0 |
| Would-remove cached Z predicates | 5,095,739 | 910,002,944 |
| Would-remove new-atom Z predicates | 335,743 | 59,966,781 |
| Eager-model XY preparations | 3,876,120 | 692,455,740 |
| Eager-model endpoint checks | 493,551 | 88,148,204 |
| Net Z predicates removed | 4,947,229 | 883,481,917 |
| Remaining Z predicates including endpoints | 5,937,466 | 1,060,350,575 |
| Eager extra XY preparations | 38,905 | 6,947,929 |

Net Z-predicate reductions are 45.451241%
and 45.450517%, respectively, after
charging every added endpoint predicate. No frame has a net Z-predicate increase.
XY preparations are unchanged in this lazy shadow model. Counting XY comparisons
plus Z comparisons plus added endpoint comparisons gives a reduction of
33.604532%
and33.600897%; this is a source predicate count,
not CPU instructions, cycles or wall time. The eager row×atom counts are a separate
upper-work model and are not proposed cache behavior. Original atom order, lazy
prefix extension and first hits remain unchanged.

## Static source operation model: not instrumented counters

The following are source-derived invocations applied to observed successful input
shapes. No allocator interception, radius lookup counter or hardware instruction
counter was enabled. Compiler elimination/vectorization and allocator reallocations
are not inferred. These rows exclude extra diagnostic recomputation.

| Source-level model | Diverse 56 | Actual 10,000 |
| --- | ---: | ---: |
| Aligned atom entries / coordinate transforms | 5,238 | 935,751 |
| Dot products (three per transform) | 15,714 | 2,807,253 |
| Radius scales and radius squares (each) | 5,238 | 935,751 |
| Parser radius-table lookups (one per input atom) | 3,097 | 553,228 |
| Grid construction Vec::with_capacity sites | 56 | 10,000 |
| Aligned Vec collection sites | 168 | 30,000 |
| Occupancy cache Vec::with_capacity sites | 168 | 30,000 |

`aligned_atoms` transforms each retained atom once per frame using three dot products,
then scales/squares its stored radius. The inner point loop performs no coordinate
transform or element radius-table lookup. XYZ parsing calls Atom::new once per input
atom; the stored radius then feeds geometry. Each occupancy frame starts a candidate
Vec with capacity min(atom_count,64), clears/reuses it across rows, and can grow if
needed. Actual allocation/reallocation counts are unknown. Other parsing, descriptor
and output allocations are outside this static table. Each original grid construction
reserves the cube capacity and retains 15,408 points; this is not one grid allocation
per row or point.

## Exactness and provenance

Two focused shadow test groups passed, covering boundary-adjacent floats, signed zero,
subnormal absorption/underflow, finite overflow, inside-range nonmonotonic rows,
reversed/interrupted rows, atom order, empty inputs, NaN and infinities. Added tests
explicitly assert endpoint outcomes; nonfinite mixed rows do not accidentally stand
in for finite-path coverage. The finite nearest-endpoint monotonic proof and required
fallbacks are recorded in PRE_EDIT_HYPOTHESIS.md before source instrumentation.
The exact unchanged predicate remains (dx*dx+dy*dy)+dz*dz <= radius_squared.

Every one of 30,168 frames matched all 15 f32 values bit for bit against unmodified
occupied_volumes: total/near/far, four quadrants, eight octants. Thus 452,520 recorded
u32 values passed runtime differential assertions. Complete stdout and stderr match
fresh frozen-native execution byte for byte, and both runs match the admitted
corrected native workload fingerprints. This does not assert equality with the old
scientifically incorrect aggregate convention: it uses the corrected baseline.

- conformers_56 stdout SHA256: `8098e16239bd7740ae4061868d83b45cd0f3989fa950bd869fdef72a16ea606b`; ordered 15-field frame-bit SHA256: `3f3e8848c8a17ea5d6882fe5343502b77a1fbec84f166c015c8af4ed0d189a00`.
- descriptors_10000 stdout SHA256: `df2815f9b1411e7d77a62de7bd1074747ff5cface781f62857e278553da37037`; ordered 15-field frame-bit SHA256: `9a2a9dc3d4d5f7229ff563faa20e8e10c23be54cb634a10951a011fbf589ecea`.

source_manifest.json records exact original/instrumented source inventories and patch.
build_manifest.json records compile/test commands, raw logs, Rust/Cargo versions and
binary identity. completion_binding.json independently rechecks build→source linkage,
zero exit codes, raw file identities, frozen workload/source inventories, each run's
pre-run binary identity and post-run native/probe identities. Executed matching bytes
are retained at bin/probe and bin/baseline_native outside disposable target/. The
runner executed the target-path probe; the archive is verified byte-identical afterward.
Two archived helper layout attempts are explained in helper_layout_receipt.json;
the first import failure happened before source preparation or execution.

Raw per-frame counts, atom hit/removal histograms, row-length histograms and bounded
elimination witnesses remain in both workload probe/counts.jsonl streams. The summary
also retains net-reduction extremes and all measured counters. Reading the receipts
and hashes is necessary when reusing this evidence; equal-looking aggregates alone
are not an admission gate.

## Conditional next step

The count reduction supports testing a C3 row-endpoint candidate after C1/C2 admission
and a fresh profile. It does not establish a 10% native improvement. Added min/max and
finiteness scans, endpoint branches, candidate memory access and cache behavior must
be measured. A candidate can scan each exact-XY row once, then preserve the original
point loop and lazy atom order; grouping the outer loop can remove its current
per-point row-change check. Retain same f32 values/grouping, equality, nonfinite fallback,
inside-range atoms, grid and regional denominators. Any optimization still requires
all science/independent-reference/engineering gates before quiet paired timings.
No production implementation or benchmark is authorized by this diagnostic result.
