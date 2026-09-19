# Exact-output performance experiments — September 16 baseline

**Accepted: 2.37× faster native descriptor batches and 3.41× faster screening**, with exact scientific outputs, preserved validation/errors and essentially unchanged memory. These ratios compare the final build with the original September 16 executable rerun on September 19. Against the original September 16 timing series, the improvements are 2.36× and 3.25× respectively. Both changes pass their end-to-end and exactness gates; no failed production experiment remains.

The screening path avoids unused calculations after an exact acceptance witness while retaining required validation. The buried-volume kernel reuses exact X/Y prefixes and rejects atoms for whole rows, reducing point–atom work without changing any scientific parameter. Detailed admission estimates, measurements and regression follow-ups follow below.

The September 16 fixed input manifest and frozen native outputs remain the source of truth. Current pre-optimization source (66 files), current executable, original September 16 executable and baseline summary hashes were frozen in `.stericx/profiling/optimization/` before computational edits. The original manifest SHA-256 remains `5f2a35784d64c12122a886db1d6bdce1e2d7b012cece31256d22b175f7d1b20c`.

## Candidate 1: screening-specific validation and Sterimol

Admission decision, recorded before implementation: **experiment justified**. Screen worker buried-volume work is 88.1186% of diagnostic process wall (elimination ceiling 8.4165×); occupancy alone is 84.6541% (6.5164×). These worker intervals overlap the main-thread library wait and are not added to it.

Full removal would change validation: `compute_from_center` rejects positive first-orientation volume with zero maximum adjacent-quadrant difference across all orientations. A simple removal is rejected. Instead compute the first orientation using the exact original kernel and adjacent-difference arithmetic. Once positive asymmetry is established, later occupancy scans cannot reverse the final zero/nonzero predicate. If first volume is not positive, that predicate cannot reject either. Continue every coordinate-basis and atom validation in original orientation order. If no acceptance witness occurs, retain the remaining exact scans and the original rejection. Remove pure unused pyramidalization and full buried-volume result assembly from screening only; keep parser, donor/axis choices, Sterimol arithmetic and conformer mean order.

With comparable orientation costs and an acceptance witness in the first orientation, skipping two of three occupancy scans removes approximately 56.44% of process wall: an estimated 2.30× ceiling for this candidate's principal change. Realistic admission estimate: 45–56% lower screen wall (1.8–2.3× speedup), allowing validation/setup and fallback scans. This is a prediction, not a benchmark result. No gain is assumed from changing density, points, radii or arithmetic. Every remaining validation/error and all outputs must match the frozen oracle; fallback cases require focused differential tests.


## Occupancy investigation (before production implementation)

The frozen baseline occupancy span is 92.51% of descriptor-batch diagnostic wall, with a 13.35× end-to-end ceiling assuming infinitely fast occupancy. In the original loop each point repeats atom X/Y differences and squares, despite consecutive points sharing identical X/Y coordinates. A prospective exact algorithm can prepare the original `(dx*dx + dy*dy)` prefix per atom per XY row, retain atoms in input order, reject only prefixes strictly greater than radius squared, and apply the unchanged `prefix + dz*dz <= radius_squared` test to survivors. No square roots, approximate bounds, point reordering, reassociation or radius changes are needed. Non-finite cases must preserve the original predicate.

Before enabling this in production, an isolated diagnostic copy of the frozen source records every original distance predicate/early exit and shadow-checks the row method's occupied boolean and first-hit atom. The shadow is measurement/feasibility analysis; native output continues to use the original kernel. This investigation preceded the measured admission decision for Candidate 2 below. Radius scaling and coordinate transformations already occur outside the point loop, so optimizing those alone is not a high-impact candidate. Grid caching's 1.035× ceiling is below the requested 10% target and is excluded.

### Screening experiment result

**Accepted:** screening falls from 1,253.096 ms to 549.458 ms, a 2.2806× native end-to-end speedup (56.15% less wall time), matching the admission estimate. The feature trace records exactly 1,000 occupancy calls instead of 3,000, while all later structural validation still runs. Native peak RSS is effectively unchanged (9,334,784 → 9,314,304 bytes).

Every full fixed-workload output fingerprint matches, including packed records. The additional frozen screening oracle exercises 70 cases at both one and six threads with exact stdout, stderr and exit status, without normalization. Both loaders/axes, failures in later conformers, donor/geometry/configuration failures, exclusions, existing explicit descriptor values, uncertainty/applicability and ranking are covered. All 259 Rust tests, Clippy and strict rustdoc pass at this checkpoint.

| Workload | Fresh control median (ms) | Candidate 1 median (ms) | Wall change |
|---|---:|---:|---:|
| `conformers_56` | 68.402 | 68.040 | -0.53% |
| `db_build_ensembles` | 69.709 | 68.880 | -1.19% |
| `descriptors_10000` | 12,061.264 | 11,936.483 | -1.03% |
| `ensemble_sdf` | 16.865 | 17.481 | +3.66% |
| `parse_ensembles_1000` | 180.551 | 183.339 | +1.54% |
| `predict_1000000` | 21.710 | 21.640 | -0.32% |
| `screen_1000` | 1,253.096 | 549.458 | -56.15% |
| `search_database` | 7.305 | 7.482 | +2.42% |
| `single_large` | 3.236 | 3.115 | -3.74% |
| `single_small` | 2.374 | 2.524 | +6.34% |

The initial seven-run six-thread database series appeared 13.06% slower. That isolated result was not discarded: a follow-up of 21 alternating baseline/candidate pairs gave 15.311 ms versus 14.928 ms, median paired ratio 0.975, with roughly 0.14 ms MAD for each binary and exact outputs throughout. The apparent regression did not reproduce. The largest positive absolute difference is parsing at 2.788 ms (+1.54%); the largest relative difference is the small molecule at 0.151 ms (+6.34%). All workloads are rechecked in the final combined candidate below.

Pure unused pyramidalization is omitted. The first complete buried-volume orientation remains required for general screening validation; later occupancy scans are skipped only after the exact rejection predicate is disproven. This is not an unconditional geometry-validation bypass.

## Candidate 2: exact XY-row filtering

Admission decision, recorded before production implementation: **experiment justified**. The isolated frozen-source probe completed both the 56-conformer and 10,000-molecule workloads with identical scientific output fingerprints and zero mismatches in occupied booleans or first-hit atom indices. Across all 30,000 orientations of the larger workload, it observed 462,240,000 integration-point visits and 11,159,539,632 original point–atom distance tests. Eager row preparation requires 692,455,740 XY-prefix evaluations and 1,943,834,272 surviving Z tests: 2,636,290,012 combined predicate evaluations, 76.38% fewer. Surviving per-point Z tests are 82.58% fewer than the original full distance tests. Every orientation improves this combined count.

The default grid still tests 32³ lattice candidates and retains 15,408 integration points per orientation, grouped into 740 consecutive XY rows. On the 56-conformer sample, unoccupied points constitute 67.2% of points and consume 86.2% of the original distance tests. The row bound rejects 83.5% of atom/row pairs, leaving a mean 5.15 candidates per row. Coordinate alignment and radius preparation already occur once per atom/orientation; their scale is tiny compared with repeated point–atom tests.

With occupancy accounting for 92.51% of original diagnostic wall, complete elimination has a 13.35× end-to-end ceiling. Treating each counted predicate as equally expensive gives a deliberately crude 3.41× projection from the 76.38% reduction. Cached prefixes also remove repeated X/Y subtraction, multiplication and addition, while candidate iteration, row detection and scratch-buffer traffic add costs. A realistic admission range is **2–4× native descriptor-batch speedup (50–75% lower wall)**, to be accepted only by end-to-end measurements. Counts and this estimate are not timing results.

Production will prepare each row lazily: cache retained atoms in their original order only as far as a point search needs to visit them. Subsequent points first search that cache, then continue through unseen atoms if necessary. Thus dense geometries that repeatedly hit the first atom avoid eagerly preparing every atom in the row. Prefix work is at most the measured eager preparation count, and surviving Z tests are unchanged. One reusable scratch vector holds candidate Z positions, squared radii and XY prefixes; no per-point allocations are needed.

The exact locked `glam` predicate is `(dx*dx + dy*dy) + dz*dz <= radius_squared`. Reject only when the original rounded XY prefix is strictly greater than squared radius; adding nonnegative squared Z cannot restore membership. Retain NaN-comparison and infinity-versus-infinity cases conservatively; an infinite prefix greater than a finite squared radius is safely rejected. Use the original arithmetic grouping for survivors. Group only consecutive points with identical X/Y bits. Atom order, point order, boundary comparison, integration density, radii, precision and final quadrant/octant arithmetic remain fixed. Independent old-kernel differential tests will cover finite boundaries and exceptional arithmetic.

The original September 16 binary/outputs and pre-optimization source remain immutable. Candidate 1 is separately frozen at `.stericx/profiling/optimization/candidate1/` with source hashes, native and diagnostic binaries, its complete seven-run native results and 70-case screening oracle comparisons. Candidate 2 must match both the original scientific outputs and these error/validation oracles.

### Production-algorithm operation audit

A second isolated probe uses the exact captured production lazy kernel, with counters and a separate first-hit index map. It checks every point against the original kernel and the earlier eager shadow. Both 56-conformer and 10,000-molecule scientific fingerprints match the frozen original; all occupied booleans and first-hit indices match. This probe supplies operation counts, not performance timings.

| Count on `descriptors_10000` | Original | Lazy row filter |
|---|---:|---:|
| Integration-point visits | 462,240,000 | 462,240,000 |
| Original full point–atom distance predicates | 11,159,539,632 | — |
| XY-prefix predicates | — | 685,507,811 |
| Surviving per-point Z predicates | — | 1,943,834,272 |
| Combined predicate evaluations | 11,159,539,632 | 2,629,342,083 |

The lazy variant eliminates **76.44% of predicate evaluations**, including preparation. It prepares 6,947,929 fewer prefixes than the eager upper-bound experiment. Successful point searches stop at exactly the same original atom. The original short-circuit already avoided 3,258,511,776 of 14,418,051,408 possible point–atom comparisons; the new algorithm removes work beyond that existing optimization. Of the 462,240,000 point visits, 149,295,103 searches (32.30%) terminate strictly before the last atom, 2,261,696 first hit the last atom, and 310,683,201 are unoccupied. Thus 151,556,799 points are occupied. These counts are summed from each orientation’s first-hit histogram; the weighted histogram independently reproduces all original distance checks. Unoccupied points account for 86.18% of original tests, explaining why rejecting atoms for an entire row is effective.

Each scan uses 23–48 eligible atoms in this corpus. Original tests average about 24.14 atoms per point. Coordinate/radius preparation processes only 935,751 eligible atom instances across 30,000 orientations: 2,807,253 relative-coordinate subtractions, 8,421,759 alignment multiplications, 5,614,506 alignment additions, and 935,751 each of radius scaling and squaring. Those operations were already outside the point loop and remain unchanged. Element-radius lookup occurs during atom construction, not during point membership tests.

The scratch vector contains three f32 values per candidate and is reused across rows. The probe observes at most 48 reserved entries (**576 bytes**) and 13 populated entries (156 bytes) per orientation, with no capacity growth in either workload. There is one candidate-vector allocation per scan, no per-point or per-row allocation. These are production-layout scratch sizes; diagnostic probe RSS includes unrelated counters and is not a memory-performance result. Native RSS and normal profiling allocator totals are reported separately below.

No integration-grid, density, radius, point-order, precision, boundary or descriptor-definition change was made. SIMD, square-root bounds and general spatial trees were not needed: this exact row partition removes most comparisons without approximate bounds or reordering.

## Final native end-to-end results

Seven measured launches after one warmup, with the unchanged September 16 manifest and original executable as control. The September 19 control accounts for conditions after the session resumed. The single-thread runs use CPU 2; six-thread runs use physical CPUs 0–5 on the same Ryzen 5 5600G, Rust 1.97.0 and Linux 7.0.0-31. Native builds have profiling disabled. Process startup and regular-file output are included; these are warm-cache runs, with no forced fsync or cold-storage claim. Builds, probes and profilers never overlap accepted native runs.

All ten single-thread workloads are retained here, including initially noisy results. Values are milliseconds; `±` is median absolute deviation, not a confidence interval. Full ranges, CPU time, page faults, RSS and per-launch records are in [optimization_summary.json](optimization_summary.json) and the evidence archive.

| Workload | Original September 16 median | September 19 control median ± MAD | Final median ± MAD | Final/control wall change |
|---|---:|---:|---:|---:|
| `conformers_56` | 70.799 | 70.056 ± 1.622 | 29.812 ± 0.024 | -57.45% |
| `db_build_ensembles` | 69.864 | 71.853 ± 0.547 | 30.902 ± 0.145 | -56.99% |
| `descriptors_10000` | 12,161.773 | 12,238.054 ± 42.343 | 5,164.105 ± 49.520 | -57.80% |
| `ensemble_sdf` | 17.098 | 17.143 ± 0.312 | 11.601 ± 0.524 | -32.33% |
| `parse_ensembles_1000` | 186.529 | 195.421 ± 1.798 | 262.844 ± 15.243 | +34.50% |
| `predict_1000000` | 21.876 | 23.339 ± 0.404 | 23.128 ± 0.653 | -0.90% |
| `screen_1000` | 1,262.674 | 1,323.825 ± 59.418 | 388.019 ± 9.309 | -70.69% |
| `search_database` | 7.478 | 8.358 ± 0.129 | 8.197 ± 0.704 | -1.92% |
| `single_large` | 3.236 | 3.334 ± 0.176 | 2.662 ± 0.201 | -20.17% |
| `single_small` | 2.487 | 3.017 ± 0.155 | 2.779 ± 0.219 | -7.89% |

The descriptor batch falls from **12.238 to 5.164 seconds (−57.80%)**; screening falls from **1.324 seconds to 388.019 ms (−70.69%)**. The 56-conformer and single-thread database workloads also improve by more than 56%. Small process launches benefit less because fixed startup/I/O costs remain. The occupancy candidate meets its pre-implementation 2–4× descriptor-batch estimate.

Six-thread sequential measurements:

| Workload | Original September 16 median | September 19 control median ± MAD | Final median ± MAD | Final/control wall change |
|---|---:|---:|---:|---:|
| `db_build_ensembles` | 15.245 | 15.378 ± 0.210 | 12.624 ± 1.753 | -17.91% |
| `predict_1000000` | 9.906 | 9.887 ± 0.356 | 12.497 ± 0.290 | +26.40% |

### Regression investigation

The sequential series showed +34.50% for parsing and +26.40% for six-thread prediction. Those samples remain in the table and raw evidence. The candidate parse series ranged from 190.18 to 285.78 ms. Follow-ups used **21 alternating original/candidate pairs**, reversing AB/BA order each pair, with every output checked and no samples removed. The slowdowns did not reproduce. Database build was included in the six-thread follow-up because its initial candidate times were also variable.

| Paired workload | Control median ± MAD (ms) | Final median ± MAD (ms) | Median candidate/control pair ratio ± MAD |
|---|---:|---:|---:|
| `parse_ensembles_1000` (1 threads) | 183.369 ± 1.636 | 186.188 ± 3.655 | 1.0125 ± 0.0207 |
| `db_build_ensembles` (6 threads) | 15.446 ± 0.260 | 8.758 ± 0.611 | 0.5628 ± 0.0585 |
| `predict_1000000` (6 threads) | 10.376 ± 0.491 | 10.435 ± 0.736 | 0.9918 ± 0.0695 |

The paired results show no meaningful repeatable regression: parsing's median pair difference is +1.25%, below its 2.07% ratio MAD; prediction's is −0.82%, with 6.95% ratio MAD. Six-thread database build falls from 15.446 to 8.758 ms. These follow-ups diagnose variability rather than replacing or suppressing the initial series.

### Memory and diagnostic attribution

Native descriptor-batch peak RSS is **19,021,824 → 19,001,344 bytes**; screening is **9,293,824 → 9,342,976 bytes**. Across the fixed sequential workloads, the largest median RSS increase is under 0.3 MB. There is no large index, duplicated integration grid or per-point allocation.

The ordinary profiling feature records 30,000 additional descriptor-batch allocations, one scratch vector per orientation. Total requested allocation traffic rises by about 11.23 MB (0.274%); peak live requested heap is 13,030,618 → 13,030,621 bytes. Screening allocation calls fall from 384,584 to 378,584 and its peak requested heap remains about 3.86 MB. Requested bytes describe allocator requests, not resident memory; the tiny residual byte differences include differing profiling metadata paths.

Separate diagnostic runs show `occupied_volumes` falling from 11.275 to 4.201 seconds in the 10,000-descriptor workload. Screening retains 1,000 occupancy calls instead of 3,000, with their summed worker time falling from 1.069 seconds to 137.273 ms. Worker spans overlap the main-thread library wait and are not added to it. These diagnostic times explain the change; the speedup claims use the native timings above.

### Instruction, memory-reference and branch audit

The final Cachegrind run uses the same complete 56-conformer and 10,000-descriptor workloads, optimized symbol builds with profiling disabled, the unchanged Valgrind executable and the original explicit cache model: 32 KiB/8-way/64-byte I1 and D1, 16 MiB/16-way/64-byte last-level cache. Both scientific fingerprints match. Raw event rows and independent `cg_annotate` function sums agree. Inlining boundaries changed, so the comparison below uses whole-process totals of equal scope.

These are **instrumented instruction/data-reference/branch counts; cache and branch misses are simulated**. Native hardware counters remain unavailable because the host denies `perf` access. `Dr`/`Dw` are not hardware load/store counters; Cachegrind counts a read-modify-write reference as a read and counts repeated string-operation iterations separately. See the [Cachegrind counting rules](https://valgrind.org/docs/manual/cg-manual.html#cg-manual.simulation).

| Event, complete `descriptors_10000` | Original | Final | Change |
|---|---:|---:|---:|
| Instruction events (`Ir`) | 238,829,039,500 | 68,100,602,668 | -71.49% |
| Data-read references (`Dr`) | 50,332,670,456 | 16,594,471,309 | -67.03% |
| Data-write references (`Dw`) | 1,264,305,814 | 3,357,986,193 | +165.60% |
| Conditional branches (`Bc`) | 25,164,808,354 | 8,849,149,349 | -64.84% |
| Indirect branches (`Bi`) | 22,402,745 | 22,452,359 | +0.22% |
| Simulated D1 misses | 124,858,590 | 125,422,650 | +0.45% |
| Simulated branch misses | 486,714,944 | 485,898,804 | -0.17% |

Scratch-buffer writes increase, while instructions, data reads and conditional branches fall substantially. Simulated cache and branch misses remain approximately flat. The evidence supports a reduction in work, not a claim of improved cache hit rates or branch prediction. Row detection, candidate-buffer writes, loop control and remaining mispredictions help explain why fewer predicates do not translate one-for-one into native speedup.

At the source-expression level, the counted distance subtractions and multiplications each fall from 33,478,618,896 to 3,314,849,894; distance additions fall from 22,319,079,264 to 2,629,342,083. These explicitly derived counts exclude row/index control and candidate writes, and are not machine-instruction counts. The archived first-hit and checks-per-point histograms, source counters, allocator diagnostics and Cachegrind totals keep those measurement scopes separate.

## Exactness and engineering gates

- Every native and diagnostic fixed-workload scientific fingerprint matches the original at its corresponding thread setting. Descriptor/search/screen JSON is compared in full, without numeric tolerance. Only explicitly named timing/resource/path metadata is removed from command summaries; requested thread settings still match. Complete v1 packed records and database CSV/manifest artifacts remain identical.
- The screening oracle freezes 70 cases and passes all of them at one and six threads with exact stdout, stderr and exit status. It covers both loaders and axes, donor/configuration/geometry failures, failures in later conformers, exclusions, explicit CSV descriptor precedence, uncertainty, applicability and ranking.
- The additional volume oracle freezes **47 complete v2 packed records (all 32 f32 fields per record) and 227 full-precision conformer rows (all nine public buried-volume fields)** across five cases. It covers all 56 real conformers with four configuration sets and three explicit centers. Both thread settings match every packed-record and conformer-CSV byte. Command stdout is compared after removing the explicitly named timing/RSS/destination-path metadata; stderr and exit status are also checked. The unrounded rows close the coverage gap left by rounded database CSV and descriptor JSON that omits octants and hemispheres.
- All four final million-value prediction exports (native/diagnostic × one/six threads) match the original 4,000,000 bytes, SHA-256 `35bc75c2dc0bffa8ba6599b0649479a30db0586f7edc4107f9f7141faa99a987`. Diagnostic exports actually enable profiling and allocation tracking, with zero dropped scopes.
- **264 Rust tests, five optimized-build kernel differential tests, 49 Python tests, Clippy with warnings denied, strict rustdoc and formatting pass.** The independent old-kernel tests compare every returned volume by f32 bits, including all 56 conformers across three orientations, adjacent float boundaries, signed zeros, subnormals, NaNs/infinities/overflow, randomized geometries, row resets and dense early hits. Study 011's frozen predictions, design lock and 57 planned splits pass verification. Ruff covers project code; pre-existing untracked `docs/media` drafts were excluded and left untouched.

The corpus is fixed and deliberately repeats chemistry: `descriptors_10000` cycles 56 real conformers, screening cycles 11 reaction rows under 1,000 identifiers, and the prediction batch cycles 11 records. It does not represent 10,000 distinct chemistries or the unavailable full Kraken quantum-conformer cache. These are measured workload-specific gains; very coarse grids may amortize row preparation less effectively. Required scientific parameters and all validation predicates remain unchanged.

## Reproduction and retained evidence

[optimization_summary.json](optimization_summary.json) retains the original and fresh control distributions, both candidates, complete scientific fingerprints, paired follow-ups, error/volume oracles, full-vector audit, operation counts and verification results. [optimization_evidence.tar.gz](optimization_evidence.tar.gz) contains the frozen source snapshots, diagnostic probe sources and patches, compressed per-orientation counters, per-launch measurements and traces, representative full workload outputs, all oracle artifacts, scripts and verification logs. Every member is hash-checked when packaging; [optimization_evidence_identity.json](optimization_evidence_identity.json) lists file hashes and the archive hash. Executables/build caches and repeated copies of identical large outputs are omitted; binary hashes and exact build commands are retained.

The original baseline is commit `37a3f0ae27b38ea4d67628fc84fd1ffc031462aa`, native binary SHA-256 `d69d2af29d8803611eef3403642b6941fdababa744d805807ba014aa42dfb01a`. The accepted final native binary is SHA-256 `b577c49d3e98b4745213fce0e70e93a55fc61120c5abe6fdb0658547c43cab94`; its source identity covers Cargo files and all Rust source/tests/examples used for the build. The original baseline and the screening-only checkpoint were frozen before subsequent computational changes. Admission decisions above were recorded before each production candidate.

The [baseline report](PROFILING.md#measurement-controls-and-reproducibility) describes input preparation, tool setup and measurement rules. Rebuild the original commit in a separate checkout and preserve its executable before building the final source. In this workspace, the frozen executable and manifest remain available; use new immutable run labels for rechecks:

```bash
cargo build --release --locked
python3 scripts/profile_stericx.py run --binary .stericx/profiling/optimization/september16_stericx --label recheck_control_t1 --reps 7 --threads 1 --affinity 2 --compare-to baseline_final_t1
python3 scripts/profile_stericx.py run --binary target/release/stericx --label recheck_final_t1 --reps 7 --threads 1 --affinity 2 --compare-to recheck_control_t1
python3 scripts/check_screening_exactness.py compare --candidate target/release/stericx --label recheck_final_t1 --threads 1
python3 scripts/check_volume_exactness.py compare --candidate target/release/stericx --label recheck_final_t1 --threads 1
```

Repeat the native database/prediction workloads and oracles with six threads and CPUs 0–5. Build diagnostics separately with `--features profiling` and `CARGO_PROFILE_RELEASE_DEBUG=1`; diagnostic elapsed time never substitutes for native latency. The full prediction audit uses `examples/profile_predictions.rs` and compares all exported bytes. Its exact commands and input hashes are archived.

For a different checkout path, run `profile_stericx.py prepare` into a fresh profiling root and create new oracle roots using each oracle's `prepare --baseline` action. Compare per-input content and workload definitions, then compare original/final outputs within that new manifest. Absolute invocation paths make the old manifest hash specific to this workspace; do not rewrite frozen evidence to accommodate a new path. The archived probe-source tarballs/patches reconstruct the instrumentation used for the counts, and their runner commands use the unchanged workload definitions.
