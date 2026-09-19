# StericX end-to-end profiling — 16 September 2026

Baseline: `37a3f0ae27b38ea4d67628fc84fd1ffc031462aa`. This report records the pre-optimization snapshot; no computational optimization had been made at this checkpoint. The `profiling` feature adds optional timing and allocation probes; the default build compiles the probes out. Recommendations below describe experiments proposed from that snapshot, not speedup claims. Formatting, Clippy, 252 Rust tests with all features, and 39 Python tests passed at this checkpoint. Subsequent changes and their exact-output gates are recorded separately in the [optimization report](OPTIMIZATION.md).

## Native process measurements

AMD Ryzen 5 5600G, Linux 7.0.0-31, Rust 1.97.0, locked release build. One Rayon thread pinned to CPU 2; seven measured launches after one warmup. Wall time includes process startup and output. CPU is user plus system time divided by wall time; 100% means one core. RSS is the exact child’s `wait4` high-water mark. MiB = 1,048,576 bytes. MAD is median absolute deviation, not a confidence interval.

| Fixed workload | Wall median ± MAD (ms) | Range (ms) | CPU (%) | Peak RSS median (MiB) |
|---|---:|---:|---:|---:|
| `conformers_56` | 70.799 ± 0.636 | 68.308–71.791 | 99.6 | 5.01 |
| `db_build_ensembles` | 69.864 ± 0.201 | 69.621–71.482 | 99.3 | 5.39 |
| `descriptors_10000` | 12,161.773 ± 30.416 | 12,131.357–12,380.137 | 99.9 | 18.34 |
| `ensemble_sdf` | 17.098 ± 0.115 | 16.904–17.470 | 98.5 | 4.92 |
| `parse_ensembles_1000` | 186.529 ± 2.812 | 181.594–190.327 | 99.8 | 4.74 |
| `predict_1000000` | 21.876 ± 0.322 | 21.487–22.411 | 99.3 | 69.51 |
| `screen_1000` | 1,262.674 ± 7.838 | 1,254.836–1,302.944 | 99.2 | 8.94 |
| `search_database` | 7.478 ± 0.076 | 7.336–23.466 | 96.1 | 6.75 |
| `single_large` | 3.236 ± 0.077 | 3.029–3.692 | 93.1 | 4.79 |
| `single_small` | 2.487 ± 0.076 | 2.375–2.836 | 91.7 | 4.92 |

Six threads pinned to six physical cores, CPUs 0–5, use the same inputs:

| Workload | One thread (ms) | Six threads (ms) | Measured speedup | Six-thread CPU (%) |
|---|---:|---:|---:|---:|
| `db_build_ensembles` | 69.864 | 15.245 | 4.58× | 499.7 |
| `predict_1000000` | 21.876 | 9.906 | 2.21× | 245.4 |

This is measured thread scaling, not an algorithm change. The repeat of the default build after adding optional probes stayed between −3.26% and +0.98% of baseline medians across the ten workloads; these short runs do not establish statistical equivalence of speed.

## Fixed workload coverage

- `single_small` and `single_large`: real 43- and 82-atom XYZ geometries, default bond-axis descriptors with buried volume and pyramidalization.
- `ensemble_sdf`: 12 conformers in one SDF; original coordinate decimal tokens are preserved in the parser’s supported whitespace-delimited V2000 representation.
- `conformers_56`: 56 actual RDKit/MMFF conformers from 11 ligands; `db_build_ensembles` groups the same 56 files into 11 database rows.
- `descriptors_10000`: 10,000 distinct files cycling the same 56 conformers. This measures scale and repeated geometry processing, not 10,000 distinct chemistries.
- `parse_ensembles_1000`: 1,000 rows cycling 11 real reaction records, with original conformer weights and 5,085 weighted geometry evaluations; exports every packed record.
- `predict_1000000`: one million packed records cycling 11 real records, with the portable model’s weights; includes memory mapping, inference, MSE validation and CLI reporting.
- `screen_1000`: 1,000 reaction-library rows cycling 11 records with unique identifiers. Missing Sterimol columns trigger one geometry per row; the screen includes uncertainty and applicability checks. This is not conformer-ensemble screening.
- `search_database`: one geometry query against all 1,541 real descriptor database rows, including JSON output.

These geometry workloads use phosphorus donors, the default bond axis, a 3.5 Å sphere and 0.01 integration density. They exclude external conformer generation and do not establish costs for non-phosphorus donors or the coordination axis.

Input creation never invokes StericX. Every source and generated file is hashed in [workloads.json.gz](workloads.json.gz); the immutable manifest SHA-256 is `5f2a35784d64c12122a886db1d6bdce1e2d7b012cece31256d22b175f7d1b20c`. Preparation reuses tracked source files. The manifest records absolute invocation paths: a checkout at another location will have a different manifest hash. Regenerate inputs there and compare per-file content hashes and workload definitions; never rewrite the accepted manifest in place. The original large Kraken quantum-conformer cache was unavailable, so these measurements must not be generalized to that full chemistry distribution.

## Wall-time attribution and theoretical limits

Five independent diagnostic launches per workload add lexical, nested wall-clock probes and requested-heap counters. Main-thread exclusive scopes partition the CLI root exactly; worker elapsed totals overlap the waiting parent and are listed separately. They are not CPU sampling, and are never added to main-thread wall shares. Inference spans include Rayon setup, allocation, demand paging and waiting as applicable (mmap setup is timed separately), rather than only SIMD arithmetic.

For each launch, `f = scope time / external diagnostic process wall`. The report takes the median of those fractions. The maximum speedup from eliminating that work is `1/(1−f)`; making it twice as fast gives `1/(1−f+f/2)`. These are ideal diagnostic-workload bounds. They assume everything else remains fixed, exclude replacement overhead, and are not achieved native speedups. Nested function ceilings cannot be added. Individual median rows need not sum to exactly 100%.

| Workload and scope | Median time (ms) | Diagnostic wall share | Infinite local speedup: end-to-end ceiling | 2× local speedup: end-to-end estimate |
|---|---:|---:|---:|---:|
| `descriptors_10000`: `buried_volume::occupied_volumes` self | 11,275.444 | 92.51% | 13.345× | 1.861× |
| `descriptors_10000`: `buried_volume::integration_grid` self | 416.568 | 3.42% | 1.035× | 1.017× |
| `descriptors_10000`: `parse_xyz` self | 148.299 | 1.22% | 1.012× | 1.006× |
| `descriptors_10000`: `sterimol::params_from_projection` self | 146.340 | 1.20% | 1.012× | 1.006× |
| `parse_ensembles_1000`: `sterimol::params_from_projection` self | 73.416 | 38.16% | 1.617× | 1.236× |
| `parse_ensembles_1000`: `parse_xyz` self | 54.392 | 29.20% | 1.412× | 1.171× |
| `parse_ensembles_1000`: `Molecule::from_xyz_file` self | 26.678 | 14.01% | 1.163× | 1.075× |
| `predict_1000000`: `model::RegressXPredictor::predict_batch` self | 15.451 | 67.66% | 3.092× | 1.511× |
| `predict_1000000`: `commands::predict::mean_squared_error` self | 2.310 | 10.08% | 1.112× | 1.053× |
| `screen_1000`: `commands::screen::load_library` self | 1,153.989 | 91.31% | 11.508× | 1.840× |
| `screen_1000`: `model::student_t_two_sided_quantile` self | 51.376 | 4.06% | 1.042× | 1.021× |
| `screen_1000`: `commands::screen::bootstrap_mean_response_interval` self | 32.178 | 2.50% | 1.026× | 1.013× |
| `db_build_ensembles`: `commands::db::featurize_wall` self | 69.226 | 93.35% | 15.030× | 1.875× |
| `search_database`: `commands::search::print_report` self | 2.156 | 19.91% | 1.249× | 1.111× |

`screen::load_library` self wall includes waiting for worker geometry calculation, which is not a nested scope on the main thread. It must not be read as 91% miscellaneous overhead. Its worker scopes spend about 1,113 ms in buried volume, including 1,069 ms in `occupied_volumes`, and about 14.6 ms in Sterimol. Likewise, database `featurize_wall` is the complete parallel phase; its worker scopes are a separate breakdown.

### Separated stage costs

| Stage | 10,000 descriptors (ms) | Parse 1,000 ensembles (ms) | Predict 1,000,000 (ms) | Screen 1,000 main thread (ms) | Screen workers, overlapping (ms) |
|---|---:|---:|---:|---:|---:|
| file_parsing | 222.714 | 93.864 | 0.008 | 1.864 | 18.094 |
| donor_bond_detection | 34.914 | 0.000 | 0.000 | 0.000 | 3.083 |
| sterimol | 150.588 | 75.571 | 0.000 | 0.000 | 14.622 |
| buried_volume | 11,716.731 | 0.000 | 0.000 | 0.000 | 1,112.830 |
| pyramidalization | 1.448 | 0.000 | 0.000 | 0.000 | 0.142 |
| conformer_processing | 21.097 | 14.037 | 0.000 | 0.000 | 1.393 |
| model_loading | 0.000 | 0.000 | 0.012 | 0.746 | 0.000 |
| model_inference | 0.000 | 0.000 | 15.451 | 1.152 | 0.000 |
| uncertainty | 0.000 | 0.000 | 0.000 | 83.766 | 0.000 |
| applicability | 0.000 | 0.000 | 0.000 | 0.420 | 0.000 |
| ranking | 0.000 | 0.000 | 0.000 | 0.573 | 0.000 |
| provenance_hashing | 0.000 | 0.000 | 0.000 | 6.465 | 0.000 |
| validation | 0.000 | 0.000 | 2.310 | 0.000 | 0.000 |
| output | 22.995 | 0.541 | 0.135 | 3.827 | 0.000 |
| orchestration | 9.227 | 0.295 | 2.690 | 1,154.630 | 0.000 |
| outside_cli_scope | 6.930 | 3.099 | 2.169 | 5.721 | 0.000 |

Zero means that this stage has no scoped work on that thread for this workload. Conformer-processing rows are exclusive conformer management/aggregation overhead; geometry costs appear in their own stages. Full per-function exclusive/inclusive times, call counts, worker breakdowns, stage ceilings and all ten workload results are in [profile_summary.json](profile_summary.json).

### Probe overhead

| Workload | Diagnostic / native wall |
|---|---:|
| `conformers_56` | 1.165× |
| `db_build_ensembles` | 1.061× |
| `descriptors_10000` | 1.002× |
| `ensemble_sdf` | 1.194× |
| `parse_ensembles_1000` | 1.044× |
| `predict_1000000` | 1.040× |
| `screen_1000` | 1.003× |
| `search_database` | 1.448× |
| `single_large` | 1.839× |
| `single_small` | 1.981× |

The geometry batch and screening traces add about 0.2–0.3%; parsing and inference about 4%. Single-molecule launches add roughly 2.4–2.7 ms, mostly outside the CLI root (including diagnostic JSON writing). Their function rankings remain informative, but their wall fractions and bounds are not reliable estimates of native launch speedups. No startup/report residual is recommended as a single optimization target.

## Allocations and memory

Counters use the existing System allocator without changing allocation requests. Counts include process startup and profiling-session metadata, but exclude construction of the diagnostic report. “Requested” sums allocation sizes and full new reallocation sizes; it is allocation traffic, not resident memory. Peak live requested bytes exclude mappings, stack, static storage and allocator bookkeeping.

| Workload | Allocations | Reallocations | Requested traffic (MiB) | Peak live requested (MiB) |
|---|---:|---:|---:|---:|
| `conformers_56` | 20,710 | 656 | 22.20 | 0.406 |
| `db_build_ensembles` | 21,111 | 773 | 22.16 | 0.424 |
| `descriptors_10000` | 3,589,984 | 104,600 | 3,912.66 | 12.427 |
| `ensemble_sdf` | 4,109 | 951 | 5.07 | 0.408 |
| `parse_ensembles_1000` | 884,204 | 2,897 | 52.64 | 0.157 |
| `predict_1000000` | 593 | 54 | 4.08 | 3.830 |
| `screen_1000` | 384,584 | 10,590 | 436.57 | 3.680 |
| `search_database` | 22,821 | 1,084 | 5.05 | 1.717 |
| `single_large` | 1,135 | 61 | 0.65 | 0.381 |
| `single_small` | 901 | 58 | 0.64 | 0.379 |

All allocation-failure counters are zero. Native RSS and requested heap answer different questions: the prediction workload maps 64 MB of input while its requested live heap is roughly 4 MB. The descriptor batch repeatedly allocates integration grids, producing about 3.82 GiB of requested traffic with about 12.4 MiB peak requested live heap. Allocation volume alone does not establish an allocation-time bottleneck.

### Independent heap attribution

DHAT runs of the frozen symbol binary preserve scientific outputs and corroborate the requested-heap counters. Its whole-process interception includes libc and shutdown; the Rust counters include diagnostic-session metadata and stop before report construction, so exact equality is not expected. DHAT reports lifetimes in instructions, not elapsed time. The global live peak is the sum of each site’s bytes at the same global peak instant, not the sum of independent site maxima. See the [DHAT manual](https://valgrind.org/docs/manual/dh-manual.html).

| Workload | Allocation/reallocation operations | Requested bytes | Global live heap peak (bytes) | Difference vs Rust: operations / requested / peak |
|---|---:|---:|---:|---|
| `conformers_56` | 21,372 | 23,277,755 | 425,394 | +6 / +1,941 / -136 |
| `single_small` | 965 | 673,472 | 396,849 | +6 / +1,942 / -135 |

For 56 conformers, the `integration_grid` allocation site requests 22,020,096 bytes in 56 blocks, 94.6% of total heap traffic. It contributes 393,216 bytes at the global peak. This confirms the grid-allocation source, while the native/diagnostic timing still limits grid-only speedups. Complete allocation stacks and cross-checks are in [heap_summary.json](heap_summary.json).

## Cache and branch simulation

Native `perf stat` and `perf record`, including software events, are denied by this host (`perf_event_paranoid=4`); noninteractive privilege escalation is unavailable. No host security settings were changed. **No native cache-miss or branch-miss measurements are claimed.** Cachegrind 3.26.0 supplies explicitly modeled counts for the exact same ten full workloads, all of which passed their native scientific-output oracle.

The model uses separate 32 KiB, 8-way I1/D1 caches and a 16 MiB, 16-way last-level cache, with 64-byte lines. These dimensions match sysfs, but Cachegrind omits the real L2 stage and modern hardware behavior such as speculative execution and prefetching. Its branch predictor is also a simplified model. Instruction counts are not elapsed time. See the [Cachegrind manual](https://valgrind.org/docs/manual/cg-manual.html).

Automatic detection initially modeled the last-level cache as direct-mapped. That partial run was stopped and explicitly excluded; every accepted case was rerun with the fixed parameters above. Each simulated process starts with its own model state. This is separate from the filesystem warm-cache policy used for native timings.

| Workload | Executed instructions | D1 misses / data references | LL misses / all references | Branch misses | Branch miss rate |
|---|---:|---:|---:|---:|---:|
| `conformers_56` | 1,337,290,690 | 0.2445% | 0.00130% | 2,742,023 | 1.9438% |
| `db_build_ensembles` | 1,337,610,619 | 0.2412% | 0.00158% | 2,749,371 | 1.9486% |
| `descriptors_10000` | 238,829,039,500 | 0.2420% | 0.00008% | 486,714,944 | 1.9324% |
| `ensemble_sdf` | 312,065,580 | 0.2359% | 0.00542% | 619,658 | 1.8787% |
| `parse_ensembles_1000` | 2,060,459,500 | 0.0365% | 0.00062% | 5,627,165 | 3.1120% |
| `predict_1000000` | 41,865,468 | 11.0415% | 3.49690% | 14,171 | 0.3896% |
| `screen_1000` | 23,537,238,010 | 0.3140% | 0.00040% | 51,674,203 | 2.1094% |
| `search_database` | 77,049,356 | 0.5277% | 0.04673% | 379,489 | 3.1881% |
| `single_large` | 35,256,880 | 0.2949% | 0.04045% | 63,289 | 1.6857% |
| `single_small` | 21,676,177 | 0.4638% | 0.06552% | 64,066 | 2.7453% |

D1 denominator = `Dr+Dw`; LL global denominator = `Ir+Dr+Dw`; branch denominator = `Bc+Bi`. The LL local denominator is instead requests that missed L1 and is retained separately in [cache_summary.json](cache_summary.json). All thirteen raw events and every numerator/denominator are preserved there.

Instruction attribution independently concentrates 98.0% of the descriptor batch in the optimized/inlined `compute_from_center` function. The timed nested scopes identify occupancy as the dominant part of that function. Parsing shows 51.7% of instructions in Sterimol projection, 12.9% in whitespace iteration, and 7.0% in decimal parsing. Prediction assigns 21.5% of instructions to `dot_avx2`; the wider Rayon mapping and CLI work account for much of the rest. Instruction percentages must not be substituted for the measured wall fractions in Amdahl calculations.

The simulated descriptor cache-miss ratios are small; prediction touches a much larger record buffer and has larger modeled miss ratios. This does not establish that native execution is compute-, bandwidth-, or branch-limited. No cache-layout or branch rewrite is recommended on those counts alone. Native hardware profiling remains the next measurement needed to choose such an implementation.

## Numerical preservation

Every measured baseline, diagnostic and current-default launch uses the identical manifest. Complete descriptor/search/screen JSON, complete packed records and database artifacts match exactly after removing only explicitly named timing/resource/path metadata. There is no numeric tolerance. Database CSV formatting rounds descriptors, so the companion 56-conformer full JSON comparison supplies the stronger geometry check.

The prediction CLI prints only five preview values and MSE. A separate [full-vector audit](prediction_audit.json) exports every prediction using [profile_predictions.rs](../../examples/profile_predictions.rs): all six baseline/current/diagnostic × one/six-thread outputs contain exactly 4,000,000 identical bytes, SHA-256 `35bc75c2dc0bffa8ba6599b0649479a30db0586f7edc4107f9f7141faa99a987`. The diagnostic audits enable real timing probes and the tracking allocator, check zero dropped scopes, and confirm inference was traced. This supplements rather than replaces the unchanged CLI workload.

## Optimization recommendations supported by this baseline

1. **Investigate a Sterimol-only screening path first.** The current loader computes full descriptors and retains only L/B1/B5 (`src/commands/screen.rs`, `load_library`). In these one-worker runs, buried-volume worker scopes occupy a median 88.12% of process wall, giving a conditional 8.42× ceiling if that work vanished. This worker interval replaces part of the overlapping caller wait; it cannot be added to the 91.31% caller fraction. The full helper also performs geometry/configuration validation and determines exclusions. Preserve those error semantics, donor and axis selection, conformer handling and exact Sterimol operations; simply deleting the call is not equivalent. Verify every screening result and exclusion against the frozen workload before accepting a change.
2. **Concentrate geometry work on `occupied_volumes`.** It takes 92.51% of the 10,000-file diagnostic wall time: a 13.35× elimination ceiling or 1.86× if this function became twice as fast. The complete buried-volume stage is 96.13% (25.82× ceiling). Investigate reducing redundant point/atom work while keeping the exact integration points, floating-point distance comparisons, boundary rules, counts and output operations. No specific SIMD, spatial-filtering or layout implementation has demonstrated a speedup here. Changing grid density, radii, floating-point precision or math order is outside the exact-results constraint.
3. **Treat grid reuse as a limited, secondary experiment.** `integration_grid` repeatedly allocates the same configured points, but accounts for only 3.42% of descriptor-batch wall: at most 1.035× end-to-end even if free. A cache must key all grid-affecting configuration and preserve point order/bit patterns. High allocation traffic is not evidence that it dominates elapsed time.
4. **For ensemble packing, investigate parsing before micro-optimizing arithmetic.** File parsing is 49.89% (2.00× ceiling); Sterimol is 39.31% (1.65× ceiling). Investigate parser allocation/buffer reuse while preserving accepted formats, errors and identical f32 conversion. Reusing already parsed conformers may help repeated files, but these fixtures repeat only 11 ligands; do not project deduplication gains to a unique-chemistry dataset. Preserve conformer weights and accumulation order.
5. **Reuse the Student-t multiplier within one screen.** `student_t_two_sided_quantile` takes 51.376 ms, 4.06% (1.042× ceiling). Inspecting the screen shows an immutable model and fixed confidence/degrees of freedom across candidates. Compute the same function once with the same arguments and return/error semantics; do not approximate its math. This is a smaller opportunity than geometry.
6. **Keep prediction claims tied to the full pipeline.** The batch-inference span accounts for 67.66% (3.09× ceiling), and MSE validation 10.08% (1.112×). The inference span includes allocation, Rayon setup, data access/page faults and result production. These traces do not justify attributing all of it to dot products. Six threads already deliver a measured 2.21× end-to-end improvement on this fixed workload. Inspect kernel/data-access costs on a machine with hardware profiling before choosing a SIMD or memory-layout rewrite. Preserve the MSE and complete prediction output.

Output, donor/bond detection and pyramidalization are low-priority targets in the descriptor batch (0.19%, 0.29%, 0.01% respectively). Search output is 19.91% of its diagnostic launch because this workload emits all 1,541 hits; its 1.249× ceiling does not apply directly to a normal top-10 search. The small search trace also has substantial probe/report overhead.

**Acceptance rule for any future optimization:** change one hypothesis at a time; keep the frozen input manifest, affinity, environment, compiler and output mode; repeat the same native workload series with `--compare-to baseline_final_t1`; inspect median, spread, CPU and RSS; compare all packed bytes and full-precision scientific outputs, plus every prediction f32. Add fixtures for changed error paths. Reject numerical differences or skipped records. Repeat diagnostics to confirm the intended cost fell, and report the measured end-to-end result even when it disagrees with an ideal bound. No optimization has been implemented or credited with a speedup in this report.

## Measurement controls and reproducibility

Inputs are hash-checked before every series and one initial launch is retained separately. Output goes to regular files without forcing storage durability. These are warm-cache repeated process launches; no cache flushing or cold-start storage claim is made. Workloads run sequentially, with no builds or profilers overlapping accepted native measurements. Boost remains enabled, the CPU governor is `powersave` under `amd-pstate-epp`, and ASLR remains enabled. Those settings and machine load limit cross-machine extrapolation.

Early runs using Python as the direct parent inherited Python’s RSS high-water mark across exec and are explicitly excluded. The accepted runs use [profile_child.c](../../scripts/profile_child.c), a small native parent with exact-child `wait4`. A touched 32 MiB allocation agreed with GNU time within 0.4%; calibration evidence is retained. No earlier partial run is used to select favorable timings.

The tools reject missing reports, dropped scopes, changed inputs, changed scientific outputs, and inconsistent native/diagnostic execution settings. Tests cover resource-launch accounting, numerical-output acceptance, stage partitioning, worker overlap, zero-filled absent samples, Amdahl bounds, simulated-event denominators and DHAT global-peak accounting.

Raw inputs, stdout/stderr, commands, metrics and stage reports remain under `.stericx/profiling/`. Run labels are immutable and must be new for each series. The [native evidence archive](native_evidence.json.gz) retains 465 JSON records: machine/compiler identity, accepted per-launch native measurements, diagnostic traces, RSS calibration, and excluded-run metadata. The [profiler evidence archive](profiler_evidence.json.gz) retains original event files, annotations, allocation profiles, commands, summaries and per-file SHA-256 hashes, including the excluded automatic-cache-detection run’s metadata.

```bash
# Freeze baseline binary before changing performance-sensitive source.
cargo build --release --locked
mkdir -p .stericx/profiling/baseline/bin
cp target/release/stericx .stericx/profiling/baseline/bin/stericx
python3 scripts/profile_stericx.py prepare
python3 scripts/profile_stericx.py run --binary .stericx/profiling/baseline/bin/stericx --label baseline_new --reps 7 --threads 1 --affinity 2

# Diagnostic build is separate from the native release used for performance.
CARGO_PROFILE_RELEASE_DEBUG=1 cargo build --release --locked --features profiling --target-dir .stericx/profiling/instrument-target
python3 scripts/profile_stericx.py run --binary .stericx/profiling/instrument-target/release/stericx --label diagnostic_new --reps 5 --threads 1 --affinity 2 --profile-env STERICX_PROFILE_PATH --compare-to baseline_new
python3 scripts/summarize_profile.py --baseline baseline_new --diagnostic diagnostic_new --output .stericx/profiling/new_report
```

The recorded baseline binary hash and compiler identity are in [baseline_identity.json](baseline_identity.json); [source_identity.json](source_identity.json) records the pre-optimization source/tool snapshot. Reproducing the original source requires commit `37a3f0a` and its locked dependencies; the symbol build was compiled from a frozen archive of that commit. The Valgrind runner defaults to that separate symbol binary and preserves exact application commands. Tool path setup and simulation parameters are recorded in its evidence; profiler wall times must not be treated as native latency.

With the recorded Valgrind installation and frozen symbol binary available, reproduce the diagnostic tools separately from native timing:

```bash
python3 scripts/profile_valgrind.py --label cache_new --tool cachegrind --compare-to baseline_final_t1 --i1 32768,8,64 --d1 32768,8,64 --ll 16777216,16,64
python3 scripts/profile_valgrind.py --label heap_new --tool dhat --compare-to baseline_final_t1 --allocator-reference instrumented_final_t1
```

`--capabilities` can point to a JSON file describing a local Valgrind installation; the required paths are `valgrind_binary`, `VALGRIND_LIB`, and `cg_annotate`. Use a no-space path for `VALGRIND_LIB` (a symlink is sufficient); on this machine the distro shell wrapper also required invoking `valgrind.bin` directly. The recorded installation is local to `.stericx/profiling/tools` and exposed through `/tmp/stericx-valgrind-1000`; recreate that symlink after reboot if needed. Build the frozen baseline source separately with `CARGO_PROFILE_RELEASE_DEBUG=1 cargo build --release --locked` for optimized symbols, or supply its path with `--binary`. Defaults preserve the reported cache model; change them deliberately only when establishing another comparison. The default timeout is unlimited, and an explicit timeout fails the run rather than shrinking its workload.
