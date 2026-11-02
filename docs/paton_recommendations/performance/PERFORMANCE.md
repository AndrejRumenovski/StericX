# Performance after the Paton descriptor extensions

The confirmed quiet-window run on **2026-11-02, 15:38–15:42 UTC** found no
median slowdown above the prospectively stated **10% investigation threshold**
in any of 28 workload/thread configurations. The largest ratio-of-medians
slowdown was **2.2505%** (six-thread database construction); its paired slowdown
was 0.8431%. All **504 attempted samples** completed with exact scientific output
fingerprints. No production numerical change was made to obtain this result.

This supports performance preservation for the measured workloads on this host.
It does not establish a universal speedup, a confidence bound on every workload,
or a global scientific pass. In particular, the independently documented strict
DBSTEP hard-sphere volume failure remains unchanged.

## Confirmed native results

Times include process startup, calculation and regular-file output on the
original workspace filesystem. Values below are medians of eight measurements
per binary; the full table also includes MAD, ranges, paired ratios, CPU time,
RSS and throughput for every configuration.

| Workload | Threads | Baseline (ms) | Final (ms) | Change |
|---|---:|---:|---:|---:|
| Small molecule | 1 | 1.943 | 1.925 | −0.9% |
| Large molecule | 1 | 1.905 | 1.899 | −0.3% |
| One 86-conformer SDF | 1 | 47.355 | 46.912 | −0.9% |
| 56 conformer files | 1 | 29.393 | 29.264 | −0.4% |
| 56 conformer files | 6 | 7.348 | 7.132 | −2.9% |
| 10,000 descriptor files | 1 | 5,029.593 | 5,055.721 | +0.5% |
| 10,000 descriptor files | 6 | 917.490 | 923.776 | +0.7% |
| 1,000 ensemble parses/packed exports | 6 | 186.159 | 187.735 | +0.8% |
| Database: 56 geometries, 11 ligands | 6 | 7.968 | 8.147 | +2.3% |

Final batch throughput is **1,978 files/s** at one thread and **10,825 files/s**
at six threads (baseline 1,988 and 10,899 files/s). This frozen workload repeats
the retained geometry corpus; it is not 10,000 independent chemical validation
cases. The ensemble-SDF throughput unit is one completed ensemble, whereas the
database throughput unit is one ligand.

Median peak child RSS for the 10,000-file workload is **21.023 MiB** final versus
20.994 MiB baseline at one thread and **26.332 MiB** versus 26.447 MiB at six.
Six-thread summed CPU time is 5.325 s versus 5.329 s. A final six-thread sample
took 1.343 s despite its 0.924 s median; it remains in the raw data and reported
range. Quiet-window measurements still contain scheduler and filesystem noise.

See [all confirmed timings](CONFIRMED_TIMING_TABLE.md),
[the native journal](native_confirmed_summary.json), and
[the complete campaign review](confirmed_review.json).

## Protocol and identities

The [protocol](PROTOCOL.md) was frozen before the first timing campaign.
It uses the existing native `wait4` launcher, immutable corrected workload
manifest, warm-cache policy, one warmup per binary and eight alternating AB/BA
pairs. Thread counts 1, 2, 4 and 6 use CPU affinities `2`, `2,3`, `0,1,2,3`, and
`0,1,2,3,4,5`. There is no cache flushing, outlier removal or sample pooling.
Both a ratio of side medians above 1.1 and a median within-pair baseline/final
ratio below 1/1.1 trigger investigation.

The host is Linux x86_64, AMD Ryzen 5 5600G (six physical cores, twelve logical
CPUs), kernel `7.0.0-31-generic`. Full CPU and controlled-environment metadata are
in the native journal. The frozen source builds use Rust/Cargo 1.97.0. The task
baseline is `c9380d2c24ba5d48cd52ec13c139ee6e57f04991`, not an earlier shipped
binary or a different optimization baseline.

| Executable | SHA-256 |
|---|---|
| Baseline native | `8ed1480f3758e6cf5320549029cf75e10dfd7754e6e9ba6184bb24f8dec95447` |
| Final native | `3b9079338386a2532b17ebe844a28d72344a7bdf91d83526d4b3d0fd00438424` |

All 53 current production source/Cargo files match the final scientific build.
Timing admission binds the byte-exact 128,021-record scientific replay, reviewed
independent references, and successful three-OS engineering CI with matching
source hashes. Executables and frozen input inventories were rehashed before
and after timing. See [scientific validation](../SCIENTIFIC_VALIDATION.md).

## Earlier attempts and storage investigation

Every earlier campaign is retained separately. Before the final run, the user
clarified that the initial idle-window request had been misunderstood. The two
directories containing `idle` therefore do **not** certify an otherwise idle
machine. Their names and unfavorable results are retained without reinterpretation.

| Campaign, in UTC order | Output storage | Attempts | Configurations triggering investigation |
|---|---|---:|---:|
| `performance_v1` | Original filesystem | 504 | 4 |
| `performance_v2_tmpfs` | RAM-backed regular files | 504 | 0 |
| `performance_idle_v1` | Original filesystem | 504 | 7 |
| `performance_idle_tmpfs_v1` | RAM-backed regular files | 504 | 2 |
| `performance_confirmed_20261102` | Original filesystem, confirmed quiet window | 504 | 0 |

All **2,520 native attempts**, including warmups, remain available. The report
controller rehashes their metrics, recomputes every raw scientific fingerprint
and reconstructs all 140 configuration summaries exactly. The earlier
`review.json` and `TIMING_TABLE.md` describe only the first storage repeat;
`confirmed_review.json` and `CONFIRMED_TIMING_TABLE.md` are authoritative for the
final conclusion. No earlier favorable run replaces an unfavorable one.

Sixteen separate intrusive `strace` diagnostics found short-process delays in
the same 678-byte output `write`: up to 84.160 ms for baseline and 14.711 ms for
final in those diagnostic samples. They are not timing samples. The
[prospective storage study](OUTPUT_STORAGE_STUDY.md) changed only output storage
to tmpfs. It removed the first campaign's triggers, but the subsequent tmpfs
campaign still showed two triggers, including 24.5% slower side medians for
four-thread database construction. Its summed CPU difference was approximately
3%, with increased scheduling variability. These effects did not recur in the
confirmed run. The evidence is consistent with host/storage contention and does
not identify a reproducible engine regression.

Confirmed-window Linux pressure snapshots show CPU `some avg10` at 0.00% before
and after, and I/O `some avg10` at 0.54% before and 0.45% after. Longer windows
retain earlier activity. This is supporting context, not proof of perfect host
isolation. Host snapshots, all traces, and byte-exact relocation receipts for
both tmpfs campaigns are retained and hashed.

## Allocations and new ensemble API

The existing tracking allocator was built separately into each frozen source
tree. All **14 baseline/final allocation comparisons** (seven workloads at one
and six threads; 28 instrumented runs) have identical allocation and reallocation
counts, zero failed allocations and exact scientific output fingerprints.
Requested-byte differences are −2, +1 or +8 bytes, consistent with differing
command/output path lengths; they are not a reduction in scientific work.
Peak live requested-byte differences are at most two bytes upward; two parallel
workloads are lower by 1,242/2,303 bytes, consistent with allocation lifetime
ordering. Instrumented timings are excluded from throughput results.

The new structured six-conformer descriptor/statistics call has no previous
counterpart. One warmed native call, including destruction of its results,
performs **871 allocations and 73 reallocations**, requesting **2,454,897 bytes**
in total with zero failed calls and zero net live-byte change. Process-wide peak
live requested bytes through the snapshot are 446,776; instrumented peak RSS is
3,772,416 bytes. Parsing and JSON encoding are excluded from this diagnostic.
These quantities measure different things and must not be added together.
See [the allocation receipt](ensemble_allocations.json).

The installed Python 3.12.13 package is `stericx 0.4.0.dev0`. Inputs are the same
six retained ligand-401 conformers and explicit options. Parsed inputs are
reused; one warmup precedes eight alternating blocks of 100 calls per layer on
CPU 2. Scientific results agree exactly before timing.

| Layer | Median per call (ms) | MAD (ms) | Range (ms) |
|---|---:|---:|---:|
| Native descriptors and statistics | 3.0183 | 0.0291 | 2.9556–3.1264 |
| Python private native JSON bridge | 3.3409 | 0.0863 | 3.2250–3.4368 |
| Public typed Python results | 4.0959 | 0.0866 | 3.9914–4.2431 |

The public-minus-bridge difference is approximately **0.7550 ms**, covering
Python options encoding, JSON parsing and typed-object construction in the same
installed wheel. Native executable and portable wheel builds have different
target settings, so their difference is not a pure FFI overhead measurement.
These are absolute new-feature costs, not a historical speedup. Weighting,
sensitivity CSV, input parsing and density analysis are outside this timed call.
See [all Python blocks and identities](python_confirmed_summary.json).

## Reproduction and evidence

From the retained workspace, choose **new, empty output paths** for each run:

```bash
.venv/bin/python scripts/benchmark_paton_extensions.py \
  --baseline .stericx/paton/baseline/build \
  --candidate .stericx/paton/final_build \
  --workloads .stericx/profiling/scientifically_validated_optimization/corrected_workloads_v2/workloads/manifest.json \
  --review .stericx/paton/final_scientific_review_v2/review.json \
  --engineering docs/paton_recommendations/validation/ci_engineering_initial.json \
  --output .stericx/paton/performance_new_run

.stericx/paton/python-build/bin/python scripts/benchmark_paton_python.py \
  --native target/release/examples/paton_ensemble_benchmark \
  --input python/examples/data/ligand_401_ensemble.xyz \
  --options python/examples/data/ligand_401_options.json \
  --iterations 100 --output .stericx/paton/python_benchmark_new_run

.venv/bin/python scripts/report_paton_performance.py \
  --output .stericx/paton/performance_review_new_run
```

The last command checks the five retained campaigns and regenerates the confirmed
tables; it does not rerun timing. The benchmark's `--allocation-baseline` and
`--allocation-candidate` options accept the separately built instrumented
executables recorded in `allocation_baseline_receipt.json` and
`allocation_final_receipt.json`. The standalone ensemble diagnostic is built
with `cargo build --release --features profiling --example paton_ensemble_allocations`.

[The final evidence receipt](confirmed_evidence_receipt.json) binds summaries,
raw confirmed-run files, allocation results, host observations and controller
sources. Older receipts and all interrupted/failed attempts remain. Large raw
campaigns and immutable input/build trees are local under `.stericx/`; committed
receipts identify them but do not imply a fresh clone contains every raw artifact.
No release or publication is admitted by this performance report.
