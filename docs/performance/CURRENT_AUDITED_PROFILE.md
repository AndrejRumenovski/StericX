# Fresh current-audited performance baseline

This is a new profile of `b515c4f`, **not an accepted optimization**. The scientific baseline prerequisite remains unresolved; see [baseline status](BASELINE_STATUS.md).

Ryzen 5 5600G; Linux 7.0.0-31; Rust 1.97.0; locked native release build. Each configuration has one warmup and seven measured launches. CPU affinity is CPU2 at one thread, CPUs2/3 at two, CPUs0–3 at four, and CPUs0–5 at six. Builds and profilers did not overlap native or diagnostic timing runs. Warm filesystem cache; process startup and regular-file output included. No samples excluded.

MAD is median absolute deviation, not a confidence interval. CPU100% means one core. RSS comes from the exact child’s wait4 high-water mark. Full CPU user/system time, page faults, context switches, distributions, fingerprints and diagnostic allocation counts are retained in [the machine-readable measurements](current_audited_measurements.json).

| Workload | Threads | Wall median ± MAD (ms) | Full range (ms) | CPU (%) | RSS median (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `conformers_56` | 1 | 30.084 ± 0.086 | 29.734–30.468 | 99.2 | 4.93 |
| `conformers_56` | 2 | 29.443 ± 0.175 | 29.203–29.907 | 99.3 | 5.03 |
| `conformers_56` | 4 | 29.763 ± 0.306 | 29.256–30.662 | 99.1 | 5.01 |
| `conformers_56` | 6 | 30.285 ± 1.133 | 29.151–50.480 | 99.3 | 4.92 |
| `db_build_ensembles` | 1 | 31.069 ± 0.120 | 30.597–32.247 | 98.9 | 5.36 |
| `db_build_ensembles` | 2 | 17.027 ± 0.448 | 16.456–18.109 | 182.2 | 5.64 |
| `db_build_ensembles` | 4 | 10.115 ± 0.102 | 9.930–10.329 | 321.2 | 6.44 |
| `db_build_ensembles` | 6 | 7.693 ± 0.135 | 7.461–9.186 | 441.9 | 7.00 |
| `descriptors_10000` | 1 | 5044.513 ± 25.539 | 5004.087–5082.803 | 99.8 | 18.13 |
| `descriptors_10000` | 2 | 5030.745 ± 2.950 | 5022.531–5061.925 | 99.8 | 18.19 |
| `descriptors_10000` | 4 | 4992.289 ± 15.158 | 4961.564–5065.330 | 99.7 | 18.04 |
| `descriptors_10000` | 6 | 5008.436 ± 20.082 | 4974.816–5081.664 | 99.7 | 17.94 |
| `ensemble_sdf` | 1 | 7.666 ± 0.089 | 7.416–8.037 | 97.1 | 5.01 |
| `ensemble_sdf` | 2 | 7.314 ± 0.041 | 7.233–7.902 | 98.5 | 5.00 |
| `ensemble_sdf` | 4 | 7.679 ± 0.259 | 7.259–9.877 | 98.4 | 4.87 |
| `ensemble_sdf` | 6 | 7.665 ± 0.205 | 7.247–7.934 | 98.3 | 4.95 |
| `parse_ensembles_1000` | 1 | 180.449 ± 0.172 | 180.145–180.655 | 99.8 | 4.75 |
| `parse_ensembles_1000` | 2 | 180.785 ± 1.163 | 179.623–189.060 | 99.9 | 4.73 |
| `parse_ensembles_1000` | 4 | 179.424 ± 0.361 | 178.290–180.322 | 99.9 | 4.75 |
| `parse_ensembles_1000` | 6 | 181.072 ± 1.874 | 178.806–184.375 | 99.9 | 4.71 |
| `predict_1000000` | 1 | 20.824 ± 0.066 | 20.616–20.891 | 99.3 | 69.50 |
| `predict_1000000` | 2 | 13.861 ± 0.128 | 13.422–13.996 | 154.6 | 69.53 |
| `predict_1000000` | 4 | 10.096 ± 0.144 | 9.942–10.335 | 218.6 | 69.30 |
| `predict_1000000` | 6 | 8.813 ± 0.205 | 8.596–10.190 | 263.9 | 69.25 |
| `screen_1000` | 1 | 323.953 ± 3.321 | 316.131–335.039 | 99.4 | 8.84 |
| `screen_1000` | 2 | 218.620 ± 1.153 | 217.467–227.467 | 150.1 | 9.20 |
| `screen_1000` | 4 | 155.138 ± 0.681 | 154.457–178.722 | 207.6 | 9.89 |
| `screen_1000` | 6 | 138.913 ± 1.812 | 136.291–153.339 | 234.3 | 10.67 |
| `search_database` | 1 | 6.240 ± 0.181 | 6.056–6.519 | 95.5 | 6.75 |
| `search_database` | 2 | 6.889 ± 0.660 | 6.089–16.065 | 95.4 | 6.67 |
| `search_database` | 4 | 6.222 ± 0.084 | 6.007–6.311 | 96.0 | 6.71 |
| `search_database` | 6 | 6.399 ± 0.058 | 6.253–6.723 | 96.1 | 6.74 |
| `single_large` | 1 | 1.869 ± 0.046 | 1.772–2.026 | 94.0 | 4.84 |
| `single_large` | 2 | 1.893 ± 0.099 | 1.778–2.364 | 94.7 | 4.88 |
| `single_large` | 4 | 1.911 ± 0.054 | 1.791–2.035 | 94.2 | 4.85 |
| `single_large` | 6 | 2.011 ± 0.139 | 1.844–2.151 | 94.2 | 4.89 |
| `single_small` | 1 | 1.907 ± 0.054 | 1.813–2.004 | 93.2 | 4.83 |
| `single_small` | 2 | 1.753 ± 0.015 | 1.736–1.793 | 94.6 | 4.86 |
| `single_small` | 4 | 1.819 ± 0.099 | 1.720–2.433 | 94.4 | 4.81 |
| `single_small` | 6 | 1.967 ± 0.058 | 1.795–2.034 | 94.4 | 4.87 |

All 40 native and 40 diagnostic configurations completed: 280 measured native launches, 200 diagnostic launches and 80 warmups. Every diagnostic scientific fingerprint matches the corresponding native thread setting. All eight native/diagnostic ×1/2/4/6-thread prediction exports contain identical 4,000,000-byte vectors, SHA-256 `35bc75c2dc0bffa8ba6599b0649479a30db0586f7edc4107f9f7141faa99a987`.

## Fresh attribution and candidates

At one thread, buried-volume occupancy takes 81.753% of descriptor-batch diagnostic wall; grid preparation takes 9.134%. Screening recomputes an invariant Student-t quantile 1,000 times, accounting for 15.651% of wall at one thread and 33.60% at six. Exact candidate hypotheses, ceilings, predicted savings and rejection gates are recorded before implementation in [candidate admission](CANDIDATE_ADMISSION.md).

Worker elapsed spans overlap caller waiting. They are retained separately and must not be added to main-thread wall shares. Diagnostic overhead reaches 1.219× for six-thread prediction and 1.105× for six-thread screening; diagnostic percentages are conditional prioritization evidence, not native speedup measurements.

The benchmark corpus deliberately repeats geometry: `descriptors_10000` cycles 56 conformers, `screen_1000` cycles 11 reaction rows, and `predict_1000000` cycles 11 records. No extrapolation to 10,000 distinct chemistries or experimental predictive validity is made.

Native hardware counters are denied by this host (`perf_event_paranoid=4`), including software task-clock. Separate Cachegrind/DHAT investigations supply modeled instruction/cache/branch and allocation evidence; they do not supply native hardware events or substitute for native timings.

## Independent counter and heap attribution

All ten complete one-thread Cachegrind workloads reproduced their native scientific
fingerprints. All 13 raw event fields were independently summed and checked against
annotation program totals. The explicit cache model uses 32 KiB, 8-way I1/D1 and
16 MiB, 16-way last-level caches, all with 64-byte lines. It omits real hardware
prefetching, speculation and the separate L2 stage. See the complete
[counter summary](current_audited_cachegrind.json).

| Workload | Instruction events | Leading self-cost shares |
| --- | ---: | --- |
| `descriptors_10000` | 68,100,568,139 | Occupancy 82.276%; grid 10.664% |
| `parse_ensembles_1000` | 2,056,089,986 | Sterimol 51.765%; whitespace parsing 12.963% |
| `predict_1000000` | 41,867,991 | Rayon mapping 52.548%; dot kernel 21.496% |
| `screen_1000` | 3,900,467,943 | Occupancy 45.361%; grid 18.618%; bootstrap sorting 8.224% |

These are instruction shares, not wall-time shares. Inlining and library costs
make their categories different from the lexical diagnostic timing scopes.

DHAT's complete 56-conformer workload requests 23,340,655 bytes and reaches a
425,706-byte live heap peak. The integration grid accounts for 22,020,096 requested
bytes in 56 allocations (94.34%). Differences from Rust counters are six
allocation operations, 1,909 requested bytes and minus 161 peak bytes, reflecting
their different interception/reporting scopes. Both this workload and
`single_small` reproduce native outputs. The [heap summary](current_audited_heap.json)
retains the complete stacks and cross-checks. Allocation traffic is not RSS or
proof of allocation-time dominance.

The separate counter seal is
`f2fd2f2c17f1fef3808ab558d02c00c3b821bd73bf53e85f4ab90a9cc40703c0`, binding 13 local
files and 97 raw profiler artifacts under the supplementary
`current_audited_b515c4f_counters` evidence directory. These runs did not change
the native measurement seal.

## Reproduction and evidence

`scripts/performance_freeze.py` records source, toolchain, CPU/OS, input and audit identities, copies fresh native/diagnostic executables, and runs the existing profiling harness under the recorded controls. New immutable labels use prefix `current_audited_b515c4f_20260921`. The fixed input manifest SHA-256 is `5f2a35784d64c12122a886db1d6bdce1e2d7b012cece31256d22b175f7d1b20c`.

The local evidence root is `.stericx/profiling/scientifically_validated_optimization/current_audited_b515c4f/`. Its measurement seal binds 201 snapshot files and 2,656 raw run artifacts; SHA-256 `75849d880c9a3e9e5965f59260bdd40143f8706e0c8d7e645fc830a5282f8aa8`. Source/executable hashes were checked again after measurement. Historical audit evidence remains unchanged.
