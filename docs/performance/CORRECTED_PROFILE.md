# Corrected baseline profiling

The corrected baseline has completed all 40 native and 40 diagnostic configurations. No optimization has been accepted and no baseline-to-candidate speedup is claimed. These measurements use commit `1ecfbd5be8b354711bdadd8a8d148417e03d8e7a`, the reviewed `accurate_baseline_v1` freeze and fresh `corrected_workloads_v2` inputs. Earlier audited-build measurements are separate historical evidence.

Native release runs have one warmup and seven measured launches per configuration; diagnostic release runs have one warmup and five. Threads/CPU affinities are 1→2, 2→2,3, 4→0–3, 6→0–5. Runs are sequential, with warm filesystem cache, process startup and regular-file output included. No measured sample is excluded. The source, executable, environment, machine and input hashes are retained in [the complete machine-readable measurements](current_corrected_measurements.json). MAD is median absolute deviation, not a confidence interval. CPU 100% means one core. RSS is the exact child's wait4 high-water mark.

The original complete one-thread native series is retained. An initial two-thread orchestrator incorrectly compared its full fingerprint with the one-thread series; the printed `rayon_threads` metadata differed. That partial failure is retained in `accurate_baseline_native_measurements_v2/initial_harness_failure.json`. New complete 2/4/6-thread labels were collected without that invalid cross-thread metadata comparison. This was a harness configuration failure, not a scientific mismatch or a tolerance change.

Every diagnostic fingerprint matches its native counterpart at the same thread setting. All eight native/diagnostic ×1/2/4/6 full prediction exports contain identical 4,000,000-byte vectors, SHA-256 `4afd2f0be859ebd8785748369fa2e157d5b38b78fa98a122dd6a55b6cdb11d97`. The rebuilt native CLI is byte-identical to the admitted native executable; headline measurements use the original admitted executable.

## Native distributions

User/system CPU time, page faults, context switches, p90 and all distribution statistics are retained for every row in the JSON. The table reports milliseconds, CPU percent, MiB and median minor/major page faults.

| Workload | Threads | Median ms | MAD ms | Range ms | CPU % | RSS MiB | Faults minor/major |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `conformers_56` | 1 | 29.911 | 0.142 | 29.466–30.417 | 99.2 | 5.14 | 307/0 |
| `conformers_56` | 2 | 29.539 | 0.292 | 29.112–32.552 | 99.3 | 5.11 | 307/0 |
| `conformers_56` | 4 | 30.453 | 0.316 | 29.312–34.116 | 99.3 | 5.09 | 307/0 |
| `conformers_56` | 6 | 30.037 | 0.089 | 29.222–31.658 | 99.3 | 5.11 | 306/0 |
| `db_build_ensembles` | 1 | 30.487 | 0.203 | 30.230–31.052 | 99.3 | 5.57 | 362/0 |
| `db_build_ensembles` | 2 | 17.078 | 0.054 | 16.891–17.180 | 182.1 | 5.88 | 512/0 |
| `db_build_ensembles` | 4 | 10.029 | 0.250 | 9.711–12.887 | 325.1 | 6.51 | 810/0 |
| `db_build_ensembles` | 6 | 7.623 | 0.263 | 7.360–8.349 | 444.4 | 7.38 | 1122/0 |
| `descriptors_10000` | 1 | 4991.033 | 13.970 | 4974.097–5033.607 | 99.7 | 21.01 | 4632/0 |
| `descriptors_10000` | 2 | 5133.013 | 71.448 | 4984.977–5213.349 | 99.8 | 20.88 | 4632/0 |
| `descriptors_10000` | 4 | 4972.785 | 18.980 | 4937.008–5019.900 | 99.7 | 20.89 | 4633/0 |
| `descriptors_10000` | 6 | 4958.143 | 4.935 | 4937.876–4963.078 | 99.7 | 20.87 | 4633/0 |
| `ensemble_sdf` | 1 | 47.033 | 0.125 | 46.719–48.267 | 99.6 | 5.89 | 597/0 |
| `ensemble_sdf` | 2 | 46.834 | 0.088 | 46.531–47.309 | 99.7 | 5.83 | 596/0 |
| `ensemble_sdf` | 4 | 46.825 | 0.057 | 46.513–47.357 | 99.7 | 5.85 | 596/0 |
| `ensemble_sdf` | 6 | 47.621 | 0.753 | 46.617–66.852 | 99.4 | 5.77 | 596/0 |
| `parse_ensembles_1000` | 1 | 184.947 | 0.724 | 183.583–187.577 | 99.8 | 4.89 | 223/0 |
| `parse_ensembles_1000` | 2 | 183.160 | 0.455 | 182.581–187.744 | 99.9 | 5.02 | 224/0 |
| `parse_ensembles_1000` | 4 | 186.805 | 2.741 | 180.007–189.546 | 99.9 | 4.85 | 223/0 |
| `parse_ensembles_1000` | 6 | 184.686 | 1.534 | 182.108–186.737 | 99.9 | 4.74 | 221/0 |
| `predict_1000000` | 1 | 21.970 | 0.132 | 21.317–22.363 | 99.3 | 69.56 | 2157/0 |
| `predict_1000000` | 2 | 14.479 | 0.168 | 14.304–16.076 | 153.2 | 69.55 | 2161/0 |
| `predict_1000000` | 4 | 10.926 | 0.210 | 10.663–11.155 | 214.7 | 69.40 | 2182/0 |
| `predict_1000000` | 6 | 10.350 | 0.565 | 9.172–11.859 | 247.4 | 69.33 | 2202/0 |
| `screen_1000` | 1 | 256.825 | 5.013 | 249.687–269.922 | 97.2 | 10.22 | 2309/0 |
| `screen_1000` | 2 | 254.398 | 3.333 | 249.669–265.958 | 99.9 | 10.24 | 2309/0 |
| `screen_1000` | 4 | 259.134 | 3.096 | 247.381–262.230 | 96.9 | 10.18 | 2309/0 |
| `screen_1000` | 6 | 252.135 | 3.442 | 246.749–260.610 | 99.9 | 10.06 | 2309/0 |
| `search_database` | 1 | 7.235 | 0.087 | 7.041–7.949 | 95.7 | 7.80 | 969/0 |
| `search_database` | 2 | 7.054 | 0.108 | 6.945–21.810 | 95.9 | 7.81 | 969/0 |
| `search_database` | 4 | 7.076 | 0.168 | 6.895–7.924 | 96.2 | 7.76 | 969/0 |
| `search_database` | 6 | 7.205 | 0.104 | 6.983–7.313 | 96.2 | 7.70 | 968/0 |
| `single_large` | 1 | 1.845 | 0.041 | 1.804–1.983 | 94.0 | 5.06 | 246/0 |
| `single_large` | 2 | 2.010 | 0.103 | 1.821–2.222 | 94.7 | 5.03 | 246/0 |
| `single_large` | 4 | 2.083 | 0.081 | 1.876–2.193 | 91.6 | 5.03 | 245/0 |
| `single_large` | 6 | 2.044 | 0.091 | 1.841–2.257 | 94.0 | 5.06 | 245/0 |
| `single_small` | 1 | 2.084 | 0.153 | 1.918–2.648 | 93.2 | 5.00 | 245/0 |
| `single_small` | 2 | 1.884 | 0.037 | 1.847–1.983 | 91.2 | 5.02 | 245/0 |
| `single_small` | 4 | 2.030 | 0.082 | 1.925–2.136 | 91.8 | 5.12 | 245/0 |
| `single_small` | 6 | 2.337 | 0.285 | 1.922–2.622 | 92.9 | 5.02 | 245/0 |

## Diagnostic components and conditional ceilings

All main-thread exclusive stage/function profiles, call counts, worker elapsed totals and allocation medians are in the JSON. Exclusive main-thread times partition the CLI root; outside-root startup/exit/report time is shown separately. Worker elapsed spans overlap caller waiting and must not be added to main-thread wall shares. Shares below use external diagnostic process wall time. An infinite local improvement has Amdahl ceiling 1/(1−f); a 2× local improvement gives 1/(1−f+f/2). These are conditional bounds, not measured native speedups. Nested/inclusive rows overlap.

| Workload / component | Scope | Threads | Wall share % | Calls | Infinite ceiling | 2× local ceiling |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `descriptors_10000` / `descriptors::descriptors_for_file` | inclusive | 1 | 99.161 | 10000 | 119.215× | 1.983× |
| `descriptors_10000` / `descriptors::descriptors_for_file` | inclusive | 2 | 99.210 | 10000 | 126.659× | 1.984× |
| `descriptors_10000` / `descriptors::descriptors_for_file` | inclusive | 4 | 99.212 | 10000 | 126.906× | 1.984× |
| `descriptors_10000` / `descriptors::descriptors_for_file` | inclusive | 6 | 99.267 | 10000 | 136.503× | 1.985× |
| `descriptors_10000` / `buried_volume::occupied_volumes` | self | 1 | 81.773 | 30000 | 5.486× | 1.692× |
| `descriptors_10000` / `buried_volume::occupied_volumes` | self | 2 | 81.774 | 30000 | 5.487× | 1.692× |
| `descriptors_10000` / `buried_volume::occupied_volumes` | self | 4 | 81.920 | 30000 | 5.531× | 1.694× |
| `descriptors_10000` / `buried_volume::occupied_volumes` | self | 6 | 82.054 | 30000 | 5.572× | 1.696× |
| `descriptors_10000` / `buried_volume::integration_grid` | self | 1 | 9.030 | 10000 | 1.099× | 1.047× |
| `descriptors_10000` / `buried_volume::integration_grid` | self | 2 | 9.030 | 10000 | 1.099× | 1.047× |
| `descriptors_10000` / `buried_volume::integration_grid` | self | 4 | 8.943 | 10000 | 1.098× | 1.047× |
| `descriptors_10000` / `buried_volume::integration_grid` | self | 6 | 8.925 | 10000 | 1.098× | 1.047× |
| `screen_1000` / `model::student_t_two_sided_quantile` | self | 1 | 20.035 | 1000 | 1.251× | 1.111× |
| `screen_1000` / `model::student_t_two_sided_quantile` | self | 2 | 20.286 | 1000 | 1.254× | 1.113× |
| `screen_1000` / `model::student_t_two_sided_quantile` | self | 4 | 20.052 | 1000 | 1.251× | 1.111× |
| `screen_1000` / `model::student_t_two_sided_quantile` | self | 6 | 20.354 | 1000 | 1.256× | 1.113× |
| `screen_1000` / `sterimol::params_from_projection` | self | 1 | 27.359 | 5085 | 1.377× | 1.158× |
| `screen_1000` / `sterimol::params_from_projection` | self | 2 | 27.268 | 5085 | 1.375× | 1.158× |
| `screen_1000` / `sterimol::params_from_projection` | self | 4 | 27.264 | 5085 | 1.375× | 1.158× |
| `screen_1000` / `sterimol::params_from_projection` | self | 6 | 27.406 | 5085 | 1.378× | 1.159× |
| `screen_1000` / `parse_coordinate_file` | inclusive | 1 | 30.237 | 5085 | 1.433× | 1.178× |
| `screen_1000` / `parse_coordinate_file` | inclusive | 2 | 30.384 | 5085 | 1.436× | 1.179× |
| `screen_1000` / `parse_coordinate_file` | inclusive | 4 | 30.353 | 5085 | 1.436× | 1.179× |
| `screen_1000` / `parse_coordinate_file` | inclusive | 6 | 30.367 | 5085 | 1.436× | 1.179× |
| `parse_ensembles_1000` / `commands::parse::conformer_row` | inclusive | 1 | 96.554 | 1000 | 29.023× | 1.933× |
| `parse_ensembles_1000` / `commands::parse::conformer_row` | inclusive | 2 | 96.601 | 1000 | 29.417× | 1.934× |
| `parse_ensembles_1000` / `commands::parse::conformer_row` | inclusive | 4 | 96.497 | 1000 | 28.546× | 1.932× |
| `parse_ensembles_1000` / `commands::parse::conformer_row` | inclusive | 6 | 96.478 | 1000 | 28.393× | 1.932× |
| `db_build_ensembles` / `commands::db::featurize_wall` | inclusive | 1 | 86.323 | 1 | 7.312× | 1.759× |
| `db_build_ensembles` / `commands::db::featurize_wall` | inclusive | 2 | 71.094 | 1 | 3.459× | 1.552× |
| `db_build_ensembles` / `commands::db::featurize_wall` | inclusive | 4 | 45.391 | 1 | 1.831× | 1.294× |
| `db_build_ensembles` / `commands::db::featurize_wall` | inclusive | 6 | 34.040 | 1 | 1.516× | 1.205× |
| `predict_1000000` / `model::RegressXPredictor::predict_batch` | inclusive | 1 | 67.549 | 1 | 3.082× | 1.510× |
| `predict_1000000` / `model::RegressXPredictor::predict_batch` | inclusive | 2 | 51.842 | 1 | 2.076× | 1.350× |
| `predict_1000000` / `model::RegressXPredictor::predict_batch` | inclusive | 4 | 35.626 | 1 | 1.553× | 1.217× |
| `predict_1000000` / `model::RegressXPredictor::predict_batch` | inclusive | 6 | 33.393 | 1 | 1.501× | 1.200× |

## Allocation counts and diagnostic overhead

Counts cover successful requested System allocations, including new reallocation sizes and process startup/session metadata, before report construction. They exclude mappings, stacks, allocator metadata and timing-report construction. Requested traffic and heap high-water are not RSS or allocation-time measurements.

| Workload | Threads | Alloc + realloc calls | Requested MiB | Peak live MiB | Diagnostic/native wall |
| --- | ---: | ---: | ---: | ---: | ---: |
| `conformers_56` | 1 | 18398 | 22.268 | 0.413 | 1.072× |
| `conformers_56` | 2 | 18398 | 22.268 | 0.413 | 1.082× |
| `conformers_56` | 4 | 18398 | 22.268 | 0.413 | 1.047× |
| `conformers_56` | 6 | 18398 | 22.268 | 0.413 | 1.060× |
| `db_build_ensembles` | 1 | 18920 | 22.247 | 0.434 | 1.107× |
| `db_build_ensembles` | 2 | 18931 | 22.254 | 0.820 | 1.234× |
| `db_build_ensembles` | 4 | 18953 | 22.267 | 1.587 | 1.749× |
| `db_build_ensembles` | 6 | 18975 | 22.281 | 2.360 | 2.378× |
| `descriptors_10000` | 1 | 3.16137e+06 | 3924.938 | 13.591 | 1.010× |
| `descriptors_10000` | 2 | 3.16137e+06 | 3924.938 | 13.591 | 0.984× |
| `descriptors_10000` | 4 | 3.16137e+06 | 3924.938 | 13.591 | 1.005× |
| `descriptors_10000` | 6 | 3.16137e+06 | 3924.938 | 13.591 | 1.007× |
| `ensemble_sdf` | 1 | 42501 | 38.883 | 1.248 | 1.054× |
| `ensemble_sdf` | 2 | 42501 | 38.883 | 1.248 | 1.060× |
| `ensemble_sdf` | 4 | 42501 | 38.883 | 1.248 | 1.048× |
| `ensemble_sdf` | 6 | 42501 | 38.883 | 1.248 | 1.039× |
| `parse_ensembles_1000` | 1 | 889109 | 55.205 | 0.158 | 1.036× |
| `parse_ensembles_1000` | 2 | 889109 | 55.205 | 0.158 | 1.037× |
| `parse_ensembles_1000` | 4 | 889109 | 55.205 | 0.158 | 1.021× |
| `parse_ensembles_1000` | 6 | 889109 | 55.205 | 0.158 | 1.043× |
| `predict_1000000` | 1 | 655 | 4.082 | 3.830 | 1.047× |
| `predict_1000000` | 2 | 666 | 4.089 | 3.837 | 1.086× |
| `predict_1000000` | 4 | 688 | 4.102 | 3.851 | 1.137× |
| `predict_1000000` | 6 | 710 | 4.116 | 3.864 | 1.070× |
| `screen_1000` | 1 | 934684 | 86.845 | 5.058 | 1.044× |
| `screen_1000` | 2 | 934684 | 86.845 | 5.058 | 1.042× |
| `screen_1000` | 4 | 934684 | 86.845 | 5.058 | 1.033× |
| `screen_1000` | 6 | 934684 | 86.845 | 5.058 | 1.045× |
| `search_database` | 1 | 23878 | 8.116 | 3.107 | 1.349× |
| `search_database` | 2 | 23878 | 8.116 | 3.107 | 1.412× |
| `search_database` | 4 | 23878 | 8.116 | 3.107 | 1.361× |
| `search_database` | 6 | 23878 | 8.116 | 3.107 | 1.416× |
| `single_large` | 1 | 1133 | 0.653 | 0.381 | 2.111× |
| `single_large` | 2 | 1133 | 0.653 | 0.381 | 1.918× |
| `single_large` | 4 | 1133 | 0.653 | 0.381 | 1.912× |
| `single_large` | 6 | 1133 | 0.653 | 0.381 | 2.020× |
| `single_small` | 1 | 935 | 0.643 | 0.379 | 1.895× |
| `single_small` | 2 | 935 | 0.643 | 0.379 | 2.027× |
| `single_small` | 4 | 935 | 0.643 | 0.379 | 2.053× |
| `single_small` | 6 | 935 | 0.643 | 0.379 | 1.749× |

## Scope and remaining profiling

The new SDF ensemble preserves original bond records. The search database was freshly computed from 1,543 predeclared whole ensembles (31,618 conformers), with all 23 exclusions recorded before native execution and zero skipped selected geometries. Screening uses explicit supplied-weight means at response temperature 353.15 K, preserving historical 298.15 K populations. Repetition to larger throughput workloads does not add new chemistry. These changes prevent direct performance comparison with the earlier invalid workloads.

Fresh hardware-counter capability checks and separately scheduled Cachegrind/DHAT runs remain pending. No historical counter permission, instruction count or allocation site percentage is carried forward. Modeled events will be labeled separately from native hardware measurements and will not provide wall-time speedups.

Scientific admission is scoped to retained equations, evidence and explicit unavailable cases. It does not establish empirical calibration, prospective validation, reaction-temperature population provenance, or resolution of the published-source conflicts. Candidate changes still require exact scientific equivalence, separate independent references and paired native measurements against this same preserved baseline.
