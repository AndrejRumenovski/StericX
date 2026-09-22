# Accepted C2 diagnostic reprofile

All ten unchanged workloads were profiled at 1, 2, 4 and 6 threads, with one warmup and five diagnostic samples per configuration. The overhead denominator is the median of the eight actual candidate samples retained in the accepted paired native benchmark. No extra native timing samples were generated or discarded.

Worker scopes measure elapsed time, not CPU time. The main-thread parallel scope is the submission/join envelope. It is neither an isolated scheduler-wait measurement nor additive to worker scopes. Per-worker CPU fractions and isolated scheduler wait are unavailable.

| Threads | Workload | Native A/B ms | Diagnostic/native | Native CPU % of one core | Native peak RSS bytes | Diagnostic requested bytes |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `conformers_56` | 28.842 | 1.098 | 99.1 | 5,423,104 | 23,349,308 |
| 1 | `db_build_ensembles` | 29.818 | 1.090 | 98.8 | 5,865,472 | 23,327,781 |
| 1 | `descriptors_10000` | 4909.577 | 1.005 | 99.6 | 22,022,144 | 4,115,596,205 |
| 1 | `ensemble_sdf` | 45.982 | 1.048 | 99.4 | 6,187,008 | 40,771,653 |
| 1 | `parse_ensembles_1000` | 187.315 | 1.016 | 99.8 | 5,234,688 | 57,886,868 |
| 1 | `predict_1000000` | 20.960 | 1.043 | 99.2 | 72,927,232 | 4,280,531 |
| 1 | `screen_1000` | 204.535 | 1.022 | 98.6 | 10,747,904 | 91,063,711 |
| 1 | `search_database` | 15.003 | 0.693 | 50.8 | 8,116,224 | 8,509,764 |
| 1 | `single_large` | 2.110 | 1.818 | 89.4 | 5,265,408 | 684,555 |
| 1 | `single_small` | 2.029 | 1.918 | 89.9 | 5,361,664 | 673,950 |
| 2 | `conformers_56` | 16.005 | 1.250 | 183.9 | 6,387,712 | 23,380,453 |
| 2 | `db_build_ensembles` | 16.608 | 1.269 | 180.8 | 6,266,880 | 23,334,845 |
| 2 | `descriptors_10000` | 2509.911 | 1.021 | 197.8 | 26,150,912 | 4,116,979,734 |
| 2 | `ensemble_sdf` | 46.251 | 1.062 | 99.5 | 6,035,456 | 40,771,653 |
| 2 | `parse_ensembles_1000` | 188.243 | 1.037 | 99.8 | 5,122,048 | 57,886,868 |
| 2 | `predict_1000000` | 13.770 | 1.072 | 155.3 | 72,906,752 | 4,287,595 |
| 2 | `screen_1000` | 207.037 | 1.028 | 97.1 | 10,641,408 | 91,063,711 |
| 2 | `search_database` | 9.160 | 1.049 | 79.9 | 8,157,184 | 8,509,764 |
| 2 | `single_large` | 2.189 | 1.754 | 90.3 | 5,300,224 | 684,555 |
| 2 | `single_small` | 2.100 | 1.933 | 90.4 | 5,249,024 | 673,950 |
| 4 | `conformers_56` | 9.540 | 1.679 | 328.4 | 7,038,976 | 23,398,725 |
| 4 | `db_build_ensembles` | 10.206 | 1.746 | 311.1 | 7,016,448 | 23,348,973 |
| 4 | `descriptors_10000` | 1301.219 | 1.042 | 389.7 | 26,656,768 | 4,116,998,006 |
| 4 | `ensemble_sdf` | 46.356 | 1.061 | 99.5 | 6,078,464 | 40,771,653 |
| 4 | `parse_ensembles_1000` | 187.939 | 1.058 | 99.8 | 5,093,376 | 57,886,868 |
| 4 | `predict_1000000` | 10.265 | 1.188 | 216.9 | 72,710,144 | 4,301,723 |
| 4 | `screen_1000` | 207.913 | 1.077 | 97.8 | 10,676,224 | 91,063,711 |
| 4 | `search_database` | 12.900 | 0.763 | 63.0 | 8,155,136 | 8,509,764 |
| 4 | `single_large` | 2.217 | 1.809 | 90.7 | 5,201,920 | 684,555 |
| 4 | `single_small` | 2.199 | 1.838 | 90.9 | 5,287,936 | 673,950 |
| 6 | `conformers_56` | 7.132 | 2.415 | 462.0 | 7,628,800 | 23,417,189 |
| 6 | `db_build_ensembles` | 7.522 | 2.327 | 441.7 | 7,651,328 | 23,363,293 |
| 6 | `descriptors_10000` | 903.662 | 1.046 | 573.8 | 27,676,672 | 4,117,016,470 |
| 6 | `ensemble_sdf` | 46.644 | 1.069 | 99.4 | 6,187,008 | 40,771,653 |
| 6 | `parse_ensembles_1000` | 187.078 | 1.050 | 99.8 | 5,021,696 | 57,886,868 |
| 6 | `predict_1000000` | 8.860 | 1.128 | 264.1 | 72,820,736 | 4,316,043 |
| 6 | `screen_1000` | 203.228 | 1.058 | 99.9 | 10,450,944 | 91,063,711 |
| 6 | `search_database` | 10.083 | 0.963 | 75.4 | 8,142,848 | 8,509,764 |
| 6 | `single_large` | 2.279 | 1.790 | 89.5 | 5,246,976 | 684,555 |
| 6 | `single_small` | 2.280 | 1.816 | 91.4 | 5,271,552 | 673,950 |

The native denominator and diagnostic numerator come from separate series. Ratios below one reflect series/build variation; they do not establish a negative instrumentation cost. Full distributions, main-thread rankings, allocations and per-worker records are retained in [current_c2_reprofile.json](current_c2_reprofile.json).

## Screen interval calculation

| Threads | Calls per screen | Inclusive diagnostic wall share | Infinite local ceiling |
| ---: | ---: | ---: | ---: |
| 1 | 1 | 0.030% | 1.000 |
| 2 | 1 | 0.029% | 1.000 |
| 4 | 1 | 0.029% | 1.000 |
| 6 | 1 | 0.030% | 1.000 |

## Parallel descriptor work

| Threads | Main join envelope, wall share | Occupancy share of summed worker elapsed | File calls per active worker, min-max | Median max/mean file count |
| ---: | ---: | ---: | --- | ---: |
| 1 | 0.000% | unavailable (serial) | none | — |
| 2 | 98.077% | 83.302% | 5000-5000 | 1.000 |
| 4 | 96.418% | 83.082% | 2475-2528 | 1.000 |
| 6 | 94.447% | 82.686% | 1562-1810 | 1.011 |

Complete scope, per-worker, allocation and native resource distributions are retained in the reprofile receipt. Worker elapsed shares do not yield a valid process-wall Amdahl ceiling. Screen's serial inclusive ceiling is a diagnostic bound, not a measured native speedup.

## Ranked serial components

These are the five largest exclusive main-thread scopes in each selected one-thread workload. Exclusive shares avoid counting child scopes twice. The ceiling assumes the selected scope costs zero while all other work stays fixed; it is not a measured native speedup.

| Workload | Function | Calls | Diagnostic wall share | Infinite local ceiling |
| --- | --- | ---: | ---: | ---: |
| `descriptors_10000` | `buried_volume::occupied_volumes` | 30,000 | 82.464% | 5.702 |
| `descriptors_10000` | `buried_volume::integration_grid` | 10,000 | 8.222% | 1.090 |
| `descriptors_10000` | `sterimol::params_from_projection` | 10,000 | 2.959% | 1.030 |
| `descriptors_10000` | `parse_xyz` | 10,000 | 2.500% | 1.026 |
| `descriptors_10000` | `parse_coordinate_file` | 10,000 | 1.445% | 1.015 |
| `screen_1000` | `sterimol::params_from_projection` | 5,085 | 34.177% | 1.519 |
| `screen_1000` | `parse_xyz` | 5,085 | 24.302% | 1.321 |
| `screen_1000` | `parse_coordinate_file` | 5,085 | 13.804% | 1.160 |
| `screen_1000` | `commands::screen::bootstrap_mean_response_interval` | 1,000 | 6.110% | 1.065 |
| `screen_1000` | `reaction::resolve_xyz_path` | 5,085 | 5.992% | 1.064 |

## Screening calls and allocation effects

Every one of the 20 measured C2 screens records 1,000 prediction-interval calls, one multiplier initialization and one Student-t quantile call. All 20 retained C1b screens record 1,000 interval and 1,000 quantile calls. The separate cache-boundary source review and scientific gates establish lazy initialization after successful leverage evaluation; aggregate scope counts alone do not establish execution order.

| Requested-heap metric | C1b, every measured screen | C2, every measured screen |
| --- | ---: | ---: |
| `allocation_calls` | 931,731 | 931,731 |
| `reallocation_calls` | 2,953 | 2,953 |
| `deallocation_calls` | 931,682 | 931,682 |
| `failed_allocation_calls` | 0 | 0 |
| `requested_bytes` | 91,063,713 | 91,063,711 |
| `live_requested_bytes` | 2,082 | 2,081 |
| `peak_live_requested_bytes` | 5,304,019 | 5,304,018 |

Allocation-call counts remain unchanged. Requested bytes include executable and profile-path metadata, whose path lengths differ between these captures; the small byte deltas do not establish memory savings from the cache. Native peak RSS for all 40 configurations is reported above and its complete sample distributions remain in the JSON.

Worker file-count ranges cover all five diagnostic samples. The max/mean column is the median of each sample's ratio; no scheduling outlier is removed. Serial batch occupancy accounts for 82.464% of diagnostic wall time. Screen quantile call counts and every sample's elapsed share are retained without rounding in the accompanying JSON. The serial diagnostic ceiling is a bound for that observed component, not a native speedup prediction.

## Provenance and validation

This report profiles the accepted frozen C2 source and executable, including lazy interval-quantile reuse, on the same corrected workloads. Its native denominators are the C2 samples from the full original-baseline campaign, while the separate incremental campaign measures C2 against C1b. This reprofile establishes the immediate baseline for the proposed occupancy candidate. A later candidate still requires exact scientific gates, its full paired comparison against the original corrected baseline, and a separately labelled incremental comparison against its accepted predecessor.

All 240 diagnostic launches completed with exact same-thread native scientific fingerprints including stderr. Stage reports recorded no dropped scopes. The native denominators retain all 320 original measured C2 candidate samples and their 40 separate warmups; this reprofile generated no native timing samples.

The accepted native executable SHA-256 is `b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`; the diagnostic executable is `693c875b9a877b64aade08bcd71b4d3c232ec589231217ce53bbe80cb9a743a7`. The reprofile manifest SHA-256 is `52e7d287332afdb1edf13348cb73e7bccc243e2c0259d187c5705dc675426187` at `.stericx/profiling/scientifically_validated_optimization/accepted_c2_reprofile_v1/manifest.json`. Its inventory binds 1,309 retained files.

The unchanged controller is `candidate_reprofile_controller_v1/run.py`, with the C2 acceptance, build and paired receipts supplied explicitly. Its historical internal report filename remains `C1_REPROFILE.md` inside the new C2 output directory; this publication labels the actual C2 inputs. The outer launch receipt is `accepted_c2_reprofile_launch_v1/result.json`. Publication is additive and verifies that the original C1 and C1b reports and compact data remain unchanged.
