# Accepted C1b diagnostic reprofile

All ten unchanged workloads were profiled at 1, 2, 4 and 6 threads, with one warmup and five diagnostic samples per configuration. The overhead denominator is the median of the eight actual candidate samples retained in the accepted paired native benchmark. No extra native timing samples were generated or discarded.

Worker scopes measure elapsed time, not CPU time. The main-thread parallel scope is the submission/join envelope. It is neither an isolated scheduler-wait measurement nor additive to worker scopes. Per-worker CPU fractions and isolated scheduler wait are unavailable.

| Threads | Workload | Native A/B ms | Diagnostic/native | Native CPU % of one core | Native peak RSS bytes | Diagnostic requested bytes |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `conformers_56` | 29.198 | 1.083 | 99.2 | 5,384,192 | 23,349,310 |
| 1 | `db_build_ensembles` | 29.944 | 1.088 | 98.9 | 5,859,328 | 23,327,793 |
| 1 | `descriptors_10000` | 4974.639 | 1.021 | 99.7 | 22,040,576 | 4,115,596,207 |
| 1 | `ensemble_sdf` | 49.471 | 1.040 | 99.3 | 6,152,192 | 40,771,655 |
| 1 | `parse_ensembles_1000` | 185.877 | 1.082 | 99.8 | 5,150,720 | 57,886,873 |
| 1 | `predict_1000000` | 20.914 | 1.109 | 99.3 | 72,962,048 | 4,280,533 |
| 1 | `screen_1000` | 258.709 | 1.012 | 97.5 | 10,811,392 | 91,063,713 |
| 1 | `search_database` | 14.764 | 0.662 | 63.1 | 8,134,656 | 8,509,766 |
| 1 | `single_large` | 2.120 | 1.826 | 88.5 | 5,269,504 | 684,557 |
| 1 | `single_small` | 2.214 | 1.786 | 88.0 | 5,230,592 | 673,952 |
| 2 | `conformers_56` | 16.321 | 1.255 | 184.2 | 6,289,408 | 23,380,455 |
| 2 | `db_build_ensembles` | 16.782 | 1.356 | 181.1 | 6,283,264 | 23,334,857 |
| 2 | `descriptors_10000` | 2570.516 | 1.029 | 196.6 | 26,015,744 | 4,116,979,736 |
| 2 | `ensemble_sdf` | 46.756 | 1.053 | 99.5 | 6,121,472 | 40,771,655 |
| 2 | `parse_ensembles_1000` | 185.058 | 1.038 | 99.9 | 5,136,384 | 57,886,873 |
| 2 | `predict_1000000` | 13.760 | 1.099 | 155.3 | 72,941,568 | 4,287,597 |
| 2 | `screen_1000` | 261.700 | 1.006 | 95.5 | 10,643,456 | 91,063,713 |
| 2 | `search_database` | 8.924 | 1.089 | 84.5 | 8,130,560 | 8,509,766 |
| 2 | `single_large` | 2.120 | 1.905 | 90.5 | 5,279,744 | 684,557 |
| 2 | `single_small` | 2.116 | 1.886 | 90.0 | 5,275,648 | 673,952 |
| 4 | `conformers_56` | 9.699 | 1.669 | 329.1 | 6,950,912 | 23,398,727 |
| 4 | `db_build_ensembles` | 10.921 | 1.576 | 291.5 | 6,893,568 | 23,348,985 |
| 4 | `descriptors_10000` | 1322.822 | 1.018 | 387.9 | 26,771,456 | 4,116,998,008 |
| 4 | `ensemble_sdf` | 46.858 | 1.069 | 99.5 | 6,100,992 | 40,771,655 |
| 4 | `parse_ensembles_1000` | 185.408 | 1.045 | 99.9 | 5,087,232 | 57,886,873 |
| 4 | `predict_1000000` | 10.363 | 1.089 | 218.1 | 72,781,824 | 4,301,725 |
| 4 | `screen_1000` | 253.440 | 1.058 | 99.5 | 10,647,552 | 91,063,713 |
| 4 | `search_database` | 7.328 | 1.356 | 96.5 | 8,130,560 | 8,509,766 |
| 4 | `single_large` | 2.272 | 1.825 | 90.4 | 5,242,880 | 684,557 |
| 4 | `single_small` | 2.262 | 1.768 | 90.4 | 5,275,648 | 673,952 |
| 6 | `conformers_56` | 7.354 | 2.385 | 457.5 | 7,677,952 | 23,417,191 |
| 6 | `db_build_ensembles` | 7.756 | 2.355 | 434.8 | 7,591,936 | 23,363,305 |
| 6 | `descriptors_10000` | 914.970 | 1.041 | 576.4 | 27,613,184 | 4,117,016,472 |
| 6 | `ensemble_sdf` | 46.946 | 1.077 | 99.5 | 6,029,312 | 40,771,655 |
| 6 | `parse_ensembles_1000` | 186.757 | 1.045 | 99.8 | 5,076,992 | 57,886,873 |
| 6 | `predict_1000000` | 9.121 | 1.131 | 255.9 | 72,769,536 | 4,316,045 |
| 6 | `screen_1000` | 260.997 | 1.021 | 95.9 | 10,586,112 | 91,063,713 |
| 6 | `search_database` | 12.153 | 0.844 | 72.0 | 8,079,360 | 8,509,766 |
| 6 | `single_large` | 2.334 | 1.879 | 90.9 | 5,281,792 | 684,557 |
| 6 | `single_small` | 2.367 | 1.846 | 89.9 | 5,283,840 | 673,952 |

The native denominator and diagnostic numerator come from separate series. Ratios below one reflect series/build variation; they do not establish a negative instrumentation cost. Full distributions, main-thread rankings, allocations and per-worker records are retained in [current_c1b_reprofile.json](current_c1b_reprofile.json).

## Screen interval calculation

| Threads | Calls per screen | Inclusive diagnostic wall share | Infinite local ceiling |
| ---: | ---: | ---: | ---: |
| 1 | 1000 | 20.246% | 1.254 |
| 2 | 1000 | 20.188% | 1.253 |
| 4 | 1000 | 20.018% | 1.250 |
| 6 | 1000 | 20.037% | 1.251 |

## Parallel descriptor work

| Threads | Main join envelope, wall share | Occupancy share of summed worker elapsed | File calls per active worker, min-max | Median max/mean file count |
| ---: | ---: | ---: | --- | ---: |
| 1 | 0.000% | unavailable (serial) | none | — |
| 2 | 98.147% | 83.414% | 5000-5000 | 1.000 |
| 4 | 96.845% | 83.208% | 2426-2548 | 1.006 |
| 6 | 94.888% | 82.867% | 1597-1718 | 1.013 |

Complete scope, per-worker, allocation and native resource distributions are retained in the reprofile receipt. Worker elapsed shares do not yield a valid process-wall Amdahl ceiling. Screen's serial inclusive ceiling is a diagnostic bound, not a measured native speedup.

Worker file-count ranges cover all five diagnostic samples. The max/mean column is the median of each sample's ratio; no scheduling outlier is removed. Serial batch occupancy accounts for 82.533% of diagnostic wall time. Screen quantile call counts and every sample's elapsed share are retained without rounding in the accompanying JSON. The serial diagnostic ceiling is a bound for that observed component, not a native speedup prediction.

## Provenance and validation

This report profiles the accepted frozen C1b source and executable, including the CLI refactor, on the same corrected workloads. It establishes the immediate baseline for a later quantile candidate; it does not isolate the CLI refactor's effect against C1. A later candidate still requires exact scientific gates, its full paired comparison against the original corrected baseline, and a separately labelled incremental comparison against its accepted predecessor.

All 240 diagnostic launches completed with exact same-thread native scientific fingerprints including stderr. Stage reports recorded no dropped scopes. The native denominators retain all 320 original measured C1b candidate samples and their 40 separate warmups; this reprofile generated no native timing samples.

The accepted native executable SHA-256 is `e7e7f00a1b4802e219144535d5c6567eb5e97bd269d1073e931b880288c5023d`; the diagnostic executable is `12d94fac2c2fda5aa67ba28966f986dbf3052c4090ca2483fb63ea17c3aa7d16`. The reprofile manifest SHA-256 is `602fb4e9c954ff0d5331c1aa282e5f0903433bcdcd668e16897965c330f7f0b5` at `.stericx/profiling/scientifically_validated_optimization/accepted_c1b_reprofile_v1/manifest.json`. Its inventory binds 1,309 retained files.

The unchanged controller is `candidate_reprofile_controller_v1/run.py`, with the C1b acceptance, build and paired receipts supplied explicitly. Its historical internal report filename remains `C1_REPROFILE.md` inside the new C1b output directory; this publication labels the actual C1b inputs. The outer launch receipt is `accepted_c1b_reprofile_launch_v1/result.json`. Publication is additive and verifies that the original C1 report and compact data remain unchanged.
