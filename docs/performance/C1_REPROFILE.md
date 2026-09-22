# Accepted C1 diagnostic reprofile

All ten unchanged workloads were profiled at 1, 2, 4 and 6 threads, with one warmup and five diagnostic samples per configuration. The overhead denominator is the median of the eight actual candidate samples retained in the accepted paired native benchmark. No extra native timing samples were generated or discarded.

Worker scopes measure elapsed time, not CPU time. The main-thread parallel scope is the submission/join envelope. It is neither an isolated scheduler-wait measurement nor additive to worker scopes. Per-worker CPU fractions and isolated scheduler wait are unavailable.

| Threads | Workload | Native A/B ms | Diagnostic/native | Native CPU % of one core | Native peak RSS bytes | Diagnostic requested bytes |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `conformers_56` | 30.176 | 1.043 | 99.3 | 5,429,248 | 23,349,308 |
| 1 | `db_build_ensembles` | 30.616 | 1.070 | 99.3 | 5,836,800 | 23,327,781 |
| 1 | `descriptors_10000` | 5043.202 | 0.984 | 99.8 | 22,009,856 | 4,115,596,205 |
| 1 | `ensemble_sdf` | 47.438 | 1.057 | 99.4 | 6,123,520 | 40,771,653 |
| 1 | `parse_ensembles_1000` | 187.995 | 1.013 | 99.9 | 5,199,872 | 57,886,868 |
| 1 | `predict_1000000` | 22.045 | 0.992 | 99.3 | 73,000,960 | 4,280,531 |
| 1 | `screen_1000` | 262.909 | 0.993 | 99.8 | 10,676,224 | 91,063,711 |
| 1 | `search_database` | 9.137 | 1.051 | 83.2 | 8,192,000 | 8,509,764 |
| 1 | `single_large` | 2.169 | 1.859 | 89.7 | 5,300,224 | 684,555 |
| 1 | `single_small` | 2.157 | 1.872 | 90.2 | 5,298,176 | 673,950 |
| 2 | `conformers_56` | 16.371 | 1.223 | 183.5 | 6,397,952 | 23,380,453 |
| 2 | `db_build_ensembles` | 17.014 | 1.252 | 180.7 | 6,283,264 | 23,334,845 |
| 2 | `descriptors_10000` | 2602.054 | 0.984 | 196.0 | 26,009,600 | 4,116,979,734 |
| 2 | `ensemble_sdf` | 47.698 | 1.025 | 99.5 | 6,086,656 | 40,771,653 |
| 2 | `parse_ensembles_1000` | 189.298 | 1.008 | 99.8 | 5,036,032 | 57,886,868 |
| 2 | `predict_1000000` | 14.758 | 1.000 | 153.3 | 72,941,568 | 4,287,595 |
| 2 | `screen_1000` | 264.510 | 1.065 | 96.9 | 10,598,400 | 91,063,711 |
| 2 | `search_database` | 12.480 | 0.741 | 70.7 | 8,155,136 | 8,509,764 |
| 2 | `single_large` | 2.212 | 1.818 | 91.1 | 5,351,424 | 684,555 |
| 2 | `single_small` | 2.166 | 1.845 | 90.0 | 5,279,744 | 673,950 |
| 4 | `conformers_56` | 9.819 | 1.638 | 327.8 | 7,081,984 | 23,398,725 |
| 4 | `db_build_ensembles` | 10.835 | 1.583 | 303.5 | 6,879,232 | 23,348,973 |
| 4 | `descriptors_10000` | 1342.070 | 0.997 | 386.9 | 26,720,256 | 4,116,998,006 |
| 4 | `ensemble_sdf` | 47.733 | 1.037 | 99.4 | 6,078,464 | 40,771,653 |
| 4 | `parse_ensembles_1000` | 185.511 | 1.033 | 99.8 | 5,042,176 | 57,886,868 |
| 4 | `predict_1000000` | 11.139 | 0.993 | 214.4 | 72,722,432 | 4,301,723 |
| 4 | `screen_1000` | 256.434 | 1.028 | 98.3 | 10,680,320 | 91,063,711 |
| 4 | `search_database` | 9.615 | 0.998 | 83.3 | 8,153,088 | 8,509,764 |
| 4 | `single_large` | 2.325 | 1.905 | 91.0 | 5,367,808 | 684,555 |
| 4 | `single_small` | 2.322 | 1.701 | 90.1 | 5,355,520 | 673,950 |
| 6 | `conformers_56` | 7.346 | 2.295 | 455.1 | 7,735,296 | 23,417,189 |
| 6 | `db_build_ensembles` | 7.912 | 2.208 | 433.2 | 7,602,176 | 23,363,293 |
| 6 | `descriptors_10000` | 940.472 | 0.998 | 568.4 | 27,545,600 | 4,117,016,470 |
| 6 | `ensemble_sdf` | 47.586 | 1.051 | 99.4 | 6,148,096 | 40,771,653 |
| 6 | `parse_ensembles_1000` | 186.021 | 1.050 | 99.9 | 5,193,728 | 57,886,868 |
| 6 | `predict_1000000` | 9.845 | 1.063 | 252.4 | 72,722,432 | 4,316,043 |
| 6 | `screen_1000` | 256.165 | 1.037 | 98.4 | 10,678,272 | 91,063,711 |
| 6 | `search_database` | 15.201 | 0.609 | 50.0 | 8,134,656 | 8,509,764 |
| 6 | `single_large` | 2.355 | 1.792 | 90.3 | 5,306,368 | 684,555 |
| 6 | `single_small` | 2.205 | 1.926 | 90.6 | 5,300,224 | 673,950 |

The native timing denominator and diagnostic numerator come from separate series. Ratios below one reflect series/build variation; they do not establish a negative instrumentation cost. The full distributions, main-thread rankings, allocations and per-worker records are in [current_c1_reprofile.json](current_c1_reprofile.json).

## Screen interval calculation

| Threads | Calls per screen | Inclusive diagnostic wall share | Infinite local ceiling |
| ---: | ---: | ---: | ---: |
| 1 | 1000 | 20.323% | 1.255 |
| 2 | 1000 | 18.869% | 1.233 |
| 4 | 1000 | 20.231% | 1.254 |
| 6 | 1000 | 20.019% | 1.250 |

## Parallel descriptor work

| Threads | Main join envelope, wall share | Occupancy share of summed worker elapsed | File calls per active worker, min-max | Median max/mean file count |
| ---: | ---: | ---: | --- | ---: |
| 1 | 0.000% | unavailable (serial) | none | — |
| 2 | 98.592% | 82.733% | 5000-5000 | 1.000 |
| 4 | 96.825% | 82.451% | 2460-2530 | 1.004 |
| 6 | 95.780% | 82.094% | 1615-1875 | 1.010 |

Complete scope, per-worker, allocation and native resource distributions are retained in the reprofile receipt. Worker elapsed shares do not yield a valid process-wall Amdahl ceiling. Screen's serial inclusive ceiling is a diagnostic bound, not a measured native speedup.

Worker file-count ranges cover all five diagnostic samples. The max/mean column is the median of each sample's ratio; no scheduling outlier is removed. For the serial batch, occupancy itself accounts for 81.967% of diagnostic wall time. Screening observes one profiled thread at every Rayon setting. These results motivate a bounded quantile reuse candidate; exact scientific gates and new paired native measurements remain required.

## Provenance and validation

This report profiles the accepted frozen C1 executable and source, using the corrected workload manifest. Admission and postflight passed unchanged despite subsequent working-tree edits. All 240 diagnostic launches completed with exact same-thread native scientific fingerprints including stderr; stage reports recorded no dropped scopes. The 40 native denominators retain all 320 original measured candidate samples and the 40 separate warmups.

The accepted native executable SHA-256 is `2f65e0aba46a502bb653c61ac3090fa0acbb0244dedcf4aba18ce06e458be358`; the diagnostic executable is `dfb6699f5eb590a02c0b9c696fcc24229a091ac70a99ff1be99cf470bb83c094`. The full reprofile manifest is `65b73ff4b7321ffd0becd589098b1049aaf092f783fadf40de0d021aae51645f` at `.stericx/profiling/scientifically_validated_optimization/accepted_c1_reprofile_v1/manifest.json`. Its inventory binds 1,309 retained files, including the native adapters, raw diagnostic outputs/resources, scopes, helper snapshots, commands and reports.

The unchanged controller is `candidate_reprofile_controller_v1/run.py`; its eight focused tests and Ruff passed before execution. The launch receipt is `accepted_c1_reprofile_launch_v1/result.json`. Publication preserves the earlier pending plan and binds the generated report; it does not alter the completed profiling evidence.
