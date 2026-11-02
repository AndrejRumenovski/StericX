# Native allocation diagnostics

Instrumented runs are excluded from timing claims. Counts below are identical for both binaries. Requested-byte and peak differences are final minus baseline.

| Threads | Workload | Allocation calls | Reallocation calls | Requested-byte difference | Peak live requested-byte difference |
|---:|---|---:|---:|---:|---:|
| 1 | conformers_56 | 17743 | 656 | -2 | 1 |
| 1 | db_build_ensembles | 18143 | 777 | 8 | 2 |
| 1 | descriptors_10000 | 3056774 | 104600 | -2 | 1 |
| 1 | ensemble_sdf | 24990 | 17512 | -2 | 1 |
| 1 | parse_ensembles_1000 | 886212 | 2897 | 1 | 1 |
| 1 | single_large | 1073 | 61 | -2 | 1 |
| 1 | single_small | 878 | 58 | -2 | 1 |
| 6 | conformers_56 | 17827 | 656 | -2 | -2303 |
| 6 | db_build_ensembles | 18198 | 777 | 8 | -1242 |
| 6 | descriptors_10000 | 3056858 | 104600 | -2 | 1 |
| 6 | ensemble_sdf | 24990 | 17512 | -2 | 1 |
| 6 | parse_ensembles_1000 | 886212 | 2897 | 1 | 1 |
| 6 | single_large | 1073 | 61 | -2 | 1 |
| 6 | single_small | 878 | 58 | -2 | 1 |
