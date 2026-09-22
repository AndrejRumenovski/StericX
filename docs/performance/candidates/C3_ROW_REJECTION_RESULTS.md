# C3 row-endpoint rejection: rejected native result

C3 passed the scientific gates but failed its predeclared performance gate. The exact accepted C2 buried-volume implementation has been restored. Against contemporaneously measured C2, the 10,000-file descriptor batch was **15.57% slower at one thread and 15.69% slower at six threads** by median paired wall cost. Both launch orders agree, and all 16 target CPU pairs are slower. No timing supplement or threshold change was needed to reject it.

The [pre-edit addendum](C3_IMPLEMENTATION_ADDENDUM.md) required at least 10% lower target batch wall time; the decision applies that threshold at both one and six threads, with support in both execution orders. The original 15–25% serial improvement hypothesis was not realized. The [operation-count evidence](../C3_OPERATION_COUNTS.md) remains valid: fewer selected source predicates did not produce a faster native executable.

## Verified paired measurements

Two quiet native campaigns retained every sample: the mandatory original-corrected-baseline matrix contains 40 configurations, 720 launches and 320 measured pairs; the fixed accepted-C2 comparison contains 14 configurations, 252 launches and 112 measured pairs. Each configuration has one warmup per binary and eight measured pairs, alternating AB/BA. The independent reviewers verified all **972 launches and 432 measured pairs**, complete scientific outputs and stderr, raw child CPU/RSS/fault records, sample order, frozen inputs/binaries/helpers, scientific gates and source/postflight hashes before C2 restoration. No samples were excluded.

Wall speedup below 1 means C3 is slower. Paired cost is the median of `100 × (C3/C2 − 1)` for each pair; with an even sample count, it need not equal the inverse of the median speedup. MAD and the full range describe the paired speedup values. The compact [review receipt](../current_c3_rejection_review.json) retains every value, CPU distribution, order effect, RSS/fault distribution and all 54 configuration rows.

| Target | C2 wall median ms | C3 wall median ms | Paired wall cost | Speedup median ± MAD | Full speedup range | AB / BA speedup | Slower wall / CPU pairs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| t1/descriptors_10000 | 4876.903 | 5657.123 | +15.570% | 0.865279 ± 0.002957 | 0.855242–0.868617 | 0.862090 / 0.866252 | 8/8 / 8/8 |
| t6/descriptors_10000 | 921.886 | 1065.886 | +15.685% | 0.864643 ± 0.080518 | 0.615657–1.251668 | 0.845478 / 0.890856 | 6/8 / 8/8 |

The six-thread target has a broad wall range and two faster pairs, all retained. Its paired median and both order medians still regress; the serial range is entirely slower. Target paired CPU cost rises 15.55% and 15.29%, respectively. The result is a clear failed benefit gate, not an inference from the original-baseline cumulative gains.

## All incremental workload costs

Positive wall/CPU cost means C3 is slower; positive RSS delta means a higher median peak RSS. Each row contains eight measured pairs. RSS changes here are small (at most 0.100 MiB in magnitude), which does not establish zero memory cost.

| Configuration | Paired wall cost | Paired CPU cost | Median RSS delta MiB | Slower wall / CPU pairs |
| --- | ---: | ---: | ---: | ---: |
| t1/descriptors_10000 | +15.570% | +15.552% | +0.012 | 8/8 / 8/8 |
| t6/descriptors_10000 | +15.685% | +15.294% | +0.010 | 6/8 / 8/8 |
| t1/conformers_56 | +13.949% | +14.199% | +0.043 | 8/8 / 8/8 |
| t6/conformers_56 | +12.439% | +14.268% | -0.004 | 5/8 / 8/8 |
| t1/ensemble_sdf | +17.511% | +17.475% | +0.006 | 8/8 / 8/8 |
| t6/ensemble_sdf | +16.428% | +16.536% | -0.016 | 8/8 / 8/8 |
| t1/db_build_ensembles | +14.230% | +14.398% | +0.041 | 8/8 / 8/8 |
| t6/db_build_ensembles | +11.362% | +14.264% | -0.020 | 8/8 / 8/8 |
| t1/screen_1000 | +1.564% | +0.873% | +0.033 | 6/8 / 6/8 |
| t6/screen_1000 | -1.319% | -0.065% | +0.045 | 4/8 / 3/8 |
| t1/parse_ensembles_1000 | +1.021% | +0.764% | -0.078 | 6/8 / 6/8 |
| t6/parse_ensembles_1000 | +1.594% | +1.693% | +0.100 | 6/8 / 8/8 |
| t1/search_database | +61.889% | +2.374% | +0.004 | 6/8 / 6/8 |
| t6/search_database | +100.308% | +1.531% | +0.004 | 6/8 / 5/8 |

Conformer, ensemble and database-build medians regress at both thread counts, with all their CPU pairs slower. Parsing has smaller negative wall and CPU medians. Screening is approximately CPU-neutral, with mixed wall results. Search has large wall dispersion and order effects: speedup ranges are 0.139–1.270 at one thread and 0.150–2.370 at six threads. Its negative medians and smaller CPU costs remain reported; no reliable search benefit or particular causal explanation is claimed.

## Complete original-baseline matrix

This table compares C3 with the original corrected baseline, so parallel descriptor/conformer and screening gains also include accepted C1/C2 changes. They cannot satisfy C3’s incremental benefit requirement. Every slower row remains visible. The review receipt provides the full ranges, MAD, order effects and faults for each row.

| Configuration | Median wall speedup | Paired CPU cost | Median RSS delta MiB | Slower wall / CPU pairs |
| --- | ---: | ---: | ---: | ---: |
| t1/conformers_56 | 0.884512 | +13.189% | -0.021 | 8/8 / 8/8 |
| t1/db_build_ensembles | 0.881512 | +13.590% | +0.025 | 8/8 / 8/8 |
| t1/descriptors_10000 | 0.872251 | +14.352% | +0.025 | 8/8 / 8/8 |
| t1/ensemble_sdf | 0.864042 | +15.790% | +0.059 | 8/8 / 8/8 |
| t1/parse_ensembles_1000 | 0.990016 | +1.789% | +0.031 | 5/8 / 6/8 |
| t1/predict_1000000 | 0.996321 | +0.394% | +0.088 | 5/8 / 4/8 |
| t1/screen_1000 | 1.254150 | -18.306% | +0.018 | 0/8 / 0/8 |
| t1/search_database | 1.699440 | +0.954% | +0.002 | 4/8 / 4/8 |
| t1/single_large | 0.957040 | +6.193% | +0.059 | 8/8 / 8/8 |
| t1/single_small | 0.976781 | +1.647% | -0.008 | 6/8 / 6/8 |
| t2/conformers_56 | 1.610723 | +16.115% | +1.006 | 0/8 / 8/8 |
| t2/db_build_ensembles | 0.887465 | +11.958% | -0.012 | 8/8 / 8/8 |
| t2/descriptors_10000 | 1.711667 | +15.893% | +3.740 | 0/8 / 8/8 |
| t2/ensemble_sdf | 0.860356 | +16.379% | +0.098 | 8/8 / 8/8 |
| t2/parse_ensembles_1000 | 0.967600 | +3.365% | +0.025 | 8/8 / 8/8 |
| t2/predict_1000000 | 1.009366 | -0.996% | -0.002 | 2/8 / 2/8 |
| t2/screen_1000 | 1.199197 | -18.142% | +0.051 | 1/8 / 0/8 |
| t2/search_database | 0.943231 | +4.444% | +0.035 | 5/8 / 6/8 |
| t2/single_large | 0.981785 | +3.118% | -0.012 | 6/8 / 7/8 |
| t2/single_small | 0.992797 | +0.347% | -0.002 | 4/8 / 5/8 |
| t4/conformers_56 | 2.786302 | +21.880% | +1.650 | 0/8 / 8/8 |
| t4/db_build_ensembles | 0.929655 | +10.817% | +0.000 | 7/8 / 7/8 |
| t4/descriptors_10000 | 3.286814 | +17.960% | +4.564 | 0/8 / 8/8 |
| t4/ensemble_sdf | 0.870860 | +14.855% | +0.025 | 8/8 / 8/8 |
| t4/parse_ensembles_1000 | 0.969843 | +3.108% | -0.031 | 7/8 / 7/8 |
| t4/predict_1000000 | 1.000955 | +1.520% | -0.002 | 4/8 / 6/8 |
| t4/screen_1000 | 1.258802 | -18.586% | -0.029 | 0/8 / 0/8 |
| t4/search_database | 0.935133 | +2.786% | +0.012 | 4/8 / 7/8 |
| t4/single_large | 0.983914 | +1.191% | +0.035 | 4/8 / 4/8 |
| t4/single_small | 0.948167 | +6.617% | +0.027 | 8/8 / 7/8 |
| t6/conformers_56 | 3.358354 | +26.339% | +2.281 | 0/8 / 8/8 |
| t6/db_build_ensembles | 0.904859 | +11.951% | +0.014 | 7/8 / 8/8 |
| t6/descriptors_10000 | 4.773350 | +22.028% | +5.475 | 0/8 / 8/8 |
| t6/ensemble_sdf | 0.851747 | +16.949% | -0.012 | 8/8 / 8/8 |
| t6/parse_ensembles_1000 | 0.959056 | +4.272% | -0.068 | 8/8 / 8/8 |
| t6/predict_1000000 | 1.044103 | -0.207% | +0.043 | 2/8 / 4/8 |
| t6/screen_1000 | 1.194426 | -17.665% | -0.090 | 0/8 / 0/8 |
| t6/search_database | 0.574751 | +7.133% | -0.018 | 5/8 / 6/8 |
| t6/single_large | 0.945854 | +4.789% | -0.031 | 5/8 / 5/8 |
| t6/single_small | 0.982888 | -0.482% | +0.029 | 5/8 / 4/8 |

Across all launches, including warmups, the full-matrix C3 maximum RSS was 69.699 MiB versus 69.781 MiB for its original-baseline side. Its maximum total child CPU was 6.107 s versus 5.126 s; maximum one-core CPU utilization was 584.1%, and maximum minor faults were 6,407 versus 4,639. In the incremental campaign, C3/C2 maximum RSS was 26.492/26.566 MiB and maximum child CPU was 6.047/5.234 s; maximum minor faults were 6,410/6,411. All 972 launches reported zero major faults. These extrema are retained observations across different workloads, not universal resource-improvement claims.

## Correct operation reduction, negative native result

The separate frozen diagnostic probe checks all 30,168 frame identities, orientations, explicit neighbor orders and 452,520 f32 fields bit for bit. For the 10,000-file workload, actual cached Z predicates are 921,865,331 and newly prepared Z predicates are 51,997,436. Charging 86,487,808 endpoint predicates gives 1,060,350,575 Z predicates versus 1,943,832,492 originally: **45.4505% fewer**. The declared combined XY/Z operation model decreases 33.6009%. Every per-frame count matches the old shadow prediction.

The native candidate also performs full-row bounds/finite scans and helper checks, and changes control flow. The probe records 462,240,000 row-scan point visits and 111,964,217 helper calls on that workload. These counts are not equal-cost hardware instructions; row discovery also replaces an earlier row-change check, so the listed row work is not all net new. The measurements establish that this implementation is slower. They do not isolate the cost of each added check or prove a particular compiler/cache mechanism. The instrumented probe’s timings and allocations are excluded from the native result.

## Scientific result and preserved evidence

Complete scientific admission passed before timing: 128,021 observations across 13 lanes; 93 CLI cases at each of four thread settings and 24 supplementary bridge cases; all eight million-value prediction exports; full independent inferred/topology geometry, 31,721 buried-volume cases, 95,163 private frames/1,141,956 bins, fixed reference union, 7,624 model checks and scoped thermodynamic checks. Engineering passed 325 Rust tests, 124 Python tests, formatting, Clippy, Ruff and documentation checks. The six focused row-endpoint test groups cover ULP/equality, absorbed addition/underflow, overflow/nonfinite fallback, late nonfinite and signed-zero rows, lazy cache/order/dense hits, and 10,000 finite-float witnesses. The unchanged direct reference checks all 15 output fields. Existing scientific-method limitations remain unchanged.

The C3 source, patch and these tests remain in the frozen candidate build and pre-edit/implementation receipts; restoration removes only the unaccepted C3 production/test additions from the live buried-volume file. No scientific expectations were rewritten. Failed tooling evidence also remains: counter-probe attempt v1 could not compile library tests because two test-only C2 fixtures were absent, before any workloads ran; v2 copies their exact archived bytes. The first independent-reference gate retains its six creation-timestamp-derived model-hash mismatches; the additive resolver proves canonical model equality after excluding only `/created/created_utc`, with all 7,624 numerical/status checks unchanged.

All paths below are under `.stericx/profiling/scientifically_validated_optimization/`. These roots retain complete source, binary, input, raw measurement and failure inventories.

| Receipt | SHA-256 |
| --- | --- |
| `candidate_c3_independent_review_v1/review.json` | `68a6ff47ddc58946ef5dc41cc769db89391800f6f097d400e69d66f2e47b222c` |
| `candidate_c3_independent_review_v1/analysis.json` | `852f2fcf9d6ac3afe98f27edd65b8a45bb6bc2674b1a991364aed755e2907172` |
| `candidate_c3_incremental_independent_review_v1/analysis.json` | `f68dc0ebe5db2887d7ef020035074684f79e7246c11601973d6b25ffbe78289f` |
| `candidate_c3_native_pairs_v1/manifest.json` | `daf5d58b442ea5694da9da59a6e2afb228cca360b81ee924d7c8d9ee363ef35e` |
| `candidate_c3_incremental_pairs_v1/manifest.json` | `244066f0b133f47f9e4f568e40fa1d8feb59b46057ee69d877cf1aad4327cc45` |
| `candidate_c3_performance_gate_v1/manifest.json` | `c4d58655c610c5feebee25baf86a12f959b963ed716ecbaf3a74ed81dd2129ff` |
| `candidate_c3_row_rejection/pre_edit_manifest.json` | `33fe1d8df0c9ea2321217225991ada5471503a8c596ad75f3cfa327ce6c172c5` |
| `candidate_c3_row_rejection/implementation_manifest.json` | `2086db59611e6799a119c2223475e695375cf825cbc3f11e0617a09348d2dfe9` |
| `candidate_c3_occupancy_counts_v2/manifest.json` | `485188e4d2aff482453effb3e36e9c2c43856d2ba6960a436daa6485c3dc8887` |
| `candidate_c3_occupancy_counts_v1/failure.json` | `2a904ecd2f98643412920c613dec1848744d7bdffecf40f4b16255046d33288e` |
| `candidate_c3_science_v1/independent_references_v1/independent_gate.json` | `b82ddd205127e0319345e8b9a5815dace5a519894d1d40ec6d43b856894ac082` |
| `candidate_c3_science_v1/independent_reference_resolution_v1/independent_gate.json` | `390870f7432e4bc5b4e49b5a54355bfadc415f66f77562ae3d74deafbfc7f311` |
| `rejected_c3_v1/decision.json` | `b07b0072c926e795bbf3517daf32d2c251ea3d2ae7892082e3ff2f8244f4340a` |
| `rejected_c3_v1/restoration.json` | `f1293e23cfeab06e4062a1ea4e3605d619ef4e6fa28998460e4b51d36899fc91` |
