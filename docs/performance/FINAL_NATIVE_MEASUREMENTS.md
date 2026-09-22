# Final native measurements

The original corrected baseline is compared directly with the final retained
C2 executable after C3 rejection and fresh scientific/reference/engineering
verification. All 720 launches and 320 measured pairs are retained. Each
configuration has one warmup per binary and eight alternating AB/BA pairs.
Inputs, release settings and affinities are identical within each pair.

Times below include startup and output, in milliseconds. MAD is median
absolute deviation, not a confidence interval. Speedup is the median paired
baseline/candidate ratio, not the ratio of the two medians. Every slower row
and every outlier remains in [the complete reviewed data](current_final_measurements.json).
The scientific fingerprints match in every launch; independent accuracy
limits remain those in the separate scientific-equivalence report.

| Configuration | Baseline median ± MAD ms | Final median ± MAD ms | Baseline range ms | Final range ms | Paired speedup ± MAD | Paired full range | Faster pairs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `t1/conformers_56` | 29.078 ± 0.049 | 28.898 ± 0.046 | 28.930–29.227 | 28.840–29.344 | 1.005× ± 0.004 | 0.986–1.013 | 6/8 |
| `t1/db_build_ensembles` | 30.716 ± 0.216 | 29.907 ± 0.161 | 29.778–30.958 | 29.424–30.551 | 1.021× ± 0.008 | 1.001–1.030 | 8/8 |
| `t1/descriptors_10000` | 4935.780 ± 9.907 | 4952.589 ± 67.567 | 4923.540–5043.049 | 4870.866–5063.421 | 1.001× ± 0.009 | 0.975–1.025 | 4/8 |
| `t1/ensemble_sdf` | 46.250 ± 0.210 | 46.205 ± 0.316 | 46.036–47.180 | 45.789–47.991 | 1.001× ± 0.007 | 0.959–1.016 | 5/8 |
| `t1/parse_ensembles_1000` | 184.256 ± 1.084 | 187.886 ± 0.674 | 182.675–185.795 | 186.177–212.061 | 0.983× ± 0.008 | 0.865–0.992 | 0/8 |
| `t1/predict_1000000` | 20.987 ± 0.101 | 20.880 ± 0.108 | 20.807–21.520 | 20.516–21.022 | 1.009× ± 0.009 | 0.994–1.028 | 6/8 |
| `t1/screen_1000` | 255.628 ± 6.846 | 207.645 ± 6.358 | 248.548–280.123 | 201.262–224.368 | 1.204× ± 0.083 | 1.108–1.385 | 8/8 |
| `t1/search_database` | 7.707 ± 0.963 | 7.459 ± 0.488 | 6.716–28.357 | 6.956–23.272 | 1.060× ± 0.726 | 0.293–4.077 | 4/8 |
| `t1/single_large` | 2.069 ± 0.030 | 2.101 ± 0.010 | 2.019–3.701 | 2.010–2.161 | 0.998× ± 0.035 | 0.960–1.769 | 4/8 |
| `t1/single_small` | 2.103 ± 0.052 | 2.123 ± 0.065 | 1.944–2.208 | 1.976–2.206 | 1.008× ± 0.021 | 0.881–1.066 | 5/8 |
| `t2/conformers_56` | 29.190 ± 0.248 | 16.026 ± 0.087 | 28.939–29.737 | 15.880–16.242 | 1.817× ± 0.027 | 1.785–1.873 | 8/8 |
| `t2/db_build_ensembles` | 17.110 ± 0.352 | 16.788 ± 0.050 | 16.685–17.766 | 16.448–214.044 | 1.014× ± 0.026 | 0.080–1.058 | 6/8 |
| `t2/descriptors_10000` | 4965.442 ± 7.594 | 2548.986 ± 26.050 | 4955.254–5064.075 | 2500.108–2579.245 | 1.954× ± 0.015 | 1.924–1.982 | 8/8 |
| `t2/ensemble_sdf` | 46.701 ± 0.162 | 46.306 ± 0.174 | 46.373–48.206 | 46.012–46.920 | 1.008× ± 0.006 | 0.997–1.047 | 7/8 |
| `t2/parse_ensembles_1000` | 184.323 ± 0.673 | 188.050 ± 0.378 | 183.308–198.618 | 186.432–194.192 | 0.982× ± 0.002 | 0.953–1.057 | 1/8 |
| `t2/predict_1000000` | 13.651 ± 0.149 | 13.693 ± 0.069 | 13.402–13.984 | 13.544–13.816 | 1.001× ± 0.016 | 0.973–1.022 | 4/8 |
| `t2/screen_1000` | 254.349 ± 4.636 | 214.041 ± 9.176 | 248.508–274.814 | 201.588–223.556 | 1.229× ± 0.059 | 1.116–1.309 | 8/8 |
| `t2/search_database` | 6.913 ± 0.095 | 17.627 ± 6.731 | 6.791–27.790 | 7.132–25.096 | 0.887× ± 0.566 | 0.290–2.239 | 2/8 |
| `t2/single_large` | 2.148 ± 0.029 | 2.156 ± 0.066 | 2.078–2.258 | 2.031–2.767 | 1.004× ± 0.050 | 0.785–1.071 | 4/8 |
| `t2/single_small` | 2.146 ± 0.038 | 2.152 ± 0.037 | 2.039–2.189 | 2.059–2.573 | 0.973× ± 0.022 | 0.831–1.046 | 2/8 |
| `t4/conformers_56` | 29.719 ± 0.237 | 9.520 ± 0.108 | 29.318–30.431 | 9.250–11.992 | 3.132× ± 0.046 | 2.478–3.241 | 8/8 |
| `t4/db_build_ensembles` | 10.939 ± 0.796 | 10.095 ± 0.087 | 9.959–14.374 | 9.825–12.705 | 1.036× ± 0.052 | 0.982–1.428 | 5/8 |
| `t4/descriptors_10000` | 4952.194 ± 14.971 | 1305.450 ± 7.535 | 4932.137–4981.562 | 1297.631–1336.612 | 3.784× ± 0.011 | 3.727–3.820 | 8/8 |
| `t4/ensemble_sdf` | 46.750 ± 0.165 | 46.609 ± 0.484 | 46.352–47.006 | 45.850–47.896 | 1.004× ± 0.009 | 0.981–1.014 | 5/8 |
| `t4/parse_ensembles_1000` | 182.949 ± 0.671 | 186.568 ± 0.613 | 182.277–186.033 | 185.928–195.647 | 0.981× ± 0.008 | 0.932–0.994 | 0/8 |
| `t4/predict_1000000` | 10.078 ± 0.098 | 10.082 ± 0.126 | 9.950–10.665 | 9.941–11.578 | 0.992× ± 0.007 | 0.867–1.073 | 2/8 |
| `t4/screen_1000` | 247.898 ± 0.987 | 206.913 ± 6.104 | 246.486–267.421 | 199.509–221.415 | 1.201× ± 0.059 | 1.119–1.323 | 8/8 |
| `t4/search_database` | 10.432 ± 3.234 | 10.998 ± 3.851 | 7.020–25.730 | 7.127–28.465 | 1.102× ± 0.746 | 0.256–3.610 | 4/8 |
| `t4/single_large` | 2.254 ± 0.054 | 2.292 ± 0.016 | 2.192–2.445 | 2.129–2.500 | 1.002× ± 0.029 | 0.877–1.129 | 4/8 |
| `t4/single_small` | 2.169 ± 0.039 | 2.227 ± 0.054 | 2.128–2.240 | 2.117–18.836 | 0.973× ± 0.034 | 0.113–1.049 | 3/8 |
| `t6/conformers_56` | 29.607 ± 0.209 | 7.161 ± 0.133 | 29.196–30.866 | 6.960–9.854 | 4.188× ± 0.095 | 3.007–4.304 | 8/8 |
| `t6/db_build_ensembles` | 8.053 ± 0.372 | 8.335 ± 0.509 | 7.519–19.844 | 7.305–12.486 | 1.004× ± 0.098 | 0.895–1.589 | 4/8 |
| `t6/descriptors_10000` | 4969.625 ± 2.637 | 911.846 ± 12.652 | 4954.735–4989.818 | 893.332–1286.438 | 5.462× ± 0.061 | 3.852–5.562 | 8/8 |
| `t6/ensemble_sdf` | 46.994 ± 0.201 | 46.355 ± 0.051 | 46.359–48.716 | 46.298–47.515 | 1.011× ± 0.005 | 0.976–1.047 | 7/8 |
| `t6/parse_ensembles_1000` | 184.547 ± 1.323 | 187.960 ± 1.964 | 182.644–188.555 | 185.261–191.093 | 0.986× ± 0.015 | 0.960–1.010 | 2/8 |
| `t6/predict_1000000` | 9.815 ± 0.165 | 8.949 ± 0.162 | 8.614–10.540 | 8.739–10.131 | 1.068× ± 0.046 | 0.863–1.128 | 7/8 |
| `t6/screen_1000` | 250.585 ± 2.147 | 206.862 ± 6.434 | 248.038–267.240 | 200.405–221.700 | 1.206× ± 0.047 | 1.132–1.333 | 8/8 |
| `t6/search_database` | 7.511 ± 0.412 | 10.227 ± 2.937 | 7.075–41.440 | 7.152–26.525 | 0.969× ± 0.585 | 0.349–5.597 | 3/8 |
| `t6/single_large` | 2.231 ± 0.066 | 2.183 ± 0.028 | 2.115–2.474 | 2.115–24.282 | 1.011× ± 0.044 | 0.095–1.170 | 5/8 |
| `t6/single_small` | 2.311 ± 0.122 | 2.265 ± 0.018 | 2.180–14.521 | 2.242–2.416 | 1.008× ± 0.061 | 0.942–6.451 | 4/8 |

## Resources and launch order

CPU cost is the median paired percent change; positive means more CPU time.
RSS columns show the median child's peak and its maximum among eight measured
samples. CPU utilization 100% is one core. Raw warmup resources, user/system
CPU, faults and context switches remain in the linked original capture.
Allocation counts and their different scope are in the C2 diagnostic reprofile.

| Configuration | AB / BA speedup | Paired CPU cost | Final CPU utilization % | Baseline RSS median / max MiB | Final RSS median / max MiB | Final minor / major faults median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `t1/conformers_56` | 1.000 / 1.005 | -0.46% | 99.3 | 5.09 / 5.16 | 5.17 / 5.32 | 302.5 / 0 |
| `t1/db_build_ensembles` | 1.019 / 1.022 | -2.00% | 99.3 | 5.53 / 5.64 | 5.54 / 5.70 | 357.5 / 0 |
| `t1/descriptors_10000` | 1.001 / 1.003 | -0.02% | 99.7 | 20.99 / 21.07 | 20.99 / 21.13 | 4632.5 / 0 |
| `t1/ensemble_sdf` | 1.001 / 0.999 | -0.38% | 99.4 | 5.84 / 5.98 | 5.85 / 5.89 | 593 / 0 |
| `t1/parse_ensembles_1000` | 0.978 / 0.988 | +1.75% | 99.8 | 4.93 / 5.02 | 4.93 / 5.09 | 218.5 / 0 |
| `t1/predict_1000000` | 1.016 / 1.003 | -0.87% | 99.3 | 69.54 / 69.61 | 69.61 / 69.65 | 2150 / 0 |
| `t1/screen_1000` | 1.204 / 1.234 | -19.01% | 97.6 | 10.26 / 10.34 | 10.22 / 10.30 | 2304.5 / 0 |
| `t1/search_database` | 0.757 / 2.443 | +6.20% | 96.1 | 7.82 / 7.90 | 7.78 / 7.91 | 964.5 / 0 |
| `t1/single_large` | 0.998 / 0.997 | +2.80% | 90.0 | 5.03 / 5.16 | 5.04 / 5.20 | 239.5 / 0 |
| `t1/single_small` | 1.018 / 0.995 | +1.48% | 90.0 | 5.04 / 5.09 | 5.08 / 5.16 | 241 / 0 |
| `t2/conformers_56` | 1.846 / 1.797 | +2.92% | 183.9 | 5.15 / 5.18 | 6.07 / 6.15 | 520.5 / 0 |
| `t2/db_build_ensembles` | 1.024 / 1.014 | -1.30% | 179.3 | 5.96 / 6.08 | 5.81 / 5.98 | 507 / 0 |
| `t2/descriptors_10000` | 1.962 / 1.952 | +1.01% | 196.4 | 20.92 / 20.97 | 24.79 / 24.91 | 5936 / 0 |
| `t2/ensemble_sdf` | 1.016 / 1.006 | -0.84% | 99.5 | 5.81 / 5.93 | 5.77 / 6.00 | 591.5 / 0 |
| `t2/parse_ensembles_1000` | 0.980 / 0.982 | +1.95% | 99.8 | 4.79 / 5.02 | 5.01 / 5.17 | 220.5 / 0 |
| `t2/predict_1000000` | 1.008 / 0.989 | +0.33% | 155.4 | 69.46 / 69.70 | 69.60 / 69.70 | 2155 / 0 |
| `t2/screen_1000` | 1.203 / 1.229 | -18.68% | 96.2 | 10.15 / 10.34 | 10.14 / 10.22 | 2305 / 0 |
| `t2/search_database` | 0.887 / 0.648 | +6.49% | 42.0 | 7.81 / 7.84 | 7.77 / 7.87 | 964 / 0 |
| `t2/single_large` | 1.019 / 1.004 | -0.37% | 91.2 | 5.00 / 5.07 | 5.05 / 5.16 | 239 / 0 |
| `t2/single_small` | 0.991 / 0.953 | +2.31% | 89.7 | 5.06 / 5.15 | 5.08 / 5.18 | 240.5 / 0 |
| `t4/conformers_56` | 3.082 / 3.177 | +6.42% | 333.1 | 5.10 / 5.25 | 6.75 / 6.91 | 818.5 / 0 |
| `t4/db_build_ensembles` | 1.007 / 1.078 | -0.90% | 314.9 | 6.56 / 6.77 | 6.64 / 6.74 | 804.5 / 0 |
| `t4/descriptors_10000` | 3.785 / 3.784 | +2.58% | 387.0 | 20.80 / 21.01 | 25.35 / 25.63 | 6155 / 0 |
| `t4/ensemble_sdf` | 1.003 / 1.004 | -0.31% | 99.5 | 5.82 / 5.93 | 5.79 / 5.90 | 592 / 0 |
| `t4/parse_ensembles_1000` | 0.981 / 0.977 | +1.97% | 99.9 | 4.78 / 4.93 | 4.76 / 4.92 | 218.5 / 0 |
| `t4/predict_1000000` | 0.992 / 1.001 | -0.10% | 220.6 | 69.41 / 69.55 | 69.34 / 69.57 | 2172 / 0 |
| `t4/screen_1000` | 1.184 / 1.243 | -18.86% | 97.5 | 10.14 / 10.29 | 10.14 / 10.34 | 2304.5 / 0 |
| `t4/search_database` | 2.086 / 0.356 | +5.00% | 74.3 | 7.69 / 7.83 | 7.76 / 7.89 | 964 / 0 |
| `t4/single_large` | 1.027 / 0.973 | -0.33% | 90.7 | 5.05 / 5.14 | 5.07 / 5.12 | 239.5 / 0 |
| `t4/single_small` | 0.930 / 1.007 | +1.27% | 89.6 | 5.02 / 5.17 | 5.03 / 5.09 | 241 / 0 |
| `t6/conformers_56` | 4.188 / 4.143 | +13.73% | 462.3 | 5.12 / 5.31 | 7.35 / 7.68 | 1118 / 0 |
| `t6/db_build_ensembles` | 0.982 / 1.063 | -2.95% | 399.3 | 7.26 / 7.49 | 7.27 / 7.38 | 1101 / 0 |
| `t6/descriptors_10000` | 5.481 / 5.445 | +5.05% | 570.3 | 20.78 / 21.14 | 26.35 / 26.48 | 6398.5 / 0 |
| `t6/ensemble_sdf` | 1.012 / 1.011 | -1.06% | 99.4 | 5.77 / 6.00 | 5.84 / 5.94 | 592 / 0 |
| `t6/parse_ensembles_1000` | 0.980 / 0.997 | +1.48% | 99.9 | 4.83 / 5.05 | 4.78 / 4.83 | 218.5 / 0 |
| `t6/predict_1000000` | 1.038 / 1.087 | +0.44% | 262.9 | 69.43 / 69.64 | 69.41 / 69.66 | 2190 / 0 |
| `t6/screen_1000` | 1.206 / 1.206 | -18.68% | 97.7 | 10.04 / 10.27 | 10.09 / 10.25 | 2304.5 / 0 |
| `t6/search_database` | 0.762 / 1.802 | +3.20% | 76.0 | 7.73 / 7.80 | 7.78 / 7.82 | 963.5 / 0 |
| `t6/single_large` | 0.951 / 1.030 | -1.85% | 89.7 | 5.00 / 5.13 | 5.02 / 5.16 | 240 / 0 |
| `t6/single_small` | 1.415 / 0.969 | -1.12% | 89.7 | 5.00 / 5.20 | 5.06 / 5.23 | 240 / 0 |

## Reproduction and interpretation

The analysis binds the final native manifest, raw streams/artifacts, child
resource records, scientific gate, full pre/postflight identities and the
unchanged independent review helper. No scientifically obsolete binary is used as
the scientific denominator. The final main report identifies accepted
workload tradeoffs and avoids a precise speedup claim for noisy search.
Desktop activity and launch-order effects remain possible; these are fixed
workload observations on one host, not cross-machine promises or prospective
chemical validation. Offline artifact verification is not a timing rerun.

Analysis SHA-256: `a428b351d87e4454314230875b02f2696448371e2c61d4c0aa707b4261895757`.
