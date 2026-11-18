# Final controlled measurements

All trials, CPU/wall time, RSS, allocation diagnostics, MAD and ranges are in [final_measurements.json](final_measurements.json). Warmups remain in the raw receipts. No outliers were removed. Fresh processes use warmed file caches.

## 10,000-file throughput

| Threads | Baseline files/s | Final files/s | Speedup | Final efficiency |
|---|---|---|---|---|
| 1 | 4,490.1 | 6,680.4 | 1.4878× | 100.00% |
| 2 | 8,701.9 | 12,766.3 | 1.4671× | 95.55% |
| 4 | 16,566.2 | 24,088.5 | 1.4541× | 90.15% |
| 6 | 23,411.7 | 34,296.4 | 1.4649× | 85.56% |
| 12 | 28,413.9 | 41,926.5 | 1.4756× | 52.30% |

## Paired kernel latency

| Mode | Baseline µs | Final µs | Final MAD µs | Speedup |
|---|---|---|---|---|
| volume | 166.737 | 81.258 | 0.078 | 2.0520× |
| regions | 166.277 | 81.135 | 0.074 | 2.0494× |
| ensemble | 1124.660 | 681.855 | 0.843 | 1.6494× |
| radius_profile | 758.524 | 327.317 | 0.422 | 2.3174× |
| sterimol | 12.472 | 12.498 | 0.005 | 0.9980× |

## Whole-process workload distributions

| Workload/threads | Baseline ms | Final ms | Final MAD ms | Final range ms | Speedup |
|---|---|---|---|---|---|
| single_small/t1 | 1.5625 | 1.6466 | 0.0453 | 1.5001–1.6963 | 0.9489× |
| single_large/t1 | 1.6174 | 1.5884 | 0.0295 | 1.5560–1.7696 | 1.0183× |
| conformers_56/t1 | 13.9553 | 9.9997 | 0.0592 | 9.9039–10.2289 | 1.3956× |
| ensemble_sdf/t1 | 23.7867 | 17.8033 | 0.0457 | 17.7213–19.7305 | 1.3361× |
| descriptors_10000/t1 | 2227.1334 | 1496.9120 | 1.9064 | 1494.1740–1505.3227 | 1.4878× |
| parse_ensembles_1000/t1 | 186.4113 | 186.4176 | 0.3449 | 185.9404–191.5590 | 1.0000× |
| db_build_ensembles/t1 | 14.6430 | 10.6698 | 0.0894 | 10.5469–12.1719 | 1.3724× |
| single_small/t2 | 1.6278 | 1.6081 | 0.0144 | 1.5811–1.7905 | 1.0123× |
| single_large/t2 | 1.6708 | 1.6576 | 0.0144 | 1.6333–2.1953 | 1.0079× |
| conformers_56/t2 | 8.1407 | 6.2279 | 0.0652 | 6.1289–6.6089 | 1.3071× |
| ensemble_sdf/t2 | 23.7635 | 17.7709 | 0.0333 | 17.7303–18.1974 | 1.3372× |
| descriptors_10000/t2 | 1149.1725 | 783.3149 | 4.2870 | 774.9609–792.5456 | 1.4671× |
| parse_ensembles_1000/t2 | 187.1927 | 187.2086 | 0.2872 | 186.5853–201.0878 | 0.9999× |
| db_build_ensembles/t2 | 8.8457 | 6.8542 | 0.1199 | 6.6122–7.0211 | 1.2905× |
| single_small/t4 | 1.7591 | 1.7302 | 0.0269 | 1.7006–1.8910 | 1.0167× |
| single_large/t4 | 1.7901 | 1.8066 | 0.0297 | 1.7429–1.8710 | 0.9909× |
| conformers_56/t4 | 5.2440 | 4.2910 | 0.0620 | 4.1124–4.4394 | 1.2221× |
| ensemble_sdf/t4 | 24.1666 | 18.2937 | 0.3334 | 17.8198–65.5252 | 1.3210× |
| descriptors_10000/t4 | 603.6372 | 415.1352 | 2.8485 | 409.9405–424.2110 | 1.4541× |
| parse_ensembles_1000/t4 | 186.5367 | 186.3350 | 0.1907 | 185.7050–187.0174 | 1.0011× |
| db_build_ensembles/t4 | 5.8277 | 4.9383 | 0.1198 | 4.8056–7.3872 | 1.1801× |
| single_small/t6 | 1.7418 | 1.7632 | 0.0285 | 1.6911–1.8879 | 0.9878× |
| single_large/t6 | 1.9030 | 1.8403 | 0.0453 | 1.7738–7.9351 | 1.0341× |
| conformers_56/t6 | 4.2228 | 3.5082 | 0.0794 | 3.3451–3.7401 | 1.2037× |
| ensemble_sdf/t6 | 24.3596 | 18.2320 | 0.2538 | 17.9166–18.6881 | 1.3361× |
| descriptors_10000/t6 | 427.1377 | 291.5758 | 4.1848 | 285.5098–411.7208 | 1.4649× |
| parse_ensembles_1000/t6 | 186.4381 | 186.9620 | 0.6362 | 185.3593–188.0462 | 0.9972× |
| db_build_ensembles/t6 | 4.8494 | 4.1125 | 0.1261 | 3.9095–36.0604 | 1.1792× |
| single_small/t12 | 1.9011 | 1.8561 | 0.0333 | 1.6486–1.9515 | 1.0242× |
| single_large/t12 | 1.9170 | 1.8645 | 0.0339 | 1.8006–1.9888 | 1.0281× |
| conformers_56/t12 | 4.2929 | 3.7098 | 0.1132 | 3.5807–3.9360 | 1.1572× |
| ensemble_sdf/t12 | 24.4872 | 18.3126 | 0.3825 | 17.7757–45.5749 | 1.3372× |
| descriptors_10000/t12 | 351.9402 | 238.5126 | 1.2487 | 235.2214–365.8811 | 1.4756× |
| parse_ensembles_1000/t12 | 186.7689 | 186.9849 | 0.5876 | 185.9661–191.9579 | 0.9988× |
| db_build_ensembles/t12 | 4.6909 | 4.1608 | 0.0686 | 3.9135–4.2705 | 1.1274× |

## competitive

| Size/workers/scope/tool | Median ms | MAD ms | Range ms |
|---|---|---|---|
| n1/t1/fresh_process/native | 1.23823 | 0.08682 | 1.14657–1.43908 |
| n1/t1/fresh_process/morfeus | 286.33656 | 0.46363 | 285.26068–287.78352 |
| n1/t1/fresh_process/dbstep | 284.89819 | 0.71584 | 284.03365–289.66283 |
| n1/t1/warm_workflow/native | 0.08283 | 0.00014 | 0.08261–0.08445 |
| n1/t1/warm_workflow/morfeus | 14.12611 | 0.02961 | 14.09511–14.18441 |
| n1/t1/warm_workflow/dbstep | 13.61238 | 0.02762 | 13.56903–13.67307 |
| n56/t1/fresh_process/native | 9.54710 | 0.05460 | 9.26089–10.22620 |
| n56/t1/fresh_process/morfeus | 1027.42081 | 9.17045 | 993.88413–1081.40484 |
| n56/t1/fresh_process/dbstep | 908.69224 | 2.24475 | 904.74991–919.84589 |
| n56/t1/warm_workflow/native | 7.49220 | 0.01171 | 7.43119–7.51296 |
| n56/t1/warm_workflow/morfeus | 751.37992 | 4.96017 | 721.85183–814.98410 |
| n56/t1/warm_workflow/dbstep | 633.19778 | 3.63565 | 628.21389–687.94216 |
| n56/t2/fresh_process/native | 6.01094 | 0.06163 | 5.92286–6.38356 |
| n56/t2/fresh_process/morfeus | 1019.42575 | 1.46360 | 1011.60269–1033.50773 |
| n56/t2/fresh_process/dbstep | 976.12152 | 11.23960 | 961.37215–1007.03640 |
| n56/t2/warm_workflow/native | 4.01239 | 0.02048 | 3.96184–5.21380 |
| n56/t2/warm_workflow/morfeus | 413.54436 | 3.35664 | 398.32708–425.49514 |
| n56/t2/warm_workflow/dbstep | 371.24897 | 4.31312 | 357.59674–393.40774 |
| n56/t4/fresh_process/native | 4.35649 | 0.10506 | 4.22999–4.56217 |
| n56/t4/fresh_process/morfeus | 879.92284 | 2.26204 | 874.70613–897.95523 |
| n56/t4/fresh_process/dbstep | 856.64653 | 2.08787 | 850.41119–947.38581 |
| n56/t4/warm_workflow/native | 2.00186 | 0.00571 | 1.98437–2.03712 |
| n56/t4/warm_workflow/morfeus | 239.01537 | 5.79448 | 230.10598–244.96013 |
| n56/t4/warm_workflow/dbstep | 216.04916 | 1.79611 | 213.95192–223.59330 |
| n56/t6/fresh_process/native | 3.67991 | 0.10344 | 3.56056–5.04600 |
| n56/t6/fresh_process/morfeus | 848.02625 | 4.68180 | 840.27849–902.72240 |
| n56/t6/fresh_process/dbstep | 822.91428 | 2.84921 | 819.65061–836.51305 |
| n56/t6/warm_workflow/native | 1.60042 | 0.05694 | 1.39624–1.65783 |
| n56/t6/warm_workflow/morfeus | 172.02767 | 2.57668 | 163.35280–177.19651 |
| n56/t6/warm_workflow/dbstep | 152.92645 | 2.89903 | 149.34480–157.32673 |
| n56/t12/fresh_process/native | 3.80674 | 0.05079 | 3.64476–3.94787 |
| n56/t12/fresh_process/morfeus | 1032.91613 | 1.54623 | 1027.84527–1035.58434 |
| n56/t12/fresh_process/dbstep | 1024.34713 | 1.18805 | 1021.56535–1034.06384 |
| n56/t12/warm_workflow/native | 1.18571 | 0.01477 | 1.14466–1.20051 |
| n56/t12/warm_workflow/morfeus | 135.90097 | 0.89942 | 134.30144–138.12115 |
| n56/t12/warm_workflow/dbstep | 126.98765 | 0.62751 | 125.24511–128.35587 |
| n1000/t1/fresh_process/native | 155.60906 | 8.67166 | 146.81617–361.11355 |
| n1000/t1/fresh_process/morfeus | 14202.72471 | 122.33676 | 13017.10811–14559.82871 |
| n1000/t1/fresh_process/dbstep | 12680.67179 | 51.24246 | 12604.00213–12746.03018 |
| n1000/t1/warm_workflow/native | 132.72347 | 0.03460 | 132.54357–133.11322 |
| n1000/t1/warm_workflow/morfeus | 13920.80874 | 118.28626 | 13484.15252–14167.36941 |
| n1000/t1/warm_workflow/dbstep | 12353.99073 | 14.46578 | 12332.14470–12448.21441 |
| n1000/t6/fresh_process/native | 37.72136 | 0.21519 | 33.59071–53.67418 |
| n1000/t6/fresh_process/morfeus | 3489.69896 | 38.78041 | 3419.76560–3637.38980 |
| n1000/t6/fresh_process/dbstep | 3296.66029 | 14.27984 | 3280.13576–3352.55356 |
| n1000/t6/warm_workflow/native | 23.63384 | 0.02634 | 23.43917–24.16940 |
| n1000/t6/warm_workflow/morfeus | 2776.24468 | 10.14345 | 2751.83745–2789.07293 |
| n1000/t6/warm_workflow/dbstep | 2598.37904 | 9.54629 | 2498.68239–2631.65622 |

## radius competitive

| Size/workers/scope/tool | Median ms | MAD ms | Range ms |
|---|---|---|---|
| n1/t1/fresh_process/native | 1.84547 | 0.04241 | 1.78595–1.91212 |
| n1/t1/fresh_process/morfeus | 315.65452 | 0.44915 | 313.97577–316.34205 |
| n1/t1/fresh_process/dbstep | 315.27341 | 1.96301 | 312.26725–322.68436 |
| n1/t1/warm_workflow/native | 0.25976 | 0.00116 | 0.25242–0.26095 |
| n1/t1/warm_workflow/morfeus | 39.74249 | 0.22200 | 39.48443–43.04969 |
| n1/t1/warm_workflow/dbstep | 40.36303 | 0.10384 | 38.24032–40.50375 |
| n12/t1/fresh_process/native | 6.95536 | 0.08823 | 6.81187–7.20922 |
| n12/t1/fresh_process/morfeus | 754.14356 | 2.28582 | 737.48754–761.25763 |
| n12/t1/fresh_process/dbstep | 716.59860 | 6.99808 | 708.17164–726.50081 |
| n12/t1/warm_workflow/native | 4.72808 | 0.01588 | 4.69046–4.76555 |
| n12/t1/warm_workflow/morfeus | 479.67950 | 2.96646 | 464.10370–491.55849 |
| n12/t1/warm_workflow/dbstep | 448.44414 | 1.36044 | 437.18333–450.24651 |
| n12/t6/fresh_process/native | 3.47590 | 0.07946 | 3.36881–3.69874 |
| n12/t6/fresh_process/morfeus | 761.84610 | 4.45094 | 750.93578–776.71974 |
| n12/t6/fresh_process/dbstep | 750.39889 | 1.68947 | 745.30657–752.23733 |
| n12/t6/warm_workflow/native | 0.95474 | 0.02536 | 0.91007–1.03930 |
| n12/t6/warm_workflow/morfeus | 103.83460 | 3.97632 | 98.49166–114.34747 |
| n12/t6/warm_workflow/dbstep | 97.83596 | 1.72065 | 93.74943–105.93457 |

## Accuracy and computational cost

Each row identifies its reference class in the linked JSON. An observed error against a fine same-algorithm molecular enclosure is not independent exact truth. Midpoint methods have no intrinsic certificate.

| Case/method | Median ms | MAD ms | Certified error Å³ | Observed error upper Å³ | Target met |
|---|---|---|---|---|---|
| methylphosphine/reference/certified_0.1 | 6.1234 | 0.0057 | 0.0999933 | 0.0267833 | True |
| methylphosphine/reference/certified_0.05 | 11.3758 | 0.0534 | 0.0499996 | 0.0153412 | True |
| methylphosphine/reference/certified_0.01 | 46.4100 | 0.0590 | 0.00999978 | 0.00420574 | True |
| methylphosphine/reference/reference_262143 | 42.7760 | 0.0613 | 1.36589 | 0.00852318 | False |
| methylphosphine/reference/reference_1048575 | 182.0771 | 0.9852 | 0.68354 | 0.00301086 | False |
| methylphosphine/reference/midpoint_128 | 0.1069 | 0.0003 | uncertified | 0.0014301 | None |
| methylphosphine/reference/midpoint_256 | 0.4178 | 0.0091 | uncertified | 0.00313962 | None |
| methylphosphine/reference/midpoint_512 | 1.6283 | 0.0227 | uncertified | 0.00123117 | None |
| dimethylphosphine/reference/certified_0.1 | 7.6411 | 0.0170 | 0.0999922 | 0.0280089 | True |
| dimethylphosphine/reference/certified_0.05 | 13.4797 | 0.0664 | 0.0499957 | 0.0157642 | True |
| dimethylphosphine/reference/certified_0.01 | 51.8116 | 0.1015 | 0.00999982 | 0.00429447 | True |
| dimethylphosphine/reference/reference_262143 | 48.9205 | 0.1085 | 1.56286 | 0.00760413 | False |
| dimethylphosphine/reference/reference_1048575 | 205.2971 | 0.1989 | 0.782342 | 0.00273454 | False |
| dimethylphosphine/reference/midpoint_128 | 0.1344 | 0.0004 | uncertified | 0.00248933 | None |
| dimethylphosphine/reference/midpoint_256 | 0.5369 | 0.0265 | uncertified | 0.00126778 | None |
| dimethylphosphine/reference/midpoint_512 | 2.0060 | 0.0062 | uncertified | 0.00134063 | None |
| trimethylphosphine/reference/certified_0.1 | 10.3333 | 0.0227 | 0.0999968 | 0.0276551 | True |
| trimethylphosphine/reference/certified_0.05 | 18.4992 | 0.0480 | 0.0499961 | 0.015536 | True |
| trimethylphosphine/reference/certified_0.01 | 71.2113 | 0.1227 | 0.00999997 | 0.0041759 | True |
| trimethylphosphine/reference/reference_262143 | 54.1140 | 0.1478 | 1.83721 | 0.0159511 | False |
| trimethylphosphine/reference/reference_1048575 | 226.7318 | 0.5716 | 0.918697 | 0.00499031 | False |
| trimethylphosphine/reference/midpoint_128 | 0.1641 | 0.0002 | uncertified | 0.0025912 | None |
| trimethylphosphine/reference/midpoint_256 | 0.6272 | 0.0085 | uncertified | 0.00131707 | None |
| trimethylphosphine/reference/midpoint_512 | 2.4448 | 0.0103 | uncertified | 0.00162138 | None |
| triethylphosphine/reference/certified_0.1 | 18.2710 | 0.0406 | 0.0999967 | 0.0267025 | True |
| triethylphosphine/reference/certified_0.05 | 29.9665 | 0.2036 | 0.0499977 | 0.0146692 | True |
| triethylphosphine/reference/certified_0.01 | 107.0425 | 0.3682 | 0.00999985 | 0.00413148 | True |
| triethylphosphine/reference/reference_262143 | 71.3154 | 0.2631 | 2.55857 | 0.0110839 | False |
| triethylphosphine/reference/reference_1048575 | 292.7526 | 1.2872 | 1.28096 | 0.00354248 | False |
| triethylphosphine/reference/midpoint_128 | 0.2644 | 0.0018 | uncertified | 0.0013023 | None |
| triethylphosphine/reference/midpoint_256 | 1.0285 | 0.0153 | uncertified | 0.0018246 | None |
| triethylphosphine/reference/midpoint_512 | 3.8471 | 0.0125 | uncertified | 0.00151106 | None |
| triisopropylphosphine/reference/certified_0.1 | 27.0270 | 0.0367 | 0.0999986 | 0.0245843 | True |
| triisopropylphosphine/reference/certified_0.05 | 44.3754 | 0.1107 | 0.0499996 | 0.0136272 | True |
| triisopropylphosphine/reference/certified_0.01 | 164.0138 | 0.4807 | 0.00999998 | 0.00399239 | True |
| triisopropylphosphine/reference/reference_262143 | 87.0782 | 0.3523 | 3.06519 | 0.0127368 | False |
| triisopropylphosphine/reference/reference_1048575 | 356.4006 | 2.2432 | 1.53561 | 0.00148723 | False |
| triisopropylphosphine/reference/midpoint_128 | 0.3556 | 0.0033 | uncertified | 0.00173602 | None |
| triisopropylphosphine/reference/midpoint_256 | 1.3614 | 0.0261 | uncertified | 0.00147176 | None |
| triisopropylphosphine/reference/midpoint_512 | 5.1486 | 0.0180 | uncertified | 0.00147162 | None |
| triphenylphosphine/reference/certified_0.1 | 30.3644 | 0.0382 | 0.0999993 | 0.0210113 | True |
| triphenylphosphine/reference/certified_0.05 | 51.4635 | 0.1259 | 0.0499995 | 0.0118389 | True |
| triphenylphosphine/reference/certified_0.01 | 185.5922 | 0.2867 | 0.00999994 | 0.00364061 | True |
| triphenylphosphine/reference/reference_262143 | 139.0516 | 0.4618 | 2.70708 | 0.0142362 | False |
| triphenylphosphine/reference/reference_1048575 | 561.0344 | 1.2024 | 1.35517 | 0.00336997 | False |
| triphenylphosphine/reference/midpoint_128 | 0.5239 | 0.0212 | uncertified | 0.00242333 | None |
| triphenylphosphine/reference/midpoint_256 | 1.8732 | 0.0053 | uncertified | 0.00201689 | None |
| triphenylphosphine/reference/midpoint_512 | 7.2885 | 0.0317 | uncertified | 0.00136432 | None |
| asymmetric_phosphine/reference/certified_0.1 | 40.3091 | 0.0756 | 0.0999996 | 0.0163265 | True |
| asymmetric_phosphine/reference/certified_0.05 | 72.6622 | 0.0924 | 0.0499988 | 0.00963343 | True |
| asymmetric_phosphine/reference/certified_0.01 | 289.9533 | 0.4598 | 0.00999999 | 0.00316793 | True |
| asymmetric_phosphine/reference/reference_262143 | 88.3175 | 0.1316 | 2.71514 | 0.0209938 | False |
| asymmetric_phosphine/reference/reference_1048575 | 361.7109 | 1.3546 | 1.35914 | 0.00288774 | False |
| asymmetric_phosphine/reference/midpoint_128 | 0.2376 | 0.0011 | uncertified | 0.0111497 | None |
| asymmetric_phosphine/reference/midpoint_256 | 0.9564 | 0.0077 | uncertified | 0.00134418 | None |
| asymmetric_phosphine/reference/midpoint_512 | 3.6260 | 0.0125 | uncertified | 0.00135685 | None |
| fluorophosphine/reference/certified_0.1 | 8.6925 | 0.0203 | 0.0999909 | 0.0263012 | True |
| fluorophosphine/reference/certified_0.05 | 15.7896 | 0.0143 | 0.0499964 | 0.0147163 | True |
| fluorophosphine/reference/certified_0.01 | 58.8021 | 0.0564 | 0.00999983 | 0.00403294 | True |
| fluorophosphine/reference/reference_262143 | 54.2248 | 0.1550 | 1.48219 | 0.0096823 | False |
| fluorophosphine/reference/reference_1048575 | 227.1508 | 0.3896 | 0.742112 | 0.00349404 | False |
| fluorophosphine/reference/midpoint_128 | 0.1681 | 0.0015 | uncertified | 0.00258205 | None |
| fluorophosphine/reference/midpoint_256 | 0.6800 | 0.0086 | uncertified | 0.00132892 | None |
| fluorophosphine/reference/midpoint_512 | 2.4846 | 0.0171 | uncertified | 0.0014614 | None |
| methylamine/reference/certified_0.1 | 6.0421 | 0.0092 | 0.0999935 | 0.0249223 | True |
| methylamine/reference/certified_0.05 | 11.1520 | 0.0169 | 0.0499952 | 0.013496 | True |
| methylamine/reference/certified_0.01 | 46.7206 | 0.0348 | 0.00999987 | 0.00391434 | True |
| methylamine/reference/reference_262143 | 42.6481 | 0.0652 | 1.17619 | 0.00670587 | False |
| methylamine/reference/reference_1048575 | 181.9150 | 0.2898 | 0.587962 | 0.00283038 | False |
| methylamine/reference/midpoint_128 | 0.0986 | 0.0001 | uncertified | 0.00468723 | None |
| methylamine/reference/midpoint_256 | 0.3832 | 0.0028 | uncertified | 0.00121957 | None |
| methylamine/reference/midpoint_512 | 1.4851 | 0.0035 | uncertified | 0.00168881 | None |
| dimethylamine/reference/certified_0.1 | 13.6030 | 0.0186 | 0.0999962 | 0.0194366 | True |
| dimethylamine/reference/certified_0.05 | 24.5472 | 0.0255 | 0.0499984 | 0.0111648 | True |
| dimethylamine/reference/certified_0.01 | 98.3331 | 0.0687 | 0.00999991 | 0.00343291 | True |
| dimethylamine/reference/reference_262143 | 48.7177 | 0.0867 | 1.56708 | 0.0101437 | False |
| dimethylamine/reference/reference_1048575 | 204.5326 | 0.3789 | 0.785251 | 0.00391588 | False |
| dimethylamine/reference/midpoint_128 | 0.1165 | 0.0007 | uncertified | 0.00102563 | None |
| dimethylamine/reference/midpoint_256 | 0.4565 | 0.0091 | uncertified | 0.00138769 | None |
| dimethylamine/reference/midpoint_512 | 1.7843 | 0.0153 | uncertified | 0.00120718 | None |
| trimethylamine/reference/certified_0.1 | 12.5802 | 0.0256 | 0.099994 | 0.0236086 | True |
| trimethylamine/reference/certified_0.05 | 22.3708 | 0.0928 | 0.0499967 | 0.0134362 | True |
| trimethylamine/reference/certified_0.01 | 87.2184 | 0.2730 | 0.0099999 | 0.0037877 | True |
| trimethylamine/reference/reference_262143 | 54.5960 | 0.0480 | 1.83818 | 0.0068573 | False |
| trimethylamine/reference/reference_1048575 | 226.4337 | 0.7843 | 0.921748 | 0.00220876 | False |
| trimethylamine/reference/midpoint_128 | 0.1592 | 0.0009 | uncertified | 0.00837883 | None |
| trimethylamine/reference/midpoint_256 | 0.6463 | 0.0186 | uncertified | 0.00178127 | None |
| trimethylamine/reference/midpoint_512 | 2.4033 | 0.0038 | uncertified | 0.00114048 | None |
| aniline/reference/certified_0.1 | 16.4574 | 0.0484 | 0.0999929 | 0.0194525 | True |
| aniline/reference/certified_0.05 | 28.7739 | 0.0563 | 0.0499962 | 0.0107798 | True |
| aniline/reference/certified_0.01 | 109.3017 | 0.1446 | 0.00999996 | 0.00336561 | True |
| aniline/reference/reference_262143 | 71.3956 | 0.2673 | 1.58963 | 0.00296663 | False |
| aniline/reference/reference_1048575 | 294.2348 | 1.8953 | 0.796235 | 0.00143586 | False |
| aniline/reference/midpoint_128 | 0.2182 | 0.0005 | uncertified | 0.00115463 | None |
| aniline/reference/midpoint_256 | 0.8644 | 0.0068 | uncertified | 0.00116342 | None |
| aniline/reference/midpoint_512 | 3.2238 | 0.0216 | uncertified | 0.00154834 | None |
| coaxial_0/axis_0/certified_0.1 | 2.6344 | 0.0097 | 0.0999797 | 0.0283776 | True |
| coaxial_0/axis_0/certified_0.05 | 4.9270 | 0.0115 | 0.0499882 | 0.0149009 | True |
| coaxial_0/axis_0/certified_0.01 | 22.5542 | 0.0487 | 0.0099999 | 0.00375714 | True |
| coaxial_0/axis_0/reference_262143 | 48.9909 | 0.1084 | 0.456735 | 0.00706866 | False |
| coaxial_0/axis_0/reference_1048575 | 206.6348 | 0.1897 | 0.228136 | 0.00152729 | False |
| coaxial_0/axis_0/midpoint_128 | 0.0971 | 0.0002 | uncertified | 0.00157745 | None |
| coaxial_0/axis_0/midpoint_256 | 0.3703 | 0.0015 | uncertified | 5.04752e-05 | None |
| coaxial_0/axis_0/midpoint_512 | 1.4854 | 0.0119 | uncertified | 0.000133271 | None |
| coaxial_0/axis_1/certified_0.1 | 2.5744 | 0.0055 | 0.0999906 | 0.0285306 | True |
| coaxial_0/axis_1/certified_0.05 | 4.8697 | 0.0185 | 0.0499909 | 0.0151132 | True |
| coaxial_0/axis_1/certified_0.01 | 22.1543 | 0.0529 | 0.00999997 | 0.0038064 | True |
| coaxial_0/axis_1/reference_262143 | 49.0320 | 0.0636 | 0.459071 | 0.0073108 | False |
| coaxial_0/axis_1/reference_1048575 | 206.3580 | 0.6916 | 0.229544 | 0.00191845 | False |
| coaxial_0/axis_1/midpoint_128 | 0.0966 | 0.0004 | uncertified | 0.00157745 | None |
| coaxial_0/axis_1/midpoint_256 | 0.3709 | 0.0011 | uncertified | 5.04752e-05 | None |
| coaxial_0/axis_1/midpoint_512 | 1.5007 | 0.0338 | uncertified | 0.000133271 | None |
| coaxial_0/axis_2/certified_0.1 | 2.3909 | 0.0141 | 0.0999606 | 0.0208634 | True |
| coaxial_0/axis_2/certified_0.05 | 4.4151 | 0.0126 | 0.0499915 | 0.0128345 | True |
| coaxial_0/axis_2/certified_0.01 | 18.7527 | 0.0554 | 0.00999945 | 0.00321525 | True |
| coaxial_0/axis_2/reference_262143 | 49.0064 | 0.0888 | 0.459235 | 0.00735551 | False |
| coaxial_0/axis_2/reference_1048575 | 206.2222 | 0.0978 | 0.229695 | 0.00193242 | False |
| coaxial_0/axis_2/midpoint_128 | 0.0967 | 0.0003 | uncertified | 0.00371181 | None |
| coaxial_0/axis_2/midpoint_256 | 0.3791 | 0.0023 | uncertified | 4.2411e-05 | None |
| coaxial_0/axis_2/midpoint_512 | 1.4886 | 0.0065 | uncertified | 0.000122119 | None |
| coaxial_8/axis_0/certified_0.1 | 19.7142 | 0.0775 | 0.0999998 | 0.0215553 | True |
| coaxial_8/axis_0/certified_0.05 | 36.0138 | 0.0854 | 0.049999 | 0.0123272 | True |
| coaxial_8/axis_0/certified_0.01 | 151.3802 | 0.2215 | 0.00999999 | 0.00327981 | True |
| coaxial_8/axis_0/reference_262143 | 96.8119 | 0.1959 | 1.94831 | 0.004604 | False |
| coaxial_8/axis_0/reference_1048575 | 394.2325 | 1.4899 | 0.9769 | 0.000244232 | False |
| coaxial_8/axis_0/midpoint_128 | 0.2321 | 0.0028 | uncertified | 0.00314195 | None |
| coaxial_8/axis_0/midpoint_256 | 0.9276 | 0.0136 | uncertified | 0.0005166 | None |
| coaxial_8/axis_0/midpoint_512 | 3.5692 | 0.0088 | uncertified | 0.000808946 | None |
| coaxial_18/axis_0/certified_0.1 | 5.3755 | 0.0072 | 0.0999891 | 0.0225678 | True |
| coaxial_18/axis_0/certified_0.05 | 10.0936 | 0.0597 | 0.0499959 | 0.0127497 | True |
| coaxial_18/axis_0/certified_0.01 | 42.3416 | 0.0717 | 0.00999998 | 0.00326768 | True |
| coaxial_18/axis_0/reference_262143 | 67.2215 | 0.1785 | 0.721477 | 0.00273531 | False |
| coaxial_18/axis_0/reference_1048575 | 278.0432 | 0.6465 | 0.361758 | 0.000948107 | False |
| coaxial_18/axis_0/midpoint_128 | 0.1654 | 0.0009 | uncertified | 0.00191347 | None |
| coaxial_18/axis_0/midpoint_256 | 0.6588 | 0.0233 | uncertified | 0.000940777 | None |
| coaxial_18/axis_0/midpoint_512 | 2.5417 | 0.0028 | uncertified | 0.000301299 | None |
| coaxial_35/axis_0/certified_0.1 | 12.3693 | 0.0146 | 0.0999935 | 0.0225826 | True |
| coaxial_35/axis_0/certified_0.05 | 22.4677 | 0.0623 | 0.049995 | 0.0126998 | True |
| coaxial_35/axis_0/certified_0.01 | 88.6346 | 0.1766 | 0.00999999 | 0.00307973 | True |
| coaxial_35/axis_0/reference_262143 | 80.2873 | 0.2778 | 1.41487 | 0.0127974 | False |
| coaxial_35/axis_0/reference_1048575 | 329.7039 | 2.1723 | 0.707216 | 0.0020424 | False |
| coaxial_35/axis_0/midpoint_128 | 0.2083 | 0.0030 | uncertified | 0.00534407 | None |
| coaxial_35/axis_0/midpoint_256 | 0.8387 | 0.0078 | uncertified | 0.0013321 | None |
| coaxial_35/axis_0/midpoint_512 | 3.1628 | 0.0122 | uncertified | 0.00010099 | None |
