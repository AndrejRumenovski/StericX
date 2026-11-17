# Final controlled measurement tables

Derived by `summarize_final.py` from final window 6 and the Sterimol supplement. Every sample, median, unscaled MAD and range is retained in [final_measurements.json](final_measurements.json). The figures below are ratios of medians. Eight alternating trials per side follow the declared warmup. No outliers were removed. Process startup uses warm file caches; there is no claim of controlled cache-cold disk performance.

## Whole CLI, milliseconds

| Workload / threads | Original median | Final median | Final MAD | Final min–max | Speedup |
|---|---|---|---|---|---|
| single_small/t1 | 1.834935 | 1.612943 | 0.076059 | 1.512454–1.710867 | 1.1376× |
| single_large/t1 | 1.889082 | 1.665176 | 0.078407 | 1.571585–1.954596 | 1.1345× |
| conformers_56/t1 | 28.634266 | 13.938790 | 0.063109 | 13.823302–16.025915 | 2.0543× |
| ensemble_sdf/t1 | 45.531233 | 24.030831 | 0.127605 | 23.650420–24.170679 | 1.8947× |
| descriptors_10000/t1 | 4807.103096 | 2220.595841 | 5.345872 | 2210.362082–2228.781827 | 2.1648× |
| parse_ensembles_1000/t1 | 187.211308 | 187.502736 | 0.504544 | 186.478730–195.825124 | 0.9984× |
| db_build_ensembles/t1 | 29.458711 | 14.648334 | 0.139117 | 14.503031–14.952766 | 2.0111× |
| single_small/t2 | 1.895253 | 1.642574 | 0.049228 | 1.576484–1.700438 | 1.1538× |
| single_large/t2 | 1.909962 | 1.709014 | 0.046538 | 1.656094–1.809762 | 1.1176× |
| conformers_56/t2 | 15.857983 | 8.414767 | 0.023841 | 8.154343–8.448206 | 1.8845× |
| ensemble_sdf/t2 | 45.829169 | 23.888423 | 0.046512 | 23.781687–23.938813 | 1.9185× |
| descriptors_10000/t2 | 2493.138564 | 1147.947787 | 5.058633 | 1134.312923–1158.477230 | 2.1718× |
| parse_ensembles_1000/t2 | 187.274872 | 187.311741 | 0.389963 | 186.540185–208.600677 | 0.9998× |
| db_build_ensembles/t2 | 16.406560 | 8.839296 | 0.021355 | 8.652680–8.945380 | 1.8561× |
| single_small/t4 | 2.014829 | 1.737672 | 0.009092 | 1.622460–1.771141 | 1.1595× |
| single_large/t4 | 2.052158 | 1.818765 | 0.045886 | 1.674238–1.900634 | 1.1283× |
| conformers_56/t4 | 9.250103 | 5.177214 | 0.061631 | 5.066978–5.414681 | 1.7867× |
| ensemble_sdf/t4 | 46.020092 | 24.154939 | 0.182077 | 23.873981–24.807746 | 1.9052× |
| descriptors_10000/t4 | 1288.302073 | 604.378175 | 2.175786 | 599.470337–796.988866 | 2.1316× |
| parse_ensembles_1000/t4 | 187.113203 | 187.024877 | 0.636327 | 186.206248–188.377009 | 1.0005× |
| db_build_ensembles/t4 | 9.867154 | 5.903942 | 0.056161 | 5.833920–6.089530 | 1.6713× |
| single_small/t6 | 2.042235 | 1.887981 | 0.059137 | 1.743008–5.831555 | 1.0817× |
| single_large/t6 | 2.053501 | 1.879053 | 0.040085 | 1.765389–1.967981 | 1.0928× |
| conformers_56/t6 | 7.034007 | 4.252861 | 0.106961 | 4.142319–5.485925 | 1.6539× |
| ensemble_sdf/t6 | 46.198352 | 23.983482 | 0.130992 | 23.784533–24.757432 | 1.9263× |
| descriptors_10000/t6 | 889.517735 | 420.455000 | 1.456112 | 418.278829–426.360314 | 2.1156× |
| parse_ensembles_1000/t6 | 187.139788 | 186.921698 | 0.400402 | 186.519638–189.959284 | 1.0012× |
| db_build_ensembles/t6 | 7.853582 | 4.703710 | 0.065974 | 4.589720–7.044646 | 1.6697× |
| single_small/t12 | 2.112563 | 1.867677 | 0.098174 | 1.737457–2.008997 | 1.1311× |
| single_large/t12 | 2.214023 | 2.008422 | 0.074240 | 1.743007–15.161849 | 1.1024× |
| conformers_56/t12 | 6.660173 | 4.496325 | 0.066966 | 4.162066–6.528366 | 1.4812× |
| ensemble_sdf/t12 | 46.161508 | 24.411952 | 0.307138 | 23.952288–25.140192 | 1.8909× |
| descriptors_10000/t12 | 717.286099 | 353.096826 | 1.822623 | 348.856583–371.159209 | 2.0314× |
| parse_ensembles_1000/t12 | 188.006342 | 187.071936 | 0.371509 | 186.597574–195.039908 | 1.0050× |
| db_build_ensembles/t12 | 6.894484 | 4.731732 | 0.040726 | 4.640656–4.919030 | 1.4571× |

## 10,000-file throughput and parallel efficiency

| Threads | Original files/s | Final files/s | Original efficiency | Final efficiency |
|---|---|---|---|---|
| 1 | 2,080.3 | 4,503.3 | 100.00% | 100.00% |
| 2 | 4,011.0 | 8,711.2 | 96.41% | 96.72% |
| 4 | 7,762.2 | 16,545.9 | 93.28% | 91.85% |
| 6 | 11,242.0 | 23,783.8 | 90.07% | 88.02% |
| 12 | 13,941.4 | 28,320.8 | 55.85% | 52.41% |

## Isolated paired kernels, microseconds

| Kernel | Original | Final | Final MAD | Speedup |
|---|---|---|---|---|
| volume | 486.487 | 166.608 | 0.501 | 2.9199× |
| regions | 486.701 | 166.797 | 0.606 | 2.9179× |
| ensemble | 2884.573 | 1126.669 | 2.285 | 2.5603× |
| radius_profile | 1650.179 | 763.852 | 2.430 | 2.1603× |
| sterimol | 12.678 | 12.478 | 0.007 | 1.0161× |

## Python paired costs, microseconds per workflow

| Mode / threads | Original | Final | Final MAD | Speedup |
|---|---|---|---|---|
| bridge/t1 | 3136.765 | 1412.374 | 1.400 | 2.2209× |
| public/t1 | 3830.302 | 2098.555 | 2.247 | 1.8252× |
| json_decode/t1 | 403.835 | 405.565 | 3.601 | 0.9957× |
| typed_construction/t1 | 157.039 | 156.465 | 0.460 | 1.0037× |
| numpy_export/t1 | 53.136 | 53.228 | 0.300 | 0.9983× |
| options_encode/t1 | 5.769 | 5.771 | 0.026 | 0.9997× |
| record_export/t1 | 39.299 | 37.484 | 0.239 | 1.0484× |
| per_file_56/t1 | 35670.195 | 20511.314 | 17.311 | 1.7390× |
| batch_56/t1 | 34757.288 | 19823.454 | 22.205 | 1.7533× |
| batch_56/t2 | 21805.592 | 14109.568 | 9.720 | 1.5454× |
| batch_56/t4 | 14617.738 | 10738.175 | 95.713 | 1.3613× |
| batch_56/t6 | 12625.038 | 10167.044 | 177.424 | 1.2418× |
| batch_56/t12 | 11684.755 | 9414.058 | 20.223 | 1.2412× |

## Final native / bridge / public layers, microseconds

| Layer | Median | MAD | Min–max |
|---|---|---|---|
| native | 1131.724 | 1.897 | 1129.627–1186.343 |
| bridge | 1414.317 | 1.296 | 1412.727–1423.800 |
| public | 2110.277 | 5.737 | 2103.840–2149.486 |

## Final descriptor/I/O kernel profile, microseconds

| Mode | Median | MAD | Min–max |
|---|---|---|---|
| sterimol | 12.525 | 0.016 | 12.509–13.459 |
| volume | 165.783 | 0.147 | 165.636–166.528 |
| regions | 166.364 | 0.317 | 166.047–166.835 |
| pyramidalization | 0.103 | 0.000 | 0.102–0.103 |
| center | 0.070 | 0.001 | 0.068–0.071 |
| snapshot | 0.781 | 0.002 | 0.759–0.784 |
| radius_profile | 764.534 | 2.381 | 758.270–766.916 |
| ensemble | 1131.649 | 9.634 | 1119.383–1141.284 |
| statistics | 9.071 | 0.011 | 9.056–9.083 |
| serialize | 39.239 | 0.174 | 39.064–39.544 |
| parse | 53.873 | 0.036 | 53.837–54.245 |
| density | 166.324 | 0.278 | 165.955–166.602 |
| density_profile | 508.583 | 0.746 | 507.837–509.811 |
| density_parse | 5759.436 | 5.797 | 5753.639–5799.541 |
| batch/t1 | 1138.006 | 0.605 | 1137.400–1142.325 |
| batch/t2 | 603.568 | 0.212 | 603.356–605.143 |
| batch/t4 | 406.823 | 7.695 | 399.128–504.361 |
| batch/t6 | 232.187 | 1.551 | 230.636–267.259 |
| batch/t12 | 355.163 | 4.671 | 350.492–365.784 |

## Matched volume workflow, milliseconds

Warm workflow: each independent trial contributes its mean internal iteration time (100 iterations for one input; three for 56; one for 1,000). All inner samples are also retained. DBSTEP is the common-grid adapter, not the default DBSTEP workflow.

| Size / workers / scope / tool | Median | MAD | Min–max |
|---|---|---|---|
| n1/t1/fresh_process/native | 1.2259 | 0.0511 | 1.0819–1.2774 |
| n1/t1/fresh_process/morfeus | 291.8688 | 0.8589 | 288.4834–293.5671 |
| n1/t1/fresh_process/dbstep | 290.4582 | 0.3772 | 289.7397–291.1169 |
| n1/t1/warm_workflow/native | 0.1615 | 0.0007 | 0.1606–0.1684 |
| n1/t1/warm_workflow/morfeus | 14.1175 | 0.0140 | 14.0942–14.1566 |
| n1/t1/warm_workflow/dbstep | 13.5528 | 0.0125 | 13.5325–13.5801 |
| n56/t1/fresh_process/native | 13.5754 | 0.0373 | 13.4105–13.6960 |
| n56/t1/fresh_process/morfeus | 1049.5759 | 23.8129 | 1021.8402–1083.2663 |
| n56/t1/fresh_process/dbstep | 940.5795 | 24.0179 | 905.5671–964.8147 |
| n56/t1/warm_workflow/native | 11.5144 | 0.0166 | 11.4817–11.5938 |
| n56/t1/warm_workflow/morfeus | 756.1068 | 1.4931 | 740.1875–810.0808 |
| n56/t1/warm_workflow/dbstep | 629.1963 | 2.6957 | 625.6336–689.2857 |
| n56/t2/fresh_process/native | 7.9779 | 0.0343 | 7.9331–8.2871 |
| n56/t2/fresh_process/morfeus | 1030.6608 | 4.8505 | 1016.6426–1046.6077 |
| n56/t2/fresh_process/dbstep | 984.6459 | 8.1650 | 974.0317–1015.3582 |
| n56/t2/warm_workflow/native | 6.0817 | 0.0102 | 6.0235–6.1420 |
| n56/t2/warm_workflow/morfeus | 406.9420 | 2.1214 | 401.6335–420.3832 |
| n56/t2/warm_workflow/dbstep | 370.0996 | 6.1371 | 357.6581–378.2125 |
| n56/t4/fresh_process/native | 5.2692 | 0.1071 | 5.0748–6.7056 |
| n56/t4/fresh_process/morfeus | 890.6372 | 5.6280 | 881.8331–900.0459 |
| n56/t4/fresh_process/dbstep | 866.0683 | 2.3560 | 862.6690–885.1828 |
| n56/t4/warm_workflow/native | 3.0887 | 0.0175 | 3.0447–3.1122 |
| n56/t4/warm_workflow/morfeus | 237.5659 | 1.9302 | 231.5205–242.8794 |
| n56/t4/warm_workflow/dbstep | 220.4014 | 1.6392 | 216.8368–229.9259 |
| n56/t6/fresh_process/native | 4.3306 | 0.1208 | 3.9028–4.9755 |
| n56/t6/fresh_process/morfeus | 853.5370 | 4.4515 | 845.2637–857.9946 |
| n56/t6/fresh_process/dbstep | 836.8913 | 5.3436 | 827.4856–847.0265 |
| n56/t6/warm_workflow/native | 2.1434 | 0.0076 | 2.1258–2.1701 |
| n56/t6/warm_workflow/morfeus | 166.7860 | 2.2796 | 163.2405–173.3498 |
| n56/t6/warm_workflow/dbstep | 155.5162 | 3.8654 | 150.3113–161.9000 |
| n56/t12/fresh_process/native | 4.3381 | 0.0784 | 4.1740–4.6296 |
| n56/t12/fresh_process/morfeus | 1053.1430 | 2.9844 | 1045.9879–1066.8440 |
| n56/t12/fresh_process/dbstep | 1044.1457 | 4.0414 | 1039.5813–1056.3689 |
| n56/t12/warm_workflow/native | 1.8913 | 0.0115 | 1.8508–1.9165 |
| n56/t12/warm_workflow/morfeus | 138.0317 | 0.5858 | 135.6644–139.1472 |
| n56/t12/warm_workflow/dbstep | 128.6745 | 0.9109 | 127.1927–130.1683 |
| n1000/t1/fresh_process/native | 218.6766 | 0.3851 | 218.0743–236.8873 |
| n1000/t1/fresh_process/morfeus | 14166.1826 | 161.4751 | 13686.9683–14472.2575 |
| n1000/t1/fresh_process/dbstep | 12621.6005 | 18.9738 | 12599.6708–12708.8851 |
| n1000/t1/warm_workflow/native | 204.6523 | 0.1877 | 204.1201–204.8996 |
| n1000/t1/warm_workflow/morfeus | 13931.3272 | 158.1027 | 13287.9644–14402.1510 |
| n1000/t1/warm_workflow/dbstep | 12344.9501 | 38.4001 | 12122.0951–12765.5244 |
| n1000/t6/fresh_process/native | 51.0622 | 0.6254 | 50.3902–56.9110 |
| n1000/t6/fresh_process/morfeus | 3482.0191 | 38.2858 | 3382.3822–3547.2528 |
| n1000/t6/fresh_process/dbstep | 3333.3873 | 40.5840 | 3254.8501–3387.4356 |
| n1000/t6/warm_workflow/native | 36.5210 | 0.0673 | 36.4376–41.6897 |
| n1000/t6/warm_workflow/morfeus | 2769.4442 | 10.6314 | 2724.8751–2807.6166 |
| n1000/t6/warm_workflow/dbstep | 2570.0559 | 19.5743 | 2499.3315–2625.2753 |

## External 56-conformer scaling efficiency

| Scope / tool | 1 workers | 2 workers | 4 workers | 6 workers | 12 workers |
|---|---|---|---|---|---|
| fresh_process/native | 100.00% | 85.08% | 64.41% | 52.25% | 26.08% |
| fresh_process/morfeus | 100.00% | 50.92% | 29.46% | 20.49% | 8.31% |
| fresh_process/dbstep | 100.00% | 47.76% | 27.15% | 18.73% | 7.51% |
| warm_workflow/native | 100.00% | 94.66% | 93.20% | 89.53% | 50.73% |
| warm_workflow/morfeus | 100.00% | 92.90% | 79.57% | 75.56% | 45.65% |
| warm_workflow/dbstep | 100.00% | 85.00% | 71.37% | 67.43% | 40.75% |

## Sterimol cost and error, twelve-molecule in-memory batch

The fixed target is 0.01 Å on every L/B1/B5. Morfeus and DBSTEP each miss one B1 target; their timings are cost/error observations, excluded from equal-accuracy speed ratios. Setup/containers/imports are outside the timer; alignment and result construction are inside.

| Method | Batch median µs | MAD µs | Max residual Å | Target passed |
|---|---|---|---|---|
| sampled361 | 121.073 | 0.156 | 0.00880934454 | True |
| continuous | 167.353 | 0.167 | 2.92530729e-12 | True |
| morfeus | 1890.307 | 9.887 | 0.0125479246 | False |
| dbstep | 25030.612 | 68.302 | 0.0155438723 | False |
