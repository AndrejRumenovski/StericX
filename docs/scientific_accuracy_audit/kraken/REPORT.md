# Independent Kraken reproduction audit

Audited system: commit `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7`. Production scientific code was not changed. Raw native outputs were frozen before independent comparison. This is an audit of numerical definitions and reproduction, not experimental validation.

## Main findings

Fresh acquisition covers all 1,566 IDs in the historical comparison universe: 1,546 have 31,721 available DFT conformers; 20 have no available DFT geometries. The audit attempted every available conformer, including cases excluded by the historical geometric prefilter. Complete native ensembles exist for 1,543 buried-volume ligands, 1,544 coordination-Sterimol ligands and all 1,546 pyramidalization ligands. No numerical outlier was removed. Eight conformers produce false fourth-neighbor inferences relative to independently reconstructed SDF topology; another valid symmetric molecule produces a buried-volume error. These failures are retained separately, never replaced with partial-ensemble minima.

The headline minimum of maximum adjacent-quadrant differences reproduces the approximate historical aggregate score: N=1,543, MAE=0.270886692 Å³, RMSE=0.490777865 Å³, one-to-one R²=0.985128630, maximum error=4.559474133 Å³. This does not establish exact convention equivalence. Every one of the twenty largest headline discrepancies becomes ≤ 0.006275429 Å³ against the published values when independently evaluated with the original Kraken center and density. The original DFT center is a normalized sum of **raw bond displacement vectors**, not an electronic localized-orbital center. StericX normalizes each bond vector before summing. Tertiary donors also differ; the two centers do not generally coincide.

Sterimol comparisons additionally differ because original Kraken substitutes H radius 1.09 Å while StericX uses 1.20 Å, and the angular scans differ. Removing those convention differences does not eliminate every implementation discrepancy: all-corpus same-center/same-radii closed-form projections reveal native L error 0.0036931285331647246 Å and B5 error 0.004176904379527002 Å in near-axis real conformers. The analytic formulas avoid the approximate quaternion rotation branch.

Pyramidalization on identical exported geometries agrees closely with independent contemporary Morfeus over all 31,721 conformers: maximum |ΔP| = 3.2010947048632943e−7 and |Δalpha| = 4.8856123654239525e−5 degrees. Yet the published ligand extrema differ by up to 0.007218747 in P and 1.155986692 degrees in alpha. Independent Morfeus reproduces those large discrepancies too. Therefore they cannot honestly be attributed solely to Rust f32 rounding; historical geometry/ensemble/reference provenance remains unresolved.

## Sources and exact conventions

- Gensch et al., JACS 144,1205–1217 (2022), [10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718), original [SI PDF](https://acs.figshare.com/articles/journal_contribution/Supplementary_material_for_A_Comprehensive_Discovery_Platform_for_Organophosphorus_Ligands_for_Catalysis_/18307112) and [SI archive](https://acs.figshare.com/articles/journal_contribution/Supplementary_material_for_A_Comprehensive_Discovery_Platform_for_Organophosphorus_Ligands_for_Catalysis_/18307121). Exact downloaded bytes, source URLs/status/timestamps/headers and hashes are in `raw_sources/`. PDF pages 8,10,12 are rendered there for equation/convention inspection.
- Original [DFT implementation at immutable commit 4eaad505c1343e6083032b4a3fda47e004e19734](https://github.com/the-matter-lab/kraken/blob/4eaad505c1343e6083032b4a3fda47e004e19734/conf_selection_and_DFT/PL_dft_library_201027.py), functions `get_conmat`, `add_valence`, `morfeus_properties`, `read_ligand`. This is evidence of the reference implementation, not physical truth.
- Fresh [MolSSI Kraken API](https://descriptor-libraries.molssi.org/api/kraken/docs) molecule metadata, published DFT descriptors and SDF conformer exports. Each response is frozen with provenance; failed responses remain preserved.
- Morfeus 0.8.0, frozen installed source and package versions in the master initial manifest. Contemporary Morfeus is an independent implementation; it is not claimed to be the exact 2021 dependency version.

| Convention | Original Kraken DFT | Frozen StericX measurement |
|---|---|---|
| Donor/geometry | First eligible trivalent P by primary connectivity; original graph rule frozen | Same mapped P chosen from independent SDF topology; SUT separately infers its neighbors |
| Virtual center | P +2.28 normalize(Σ(P − neighbor)) | P +2.28 normalize(Σ normalize(P − neighbor)) |
| Vbur | Sphere 3.5 Å; Bondi ×1.17; H excluded; density 0.001 | Same sphere/radii/H; native default density 0.01 |
| Frame | Three substituent planes; original neighbor order; extrema across orientations | Reference first then detected neighbors; extrema across orientations |
| Total/near/far | Last orientation retained | First orientation retained |
| Sterimol | Dummy Pd→P; Bondi except 1.20→1.09; 3,600 directions; L +0.40 Å | Dummy→P; H 1.20; 360 directions; L +0.40 Å |
| Pyramidalization | Three nearest neighbors excluding dummy Pd; P dimensionless, alpha degrees | Explicit SDF-bonded three neighbors in observation adapter; automatic donor/frame audit is reported separately |
| Ensemble | min/max/delta and descriptor at minimum-Vbur conformer; Boltzmann reduction uses source free energies | Exact full available ensembles for four unweighted reductions; source weights unavailable |

SDF connectivity was reconstructed by the API export tool (Open Babel), so it is independent of StericX but still an inference, not an experimental bond determination. Original Kraken and SDF neighbors agree in the diagnostic records; current automatic StericX connectivity can disagree. The export coordinates have four decimal places. Index mapping is unchanged; dummy Pd is appended only in independent calculations and excluded from occupancy.

## Freeze, reproducibility and coverage

`primary/acquisition_plan.json`, `primary/ligands_manifest.json`, `primary/conformers_manifest.json`, `prepared_all/manifest.json`, `sut/manifest_frozen.json`, `morfeus_outliers/selection_frozen.json`, `morfeus_outliers/manifest_frozen.json`, `full_analytic_reference/manifest_frozen.json` and `interpretation/manifest.json` form the phase manifest chain. Native values include exact observed f32 values and IEEE bit patterns. Reference calculations retain full available Python float precision.

`analysis/comparisons.csv` contains molecule-by-molecule results for 56 descriptor/reduction combinations; `analysis/conformers.csv` contains every native conformer result. Complete property ensembles are required. `analysis/excluded_comparisons.json` lists every unavailable input or failed property; `analysis/operation_errors.json` contains all 17 operation errors on 9 conformers. `analysis/top20_by_descriptor.json` retains 20 worst entries for every reported metric. `delta_interpretation/all_top20_outlier_dossiers.json` joins each entry to independent diagnostics, exact input IDs and the scope/limitations of the causal interpretation.

The published SI workbook contains 1,558 rows; the newer Ni-hDA comparison table contains 1,566. The extra IDs are 724,1057,1058,2059,2062,2063,2064,2067. Across 190 shared numeric columns the largest SI/API-table difference is 1.000000082740371e−8. `primary_si/` preserves exact XML-cell extraction and every difference. The corpus expansion and independently retained failures explain why this audit reports different coverage from the historical 1,541/31,611 headline; it did not select IDs to force that count.

Reproduction scripts are under `scripts/`. Run acquisition/preparation/native observation before reference interpretation in a fresh audit destination; raw phase directories intentionally refuse overwriting. `kraken_run_sut.py` observes the frozen production implementation through the audited adapter. `kraken_morfeus_outliers.py` performs the untuned convention matrix; `kraken_full_analytic_reference.py` evaluates all 31,721 conformers; `kraken_interpret.py` performs independent dot-product L/B5 projections and statistics. Source/method changes require a new phase directory. The master reproducibility procedure documents complete commands and environment.

## Full numerical comparison

Each row below uses published values as x and native values as y. R² is **1−Σ(y−x)²/Σ(x−mean(x))²**, not squared correlation. Slope and intercept are diagnostic OLS y~x only; they are never applied to correct native descriptors. `analysis/metrics.json` preserves full precision. Absolute BV fields use Å³; Sterimol uses Å; P is dimensionless; alpha uses degrees; percent Vbur uses percentage points.

| Descriptor/reduction | N | MAE | RMSE | Max | Median | R² | Slope | Intercept |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| buried_volume_delta | 1543 | 0.340515831 | 0.634084844 | 9.7993614 | 0.153508182 | 0.999092477 | 0.997094506 | 0.00435254255 |
| buried_volume_max | 1543 | 0.438819604 | 0.763159272 | 9.85005516 | 0.222343617 | 0.99905338 | 0.992273634 | 0.613161419 |
| buried_volume_min | 1543 | 0.266436872 | 0.49792411 | 4.47994984 | 0.137182132 | 0.998165385 | 0.984409054 | 0.936689681 |
| buried_volume_vburminconf | 1543 | 0.266436872 | 0.49792411 | 4.47994984 | 0.137182132 | 0.998165385 | 0.984409054 | 0.936689681 |
| far_vbur_delta | 1543 | 0.313655915 | 0.649715706 | 9.74494485 | 0.0941483778 | 0.998133203 | 0.995527332 | -0.0218152547 |
| far_vbur_max | 1543 | 0.376659527 | 0.768070021 | 9.71412137 | 0.115855688 | 0.997671505 | 0.991311339 | -0.0263496396 |
| far_vbur_min | 1543 | 0.0896633414 | 0.369685196 | 4.64267888 | 0 | 0.994030729 | 0.95383944 | 0.00113641001 |
| far_vbur_vburminconf | 1543 | 0.100627329 | 0.391899524 | 4.64267888 | 0 | 0.993414569 | 0.954366071 | 0.00765325001 |
| max_delta_qvbur_delta | 1543 | 0.334647423 | 0.626877007 | 7.63063705 | 0.143850031 | 0.994243285 | 0.998621299 | 0.00925293148 |
| max_delta_qvbur_max | 1543 | 0.385909256 | 0.674079946 | 4.71705835 | 0.193699041 | 0.994145734 | 0.991200641 | 0.122498077 |
| max_delta_qvbur_min | 1543 | 0.270886692 | 0.490777865 | 4.55947413 | 0.111850105 | 0.98512863 | 0.970133304 | 0.139005963 |
| max_delta_qvbur_vburminconf | 1543 | 0.347538351 | 0.646200613 | 6.31337608 | 0.148191936 | 0.977900052 | 0.966562003 | 0.184685305 |
| near_vbur_delta | 1543 | 0.176115256 | 0.337257441 | 4.18997284 | 0.096162624 | 0.998474147 | 0.99771317 | 0.0482300608 |
| near_vbur_max | 1543 | 0.247061038 | 0.486672951 | 4.6206651 | 0.141613417 | 0.998197652 | 0.993924194 | 0.563182971 |
| near_vbur_min | 1543 | 0.225601961 | 0.416966369 | 4.47994221 | 0.132283874 | 0.997423934 | 0.995873936 | 0.376892341 |
| near_vbur_vburminconf | 1543 | 0.236297933 | 0.449815722 | 4.47994221 | 0.133402482 | 0.997061871 | 0.995534747 | 0.390473126 |
| ovbur_max_delta | 1543 | 0.178084391 | 0.352457418 | 3.9858847 | 0.0740937876 | 0.985254306 | 0.993855196 | 0.0364830469 |
| ovbur_max_max | 1543 | 0.161536661 | 0.373698596 | 4.25148796 | 0.0549356542 | 0.987068195 | 0.975054238 | 0.531983966 |
| ovbur_max_min | 1543 | 0.187427135 | 0.353675956 | 2.78333109 | 0.0810804913 | 0.984459016 | 0.977109045 | 0.397688312 |
| ovbur_max_vburminconf | 1543 | 0.224489949 | 0.404354127 | 2.97263984 | 0.103337057 | 0.980790153 | 0.976512168 | 0.415024524 |
| ovbur_min_delta | 1543 | 0.0244564498 | 0.106561877 | 2.06233883 | 0 | 0.996715174 | 1.00372992 | -0.0051510329 |
| ovbur_min_max | 1543 | 0.0243389265 | 0.103592129 | 2.06233883 | 0 | 0.997082292 | 1.00287044 | -0.0056115878 |
| ovbur_min_min | 1543 | 0.000916890209 | 0.0167609275 | 0.614029811 | 0 | 0.998185317 | 0.988651735 | -0.000562952265 |
| ovbur_min_vburminconf | 1543 | 0.00210816764 | 0.0486000576 | 1.81577964 | 0 | 0.987703666 | 1.00237867 | 0.000435056806 |
| percent_buried_volume_delta | 1543 | 0.189602669 | 0.353065021 | 5.456388 | 0.0854695595 | 0.999092477 | 0.997094534 | 0.00242355611 |
| percent_buried_volume_max | 1543 | 0.244339455 | 0.42493495 | 5.48461464 | 0.123799699 | 0.99905338 | 0.992273666 | 0.341414366 |
| percent_buried_volume_min | 1543 | 0.148355162 | 0.277249412 | 2.4944832 | 0.0763876636 | 0.998165383 | 0.984409076 | 0.521558537 |
| percent_buried_volume_vburminconf | 1543 | 0.148355162 | 0.277249412 | 2.4944832 | 0.0763876636 | 0.998165383 | 0.984409076 | 0.521558537 |
| pyr_alpha_delta | 1546 | 0.00708419241 | 0.0381641803 | 1.15884495 | 0.00284351171 | 0.999906291 | 0.998977969 | 0.0113543326 |
| pyr_alpha_max | 1546 | 0.003777242 | 0.0304185768 | 1.15598669 | 0.00175196729 | 0.999967911 | 0.999727664 | 0.00841305235 |
| pyr_alpha_min | 1546 | 0.00445035524 | 0.0224271183 | 0.719327362 | 0.00191772205 | 0.999979269 | 0.999733974 | 0.000908332615 |
| pyr_alpha_vburminconf | 1543 | 0.0900094593 | 0.627665909 | 8.87123298 | 0.0019781169 | 0.985429689 | 0.994696462 | 0.125158377 |
| pyr_p_delta | 1546 | 5.46318007e-05 | 0.000268356805 | 0.00723212957 | 2.04179883e-05 | 0.999948189 | 0.999367538 | 7.16451604e-05 |
| pyr_p_max | 1546 | 3.20701927e-05 | 0.000170942074 | 0.00481712251 | 1.19241416e-05 | 0.999977619 | 0.99955338 | 0.000445611415 |
| pyr_p_min | 1546 | 3.07904868e-05 | 0.000201114618 | 0.00721874721 | 1.3845446e-05 | 0.999983938 | 0.999913263 | 5.72160861e-05 |
| pyr_p_vburminconf | 1543 | 0.000746124863 | 0.00535039001 | 0.0844949813 | 1.46332716e-05 | 0.983025693 | 0.996907267 | 0.00263421171 |
| qvbur_max_delta | 1543 | 0.30252544 | 0.588447535 | 6.60837613 | 0.119328839 | 0.995949632 | 0.996745844 | 0.0440367488 |
| qvbur_max_max | 1543 | 0.358552912 | 0.638657758 | 4.27092068 | 0.176774669 | 0.996235526 | 0.990211219 | 0.272153349 |
| qvbur_max_min | 1543 | 0.243719391 | 0.451200573 | 4.27092068 | 0.104228558 | 0.992540046 | 0.970356482 | 0.495292713 |
| qvbur_max_vburminconf | 1543 | 0.290958052 | 0.529261635 | 5.28085605 | 0.134297688 | 0.990241923 | 0.970379181 | 0.507602251 |
| qvbur_min_delta | 1543 | 0.123950799 | 0.263107135 | 3.77152352 | 0.0553473957 | 0.995536753 | 1.00182228 | 0.0371110855 |
| qvbur_min_max | 1543 | 0.134859023 | 0.253092755 | 3.59669799 | 0.0660947265 | 0.996223744 | 1.00602643 | -0.0330350281 |
| qvbur_min_min | 1543 | 0.0890655601 | 0.16743433 | 1.97809971 | 0.0444897659 | 0.990015371 | 1.00811381 | -0.0800297477 |
| qvbur_min_vburminconf | 1543 | 0.126650081 | 0.251594208 | 3.10782459 | 0.0595721843 | 0.979802171 | 1.00710805 | -0.0658026993 |
| sterimol_b1_delta | 1544 | 0.0330976834 | 0.0555472663 | 0.499836171 | 0.0189847116 | 0.99317081 | 1.00469067 | 0.0112949301 |
| sterimol_b1_max | 1544 | 0.106384969 | 0.11392858 | 0.366607011 | 0.114507216 | 0.982570634 | 1.01636646 | 0.0331605811 |
| sterimol_b1_min | 1544 | 0.0963253418 | 0.105101631 | 0.437885613 | 0.10921397 | 0.981465952 | 1.01113338 | 0.0490739898 |
| sterimol_b1_vburminconf | 1543 | 0.117035402 | 0.149929192 | 1.06658566 | 0.113947736 | 0.963334579 | 1.01164988 | 0.0490038173 |
| sterimol_b5_delta | 1544 | 0.0305668627 | 0.0607627963 | 0.648117579 | 0.00982002275 | 0.997305909 | 0.996748568 | 0.0107205604 |
| sterimol_b5_max | 1544 | 0.104618244 | 0.115741522 | 0.399393848 | 0.109653099 | 0.995548325 | 0.996413074 | 0.130244146 |
| sterimol_b5_min | 1544 | 0.102531159 | 0.118760709 | 0.777098718 | 0.107344321 | 0.992708097 | 0.992719233 | 0.144040347 |
| sterimol_b5_vburminconf | 1543 | 0.123430171 | 0.186987617 | 2.73796523 | 0.108833812 | 0.985116635 | 0.990359449 | 0.170902256 |
| sterimol_l_delta | 1544 | 0.0575956319 | 0.100624442 | 0.952176094 | 0.0275007033 | 0.996433735 | 0.998932579 | 0.0115945051 |
| sterimol_l_max | 1544 | 0.112372573 | 0.135827639 | 0.986380092 | 0.109156596 | 0.993568458 | 1.00502317 | 0.0516730224 |
| sterimol_l_min | 1544 | 0.110917866 | 0.132462238 | 0.759902644 | 0.109507274 | 0.986406818 | 1.00568952 | 0.0460905711 |
| sterimol_l_vburminconf | 1543 | 0.143046385 | 0.250159927 | 3.0846052 | 0.109773015 | 0.971915021 | 0.998618885 | 0.100949928 |

Residual figures: [buried-volume minima](analysis/buried_volume_min_residuals.png), [Sterimol extrema](analysis/sterimol_extrema_residuals.png), [pyramidalization extrema](analysis/pyramidalization_extrema_residuals.png). The Sterimol plots show the systematic positive radius-convention shift clearly; a large R² does not remove that bias.

## Direct Morfeus and convention isolation

The 590-record selection was frozen before calculating the reference variants. It includes all conformers of the 20 largest headline outliers, every available SUT extremizer for all per-descriptor top 20 min/max/vburminconf entries, and every native operation failure. This is a deliberately adverse diagnostic selection, not a random estimate of overall error. Native-success denominators differ because errors are retained, not imputed.

An initial matched-radii adapter accidentally omitted explicit 1.17 scaling because Morfeus does not rescale supplied radii. That failed campaign is preserved in `morfeus_outliers_attempt1/INVALID_MATCHED_VARIANT.json`; it is **not scientific evidence against StericX**. The corrected campaign explicitly scales the supplied radii and was rerun before interpretation. An initial analysis inventory parser also failed before comparisons; its failure record is preserved in `interpretation_attempt1/`.

With center, radii, coordinates, density and first-orientation convention matched, direct Morfeus total BV differs by at most 1.095939546758018e−5 Å³ over 581 native-success cases. The corresponding maximum percent Vbur difference is 8.34638332491977e−6 percentage points, including native percent rounding (`interpretation/matched_default_percent.json`). Quadrant/octant extrema contain one grid-cell discrepancy: max 0.011656810711858867 Å³ in adjacent-quadrant difference. Thus the total-volume kernel shows strong agreement on those cases while finite boundary/frame arithmetic remains measurable. Matching default density is essential: changing only the independent density to 0.001 gives total BV difference up to 0.45010229210903674 Å³. That is convergence sensitivity, not an exactness failure by itself.

For identical coordinates/center/radii, Morfeus B1 differs by up to 0.044101899989216875 Å on 582 diagnostic cases. Its 3,600-direction scan and the native 360-direction scan do not define identical numerical approximations. L/B5 maxima on that subset are 0.0017224155361290627/0.002769704392543204 Å. The independent all-corpus closed-form campaign broadens the observation to 31,713 valid native axes and finds larger deviations:

| Descriptor | Max error Å | Native | Independent | Exact case |
|---|---:|---:|---:|---|
| sterimol_b5 | 0.0041769043795270022 | 8.3472118377685547 | 8.3513887421480817 | KRAKEN:1407:54385 |
| sterimol_l | 0.0036931285331647246 | 7.6446890830993652 | 7.6409959545662005 | KRAKEN:16:31923 |

The formulas are L=max_i[(x_i−c)·u+r_i]+0.40 and B5=max_i[||(x_i−c)−((x_i−c)·u)u||+r_i], where u=(x_P−c)/||x_P−c||. They require no quaternion and do not call a scientific StericX kernel. Frozen coordinates, radii and center are used verbatim, promoted to float64. Most errors are much smaller (median about 4–5e−7 Å), but the near-parallel/antiparallel rotation approximation in `glam::Quat::from_rotation_arc` is an implementation limitation that a mean would conceal. See `focused_diagnosis/axis_counterexamples.json` for direct Morfeus corroboration and complete failing inputs.

## Twenty largest headline discrepancies

All 20 complete ensembles were independently recalculated under both centers and both densities. The table gives the original published reference, native default result and untuned primary-convention Morfeus result. Every input, donor, neighbor list, exact center, four variant values and proposed correction is retained in `interpretation/headline_outlier_dossiers.json`.

| Kraken ID | P–H count | Published Å³ | Native Å³ | Error Å³ | Primary Morfeus Å³ | Remaining error Å³ |
|---|---:|---:|---:|---:|---:|---:|
| 1796 | 0 | 19.71216502 | 15.15269089 | -4.559474133 | 19.7079814 | -0.004183615617 |
| 1455 | 0 | 6.741901403 | 2.948947906 | -3.792953497 | 6.735625975 | -0.006275428143 |
| 1454 | 0 | 2.476702222 | 6.165979385 | 3.689277163 | 2.475656317 | -0.001045904844 |
| 1456 | 0 | 6.777462161 | 3.648302078 | -3.129160083 | 6.774324447 | -0.003137713928 |
| 254 | 0 | 5.745154267 | 2.622580528 | -3.122573739 | 5.746200172 | 0.001045904718 |
| 1712 | 0 | 3.029985784 | 5.932861328 | 2.902875544 | 3.031031689 | 0.001045904686 |
| 1654 | 0 | 5.668803227 | 2.867355347 | -2.80144788 | 5.666711418 | -0.002091809023 |
| 1064 | 0 | 22.67625881 | 20.14142609 | -2.534832724 | 22.67730472 | 0.001045908454 |
| 1675 | 0 | 4.023595206 | 1.608517647 | -2.415077559 | 4.025687015 | 0.002091809097 |
| 1488 | 2 | 6.231499931 | 8.333979607 | 2.102479676 | 6.231499931 | 3.969322648e-10 |
| 1486 | 2 | 4.851951692 | 6.888647079 | 2.036695387 | 4.851951692 | 1.366844415e-10 |
| 912 | 0 | 5.569442285 | 3.543399811 | -2.026042474 | 5.569442285 | 1.10547127e-10 |
| 1298 | 2 | 3.780945326 | 5.734710693 | 1.953765367 | 3.779899421 | -0.001045904706 |
| 424 | 0 | 9.6850771 | 7.75118351 | -1.93389359 | 9.686123005 | 0.001045905147 |
| 1449 | 0 | 5.042306339 | 6.94692421 | 1.904617871 | 5.042306339 | 2.522000386e-10 |
| 1296 | 2 | 3.779899421 | 5.664776802 | 1.884877381 | 3.778853517 | -0.001045904361 |
| 575 | 2 | 4.005814827 | 5.874582291 | 1.868767464 | 4.005814827 | -3.383258118e-10 |
| 1019 | 0 | 10.36177741 | 12.20374584 | 1.841968432 | 10.36177741 | 1.941822703e-09 |
| 1490 | 2 | 3.984896734 | 5.792991638 | 1.808094904 | 3.983850829 | -0.001045905082 |
| 1524 | 0 | 3.683676193 | 5.385033607 | 1.701357414 | 3.681584384 | -0.002091809227 |

Many of the largest cases contain no donor-bound hydrogens. For example ID 1796 has center directions separated by about 9.715–9.718 degrees and native headline error−4.559474133 Å³; the primary convention gives 19.707981404 Å³ versus published 19.712165020 Å³. Raw-versus-unit bond weighting also matters for heterogeneous substituent lengths and elements in tertiary donors. Labeling all such differences as an electronic P–H LMO effect is contradicted by the original DFT code and these calculations.

## Completed investigation of all delta outliers

All 280 top-20 delta entries across 14 descriptors now identify both native extremizers, exact input hashes, independent per-conformer values and complete reference-ensemble reductions. An additional frozen campaign covers every conformer of all 77 ligands among the buried-volume delta outliers: 2,009 conformers, with 218 exact records reused from the previous frozen campaign and 1,791 newly computed. The authoritative combined 1,120-entry dossier is `delta_interpretation/all_top20_outlier_dossiers.json`; the earlier derivative is preserved. No reported main-comparison input was removed.

Across the union of both diagnostic campaigns, matched default-grid direct Morfeus comparisons cover 2,372 native-success cases. Additional boundary-sensitive cases broaden the numerical limits beyond the earlier 581-case sample:

| Field | N | Maximum matched error | Exact case |
|---|---:|---:|---|
| buried_volume | 2372 | 0.011661367234964359 | KRAKEN:1769:61173 |
| far_vbur | 2372 | 7.1303596840266437e-06 | KRAKEN:51:32718 |
| max_delta_qvbur | 2372 | 0.011656810711858867 | KRAKEN:1675:58993 |
| near_vbur | 2372 | 0.011655667877072062 | KRAKEN:591:42889 |
| ovbur_max | 2372 | 0.011656053236984931 | KRAKEN:1675:58993 |
| ovbur_min | 2372 | 1.0496953315453084e-06 | KRAKEN:1233:52788 |
| percent_buried_volume | 2372 | 0.0064943656248033221 | KRAKEN:591:42889 |
| qvbur_max | 2372 | 0.011656053236984931 | KRAKEN:1675:58993 |
| qvbur_min | 2372 | 0.011655196130472945 | KRAKEN:591:42889 |

Absolute volumes use Å³; percent buried volume uses percentage points. The largest total and percent differences correspond to approximately one default integration-cell weight. These cases remain part of the quantified arithmetic/frame-boundary limitation, rather than being discarded as outliers.

Ligand 369 illustrates a different problem. Native total-volume delta is 76.17142486572266 Å³ versus published 66.37206347 Å³. Independent Morfeus with the native convention gives 76.17142221556755 Å³, whereas the original center/density gives 76.73697859460583 Å³. Thus its large published discrepancy is not explained by the native occupancy algorithm. All 84 exported geometries have the same canonical connectivity graph as the API SMILES when unspecified stereotags are excluded. SDF conversion adds stereotags to adamantyl cages that the API SMILES leaves unspecified. The API IDs form two contiguous blocks: 38367–38410 (44 conformers) and 64409–64448 (40). The latter contains conformer 64422, which sets the larger maximum. As a diagnostic only, the first block gives primary total delta 66.37310937335505 Å³ and far-hemisphere delta 56.1650799456219 Å³, within 0.001046 and 0.004184 Å³ of the published values. All 84 remain in the main comparison. This supports an ensemble mismatch; ID ordering alone does not prove insertion dates or the original reason for inclusion/exclusion. Exact geometry, graph, donor, center and block evidence is in `delta_interpretation/ligand_369_dossier.json`.

Other remaining primary-reference delta discrepancies, including ligands 1290 and 1805, remain explicitly UNCERTAIN in the dossiers. Four untuned complete-ensemble variants quantify center/grid effects separately. The primary `read_ligand` source removes calculation errors/imaginary-frequency structures and energy/RMSD duplicates before aggregation; per-property missing values can also reduce its extrema set. The current public exports lack the original selection/error/energy records. Our strict complete-available-geometry audit is therefore reproducible but cannot assert identity with every original retained calculation ensemble.

## Remaining published-value outliers

Primary-convention contemporary Morfeus was run on every available conformer for all three Sterimol values and both pyramidalization values. It substantially reduces overall Sterimol error, but residual outliers persist: L max for ID 1290 differs 0.9550323238 Å, B5 min for ID 1805 differs 0.6492286866 Å, B1 max for ID 560 differs 0.2742709684 Å. Pyramidalization alpha max for ID 1290 differs 1.156001260 degrees. These are independent disagreements with the published aggregate using currently exported geometries. They are not evidence that native code alone caused the published mismatch.

`interpretation/primary_ensemble_vs_published.csv` and `primary_ensemble_metrics.json` retain every complete reference ensemble without trimming. Independent nearest-neighbor and SDF-neighbor pyramidalization agree throughout this corpus; native differences on the same geometries are orders of magnitude below the largest published discrepancy. `focused_diagnosis/pyramidalization_coordinate_rounding.json` examines all 4,096 corners of the donor/three-neighbor ±0.00005 Å rounding box in major outlier geometries. Across 251 conformers and 1,028,096 corner evaluations, the largest alpha deviation is 0.0103514124 degrees; for ID 1290 it is at most 0.0089272857 degrees, far below its 1.156-degree published discrepancy. It is a sensitivity experiment, not a rigorous global bound over the box. A separate conservative determinant bound (`focused_diagnosis/pyramidalization_rounding_bound.json`) covers arbitrary continuous coordinate perturbations in that box. All 251 cases have positive-alpha branch and noncollinearity margins. The bound is at most 0.000623445508668 in P, including both ID 1290 conformers, so its 0.007218735 published P discrepancy cannot arise from ordinary 4-decimal rounding alone if ensemble and atom identity are unchanged. The defensible unresolved explanations include exported-versus-original geometry/ensemble provenance and historical reference dependency behavior. Raw Gaussian logs, original full-precision conformer coordinates, per-conformer computed properties and exact historical dependency versions are needed to choose among them. Neither dataset is declared wrong solely because they disagree.

## Failures and chemical eligibility

Ligand 1281 conformers 54235/54241 and ligand 1907 conformers 63291–63296 receive a fourth geometrically close neighbor in StericX; independent SDF topology and original primary connectivity have three. Eight BV and eight coordination-Sterimol errors are preserved. Atom identity, donor distances, radius thresholds, independent neighbor lists and outputs are in `focused_diagnosis/operation_failure_dossiers.json`. The historical prefilter mirrored the same native cutoff, so removing those cases before validation would conceal this failure mode. Explicit topology is preferable when trustworthy; a distance cutoff alone cannot establish chemical bonding.

Ligand 1299 conformer 54318 has three bonded substituents and yields a valid symmetric zero adjacent-quadrant difference in an independent calculation. Native BV rejects it as a degenerate frame. Zero anisotropy is not by itself chemically impossible or mathematically undefined. The exact molecule and all failed/reference records are retained; no algorithm was changed.

## Boltzmann data availability and limits

Full-library Boltzmann reductions cannot be independently reconstructed from the acquired historical inputs. This is an evidence limitation, not a successful reproduction. `energy_source_search/REPORT.md` records checks of both current/v2 API schemas and endpoints, five diverse conformers, alternate exports, the original immutable repository and releases, the complete SI archive contents, and the updated Sigman repository. API conformer records contain `data:null`; neither geometry export nor the official SI tables supply source per-conformer free energies/weights. The updated repository contains 33 regenerated validation ensembles, which cannot be substituted for historical source weights. Published aggregate Boltzmann descriptors are outcomes, not legitimate inputs for reproducing themselves.

## Classification and proposed follow-up

- **Supported approximate reproduction:** the headline published-reference score and the listed subset of geometric families, with all errors and missing inputs quantified.
- **Verified with numerical limits on tested geometry:** direct same-convention BV kernel agreement (including measured cell-boundary differences) and native pyramidalization versus Morfeus over the full exported corpus.
- **Incorrect claims:** electronic-LMO explanation of original Kraken DFT center; assertion that raw and normalized bond sums generally coincide for tertiary phosphines; exact axis rotation for all native L/B5 cases; blanket chemical validity of current geometric donor-neighbor inference; rejection of all symmetric zero-anisotropy cases as invalid.
- **Terminology issue:** saying the “entire vbur family” was validated when only eight absolute geometric fields and selected ensemble reductions were tested. Distal/total molecular volumes, ratios and source-energy Boltzmann reductions remain outside that numerical claim.
- **Uncertain:** why a minority of published extrema disagree with both native and independent calculations on currently exported geometries; exact historical per-conformer thermochemistry; predictive/experimental usefulness of these descriptor agreements.

First correct the documentation and qualify convention-specific claims. A later scientific code change could address the demonstrated axis approximation, false bonded neighbors and symmetric-zero rejection only after its desired convention is explicit and separately reviewed. Reproducing original Kraken values requires deciding whether to adopt its center/radius/density conventions; that is a descriptor-definition choice, not a harmless performance fix. Obtain original calculation records to resolve remaining published-value provenance. Do not remove disagreeing molecules or infer experimental generalization from descriptor parity.
