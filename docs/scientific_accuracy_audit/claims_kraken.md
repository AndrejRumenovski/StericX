# Independent Kraken claim inventory

These are final classifications for the frozen audit evidence, not universal correctness labels. Each claim is scoped separately. Numerical reproduction never proves chemical predictive validity. Full data, negative findings and limitations: [Kraken report](kraken/REPORT.md).

## K01 — Full-library reproduction coverage

1. **StericX claim:** Study 004/README: 1,541 validated ligands and 31,611 DFT conformers from 1,566 reference IDs.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Every available conformer is attempted; failed/native-incomplete ensembles are reported, never reduced from successful subsets.
4. **Implementation location:** `studies/study_004_scaled.py; data/ligand_db/kraken_phosphines.manifest.json`
5. **Independent validation:** Fresh acquisition and independent SDF topology inventory, full native observation and complete-ensemble counting.
6. **Result — SUPPORTED:** Current public snapshot provides 1,546 ligands/31,721 conformers; 20 IDs lack DFT geometries. Full successful ensembles: BV 1,543, Sterimol 1,544, pyramidalization 1,546. Historical exact coverage count is a historical artifact, not the current independently eligible universe.
7. **Confidence:** High for frozen coverage.
8. **Limitations:** SDF connectivity comes from independent OpenBabel inference; available public snapshot can differ from original calculation archive.

## K02 — Headline buried-volume reproduction

1. **StericX claim:** Minimum of max adjacent-quadrant differences agrees with Kraken at approximately R² = 0.9852.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Published primary 3.5 Å sphere, 2.28 Å raw bond vector center, Bondi 1.17, H excluded, density 0.001, three orientations and complete-ensemble minimum.
4. **Implementation location:** `src/geometry/buried_volume.rs; studies/study_004_scaled.py`
5. **Independent validation:** Fresh native outputs against published values; all conventional error statistics, residuals, full 20 worst ensembles under four untuned Morfeus variants.
6. **Result — SUPPORTED:** N 1,543, MAE 0.270886692 Å³, RMSE 0.490777865 Å³, max 4.559474133 Å³, R² = 0.985128630. Original-convention Morfeus reduces all 20 headline residuals to ≤ 0.006275429 Å³. Approximate agreement is reproduced; exact convention equivalence is not.
7. **Confidence:** High for measured results.
8. **Limitations:** Default density and center differ; aggregate accuracy does not establish physical or predictive validity.

## K03 — Entire buried-volume family validated

1. **StericX claim:** Study 004 describes the entire Kraken vbur family as reproduced.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Every claimed field/reduction must be checked separately; original source includes distal/total volumes, ratios and Boltzmann reductions beyond eight geometric fields.
4. **Implementation location:** `studies/study_004_vbur_family.py::FAMILY; primary source boltzproperties/mmproperties`
5. **Independent validation:** Enumerate primary properties and compare every implemented matching unweighted reduction: 56 metrics across 14 fields × 4 reductions.
6. **Result — TERMINOLOGY ISSUE:** Eight absolute BV fields are a subset; even 56 current comparisons do not include source-energy Boltzmann or unimplemented distal/total families.
7. **Confidence:** High.
8. **Limitations:** Unimplemented source properties are OUT OF SCOPE for numerical SUT reproduction; missing weights remain UNCERTAIN rather than verified.

## K04 — Matched published Sterimol convention

1. **StericX claim:** Study 004 reports six Sterimol extrema under matched Kraken coordination-axis convention.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Original virtual Pd→P axis with raw bond vector center, H radius 1.09 Å, 3,600 directions, +0.40 Å L.
4. **Implementation location:** `src/geometry/sterimol.rs; src/descriptors.rs; studies/study_004_sterimol.py`
5. **Independent validation:** All 31,721 conformers in independent primary-convention Morfeus; explicit same center/radii and radius-only variants; full minima/maxima/other reductions.
6. **Result — INCORRECT:** Strict matched-convention wording is false: native H radius 1.20 Å, unit bond vector center and 360 directions differ. Six native extrema MAEs 0.0963–0.1124 Å; max L error 0.98638 Å against published values.
7. **Confidence:** High for conventions and measurements.
8. **Limitations:** A different documented descriptor convention can be legitimate; choosing one requires scientific intent, not tuning against these outputs.

## K05 — Pyramidalization formula on exported geometries

1. **StericX claim:** Native P/alpha reproduce the intended Morfeus construction for ordinary three-neighbor Kraken geometries.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Same donor, coordinates and three selected neighbors; P dimensionless and alpha degrees.
4. **Implementation location:** `src/geometry/pyramidalization.rs`
5. **Independent validation:** Independent contemporary Morfeus on all 31,721 exported conformers; nearest-three and SDF-topology variants independently agree in this corpus.
6. **Result — VERIFIED WITH NUMERICAL LIMIT:** Max identical-geometry |ΔP| = 3.2010947048632943e−7; |Δalpha| = 4.8856123654239525e−5 degrees.
7. **Confidence:** High within tested corpus.
8. **Limitations:** Does not establish undocumented degeneracies or every automatic neighbor policy; see geometry audit. Published-value extrema are a separate claim K09.

## K06 — Electronic-LMO explanation of P–H residuals

1. **StericX claim:** STUDY 004_RESIDUAL attributes the Kraken DFT center to an xTB localized orbital and says it coincides with the geometric center for tertiary phosphines.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Actual published DFT code sums raw P − neighbor vectors then normalizes; native sums individually normalized vectors. Neither is an electronic orbital center in this DFT path.
4. **Implementation location:** `docs/study_004/STUDY_004_RESIDUAL.md; src/geometry/buried_volume.rs::lone_pair_direction`
5. **Independent validation:** Immutable source inspection plus original center/native center explicit comparisons across 590 diagnostic conformers, including all 20 headline ensembles.
6. **Result — INCORRECT:** Largest headline case ID 1796 has zero P–H bonds, center angle≈ 9.715 degrees and 4.55947 Å³ discrepancy; original center/density resolves it to 0.00418 Å³. Many major cases have no P–H bonds.
7. **Confidence:** High.
8. **Limitations:** Finding is specific to this DFT path; other xTB/LMO workflows are not inferred identical.

## K07 — Orientation and finite-grid identity

1. **StericX claim:** Matched Kraken convention implies identical total, near/far and finite-grid outputs.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Original source retains last of three frame evaluations for total/near/far; native retains first. Quadrant/octant extrema combine all three.
4. **Implementation location:** `src/geometry/buried_volume.rs::compute_from_center; primary morfeus_properties`
5. **Independent validation:** Direct reference variants explicitly match first orientation before native comparison, then change only density or original center.
6. **Result — VERIFIED WITH NUMERICAL LIMIT:** Expanded matched-default comparison covers 2,372 cases: total BV max error 0.01166136723496436 Å³; percent BV max error 0.006494365624803322 percentage points; adjacent-quadrant max error 0.011656810711858867 Å³. These are approximately one grid cell. The earlier 581-case sample had a smaller total error; expansion is retained.
7. **Confidence:** High for measured finite-grid effects.
8. **Limitations:** Matched frame/density is necessary; first-versus-last convention is scientifically equivalent only in the continuum, not bitwise on finite grids.

## K08 — Shipped database identity

1. **StericX claim:** Committed table contains native descriptors and declares bond-axis Sterimol, rather than containing the published coordination-axis values.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Artifact provenance/declared axis must be distinguished from scientific numerical accuracy.
4. **Implementation location:** `data/ligand_db/kraken_phosphines.manifest.json; src/commands/db.rs`
5. **Independent validation:** Inspect frozen manifest, source table and consumers separately from all fresh scientific comparisons.
6. **Result — VERIFIED:** Frozen manifest explicitly declares computed native values and sterimol_axis=bond. This establishes identity/declared convention only.
7. **Confidence:** High for artifact metadata.
8. **Limitations:** Does not independently verify every historical packed record or authorize exchanging bond-axis and coordination-axis model descriptors.

## K09 — Published pyramidalization residual is coordinate/precision only

1. **StericX claim:** Study 005 interprets four extrema residuals as essentially export-coordinate/Rust-precision error.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Any attribution must explain molecule-specific residual magnitude with independent matched-geometry calculations.
4. **Implementation location:** `docs/study_005; studies/study_005_pyramidalization.py`
5. **Independent validation:** Full corpus Morfeus and native; 1,028,096 rounding corner calculations on 251 major outlier conformers.
6. **Result — UNCERTAIN:** Native/published max alpha discrepancy 1.155986692 degrees is reproduced independently (1.156001260 degrees). Same geometry native error≤ 4.89e−5 degrees; sampled rounding deviation≤ 0.010351413 degrees overall and 0.008927286 degrees for worst ID 1290. Simple precision attribution is unsupported. A conservative continuous rounding-box bound gives maximum |ΔP| of 0.000623445508668 across all 251 examined conformers, including both ID 1290 cases; its published P discrepancy is about 0.007218735, assuming unchanged ensemble and atom identity.
7. **Confidence:** High that tested numerical effects are insufficient; provenance cause uncertain.
8. **Limitations:** Corner sensitivity is not a rigorous interior bound. Original full precision logs/ensembles and historical software needed; do not simply blame published dataset.

## K10 — Sterimol L/B5 exact axis support

1. **StericX claim:** Comments describe B5 as exact maximum radial support and imply an exact alignment for all valid finite axes.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** L=max axial support+0.40; B5=max perpendicular distance+radius on intended axis.
4. **Implementation location:** `src/geometry/sterimol.rs::compute_with_dummy; glam0.30.10::Quat::from_rotation_arc`
5. **Independent validation:** Independent dot/perpendicular formulas on 31,713 successful native axes, corroborated with direct Morfeus on both worst real conformers.
6. **Result — INCORRECT:** L error 0.0036931285331647246 Å at KRAKEN:16:31923; B5 error 0.004176904379527002 Å at KRAKEN:1407:54385. Near parallel/antiparallel approximation discards a small real tilt.
7. **Confidence:** High; exact original failing inputs and independent results retained.
8. **Limitations:** Most errors≈ 5e−7 Å; magnitude is small relative to many ligand differences but nonzero and larger than claimed exactness. No production fix made.

## K11 — Geometric donor-neighbor inference matches chemical topology

1. **StericX claim:** Geometric proximity reliably supplies three bonded substituents for intended library phosphines.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Known molecular topology should not gain nonbonded nearby P/Si atoms solely because of generous radius cutoff.
4. **Implementation location:** `src/geometry/buried_volume.rs::donor_neighbor_indices; studies/study_004_scaled.py prefilter`
5. **Independent validation:** Fresh independent SDF graph and original Kraken connectivity compared with native; distances and element thresholds preserved for failures.
6. **Result — INCORRECT:** Eight conformers from IDs 1281/1907 get a fourth neighbor (P or Si); SDF and original primary connectivity have three. Both BV and coordination Sterimol error.
7. **Confidence:** High for explicit disagreement with two independent connectivity conventions.
8. **Limitations:** Chemical bond assignment still inferred; no experimental bond-order truth claimed. The historical prefilter mirrored native cutoff and hid this failure mode.

## K12 — Symmetric zero adjacent-quadrant difference implies invalid geometry

1. **StericX claim:** Native BV rejects nonzero occupancy with global max_delta_qvbur=0 as degenerate.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** A symmetrically occupied valid sphere can have zero anisotropy; finite grid zero does not prove an invalid frame.
4. **Implementation location:** `src/geometry/buried_volume.rs::compute_from_center`
5. **Independent validation:** Frozen actual PH3-like ligand 1299 conformer 54318; independent geometry/frame and Morfeus occupancy comparison.
6. **Result — INCORRECT:** Native raises symmetric-zero error for valid three-neighbor input; independent calculation succeeds.
7. **Confidence:** High with actual input and reference output.
8. **Limitations:** Do not infer every zero is valid; validate geometric degeneracy independently. Synthetic symmetric cases examined in geometry audit.

## K13 — Full Kraken Boltzmann reproduction

1. **StericX claim:** Full reproduction could extend published ligand-level Boltzmann values from available geometries.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Need per-conformer source corrected free energies/populations, temperature/degeneracy, selection and identity mapping. Published aggregates are comparison outputs only.
4. **Implementation location:** `studies/study_004_*; original read_ligand weighting`
5. **Independent validation:** Current/v2 API schemas/endpoints on five conformers, alternate exports, original repo/release tree, full SI archive, updated Sigman repo checked and frozen.
6. **Result — UNCERTAIN:** Historical weights not recovered. API returns data:null; SI contains aggregate tables; updated repo 33 new ensembles aredifferentcalculations. Full corpus Boltzmann not reproduced.
7. **Confidence:** High in checked-source limitation, not universal absence.
8. **Limitations:** Obtain original source energies; synthetic Boltzmann verification elsewhere does not validate historical populations.

## K14 — Independent reference reproduces every published Sterimol extremum

1. **StericX claim:** Matching primaryconventions suffices to exactly reproduce all published Sterimol values fromcurrentAPI geometries.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Same original geometry/ensemble, axis, radii and historical dependency implementation are required.
4. **Implementation location:** `Primary API and original Kraken/SI artifacts; a provenance claim, not exclusively a native-code claim`
5. **Independent validation:** Contemporary Morfeus primary convention on all 31,721 geometries, complete ensembles, all outliers retained.
6. **Result — UNCERTAIN:** Remaining largest errors: L max 0.955032324 Å (ID 1290), B5 min 0.649228687 Å (ID 1805), B1 max 0.274270968 Å (ID 560). Bothnativeandindependentreferencedisagreewithsomepublished extrema.
7. **Confidence:** High in measured disagreement; attribution uncertain.
8. **Limitations:** Need original calculation provenance before declaring either the public table or current implementation wrong.

## K15 — Native pyramidalization machine-precision claim

1. **StericX claim:** Study 005 says the native kernel agrees with Morfeus to 4.4e−16 for P and 2.8e−14 degrees for alpha.
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Expected convention:** Those numbers must describe native Rust f32 output if attributed to the native kernel, not a Python float64 identity check.
4. **Implementation location:** `docs/study_005/STUDY_005.md:21; src/geometry/pyramidalization.rs`
5. **Independent validation:** Actual native values and IEEE bits frozen for all 31,721 conformers, compared independently with Morfeus on the identical exported geometry.
6. **Result — INCORRECT:** Measured maxima are 3.2010947048632943e−7 for P and 4.8856123654239525e−5 degrees for alpha. A Python algebraic identity at float64 precision cannot support the stated native-output precision.
7. **Confidence:** High, with full corpus native and independent records.
8. **Limitations:** The measured native errors are small in ordinary geometries; this finding corrects numerical-precision attribution, not every aspect of the underlying definition.
