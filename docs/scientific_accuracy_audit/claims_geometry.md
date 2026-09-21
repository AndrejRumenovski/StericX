# Geometry claim inventory

Each result is independently classified. An implementation-level agreement is not evidence of chemical predictive validity. Detailed evidence is in [geometry/REPORT.md](geometry/REPORT.md), with equations and primary-source limitations in [geometry/METHODS.md](geometry/METHODS.md).

## G01: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Raw Sterimol L is the maximum axial atomic-sphere extent. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact expected convention:** max_i[(r_i-o)·u + rho_i], excluding a real dummy. Å.
4. **Implementation:** `src/geometry/sterimol.rs::compute/params_from_projection`
5. **Independent validation:** Analytic support equation in independent float64 + Morfeus, 31 finite base cases.
6. **Result:** Max error 6.994274599492201e-7 Å for bond-axis base cases.
7. **Confidence:** High for controlled valid inputs
8. **Limitations:** Restricted base cases are not a global accuracy bound. Real Kraken near-axis cases expose a separate 0.00172232 Å L alignment error (G39). Does not validate automatic axis chemistry, unknown radii, large-coordinate overflow, or historical original radius parameterization.

## G02: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** The documented one-degree B1 estimate approximates the continuous minimum support radius within the measured errors. Scoped numerical audit of the one-degree scan documented in the kernel; not a claim that the source promises exact continuous minimization.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact expected convention:** min over all transverse unit directions of max_i(p_i·v+rho_i).
4. **Implementation:** `src/geometry/sterimol.rs::params_from_projection`
5. **Independent validation:** Analytic finite candidate enumeration of support-envelope stationary points/intersections; minimized four-atom witness.
6. **Result:** SUT 1.7261793613433838 Å vs exact 1.7000000000000002 Å; a 0.0261793613433836 Å overestimate from one-degree sampling.
7. **Confidence:** High
8. **Limitations:** Approximation is documented in kernel comments. Exact B1 equality/strict invariance is false; error scales with transverse size. A retained full Kraken case has error 0.0416578313 Å; the minimized witness is not the audit-wide maximum. See G39 for the distinct alignment contribution.

## G03: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** B5 is the maximum transverse atomic-sphere extent. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact expected convention:** max_i(||p_i||+rho_i), Å.
4. **Implementation:** `src/geometry/sterimol.rs::params_from_projection`
5. **Independent validation:** Independent analytic radial support; compare Morfeus sampled maximum separately.
6. **Result:** Max bond-axis base error 5.23377122085833e-7 Å; Morfeus B5 itself uses an approximate angular scan.
7. **Confidence:** High
8. **Limitations:** The base maximum is not a global bound: a real near-antiparallel Kraken case has B5 error 0.00276961144 Å (G39). Extreme coordinates also require separate limits.

## G04: SUPPORTED

1. **StericX claim / audited proposition:** Coordination-axis L applies the historical +0.40 Å convention. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html ; https://aarontools.readthedocs.io/en/latest/api/substituent.html
3. **Exact expected convention:** Raw API L plus 0.40 Å in CLI coordination mode; bond mode raw.
4. **Implementation:** `src/descriptors.rs::sterimol_for_conformer; src/cli.rs::STERIMOL_L_CORRECTION`
5. **Independent validation:** Read frozen source and compare raw API / Morfeus L_value_uncorrected and L_value.
6. **Result:** Arithmetic/convention placement agrees; raw API and CLI use distinct conventions.
7. **Confidence:** High implementation; moderate historical attribution
8. **Limitations:** Original Verloop chapter not fully accessible. The correction is a historical parameter convention, not a universal physical requirement for arbitrary metal dummies.

## G05: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Explicit hydrogen spheres contribute to Sterimol and real dummy is excluded. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact expected convention:** All non-dummy atoms with their chosen radii; virtual dummy includes real donor.
4. **Implementation:** `src/geometry/sterimol.rs::compute/compute_with_dummy`
5. **Independent validation:** Independent explicit-H graph molecules and envelope calculation.
6. **Result:** Ordinary geometries agree within reported B1 scan/f32 limits.
7. **Confidence:** High
8. **Limitations:** Kraken DFT code uses H radius1.09 Å; StericX uses1.20 Å, so published reproduction is separate.

## G06: INCORRECT

1. **StericX claim / audited proposition:** Default bond-axis descriptors are independent of atom ordering. 
2. **Primary/reference source:** Geometric invariance requirement; https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact expected convention:** Same chemical structure and chosen chemical axis should produce the same descriptor.
4. **Implementation:** `src/descriptors.rs::detect_donor; src/geometry/buried_volume.rs::bonded_neighbors`
5. **Independent validation:** Frozen actual CLI on all six permutations of three equally distant chemically distinct substituents.
6. **Result:** L changes6.0699997→3.5 Å and B5 changes3.5→6.0699997 Å because the tie selects the first atom index.
7. **Confidence:** High
8. **Limitations:** The underlying explicitly selected-axis kernel is invariant to reordering within f32 limits. This falsifies an unrestricted auto-axis invariance claim, not the mathematics for a specified axis.

## G07: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Sterimol is rigid-rotation/translation invariant. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact expected convention:** Physical continuous descriptors invariant for fixed chemical axis.
4. **Implementation:** `src/geometry/sterimol.rs::compute/compute_with_dummy`
5. **Independent validation:** 6 representatives×100 seeded random rotations plus X/Y/Z rotations and random translations.
6. **Result:** B1 max deviation0.019944190979003906 Å; L≤1.0013580322265625e-5 Å, B5≤4.76837158203125e-6 Å.
7. **Confidence:** High for measured range
8. **Limitations:** Rotations include modest translations ≤50 Å; no universal bound. B1 scan phase depends on alignment. A targeted exact 90-degree rotation of real near-axis Kraken cases changes L by 0.00172233582 Å and B5 by 0.00276851654 Å (G39), despite exact independent invariance.

## G08: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Buried volume integrates the union of atomic spheres inside a sphere. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html ; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** V=integral of union membership inside R; %Vbur=100V/(4πR³/3).
4. **Implementation:** `src/geometry/buried_volume.rs::occupied_volumes`
5. **Independent validation:** Independent float64 union-of-balls integration and Morfeus under matched center/radii/grid convention.
6. **Result:** 23 finite base values max total-volume difference3.6215270000639066e-6 Å³; default lattice includes15,408 points.
7. **Confidence:** High for compatible default cases
8. **Limitations:** Numerical quadrature error is larger than implementation residual; invalid symmetry rejection and unknown radii are distinct failures.

## G09: VERIFIED

1. **StericX claim / audited proposition:** The default sphere radius3.5 Å follows a published descriptor convention. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html
3. **Exact expected convention:** Sphere centered on selected real/virtual coordination point, R3.5 Å.
4. **Implementation:** `src/geometry/buried_volume.rs::BuriedVolumeConfig; src/cli.rs`
5. **Independent validation:** Official SambVca manual and Morfeus documentation.
6. **Result:** Convention agrees.
7. **Confidence:** High
8. **Limitations:** This empirically useful convention does not prove a universal physically optimal sphere radius.

## G10: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Bondi radius scaling1.17 follows the reference convention. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html ; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** Every included vdW radius multiplied by1.17 before intersection.
4. **Implementation:** `src/geometry/buried_volume.rs::aligned_atoms`
5. **Independent validation:** Official manual plus independent scaled-radius calculations.
6. **Result:** Scaling agrees on supported atoms and finite ordinary values.
7. **Confidence:** High
8. **Limitations:** Input radii table and overflow behavior must be validated separately.

## G11: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Hydrogens are excluded from buried-volume occupancy by default. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html ; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** H may define connectivity/frame but contributes no sphere unless include_hydrogens=true.
4. **Implementation:** `src/geometry/buried_volume.rs::aligned_atoms`
5. **Independent validation:** Explicit-H methylphosphine/trimethylamine protocol pairs and PH3 analytic lens witness.
6. **Result:** Exclusion/inclusion matches controlled Morfeus reference.
7. **Confidence:** High
8. **Limitations:** Removing donor-bound H from the frame would be a different, incorrect chemical operation.

## G12: INCORRECT

1. **StericX claim / audited proposition:** All supported radii are exactly Morfeus Bondi radii. 
2. **Primary/reference source:** https://doi.org/10.1021/j100785a001 ; frozen/reference_sources/morfeus/data.py
3. **Exact expected convention:** Declared table equality, not only the label Bondi-style.
4. **Implementation:** `src/geometry/xyz.rs::van_der_waals_radius`
5. **Independent validation:** Read independent table; per-element synthetic witnesses.
6. **Result:** Br is1.85 Å in SUT vs1.83 in Morfeus; B1.92 in SUT vs Morfeus fallback2.0.
7. **Confidence:** High table difference
8. **Limitations:** The original Bondi table was not available for direct inspection, so the Br discrepancy was not adjudicated against that table; B has no Morfeus Bondi entry; Mantina et al., DOI 10.1021/jp8111556, explicitly gives the later 1.92 Å extension. Difference alone does not identify physical truth.

## G13: UNCERTAIN

1. **StericX claim / audited proposition:** Unknown-element radius fallbacks are chemically conservative. 
2. **Primary/reference source:** https://doi.org/10.1021/j100785a001
3. **Exact expected convention:** An arbitrary fallback must not be equated with an element-specific physical radius.
4. **Implementation:** `src/geometry/xyz.rs::van_der_waals_radius/covalent_radius`
5. **Independent validation:** He/Xe/Na/Se/Fe toy witnesses; official/reference tables.
6. **Result:** SUT1.8 Å can be both too large(He) and too small(Na/Xe); max observed BV difference6.072736800747158 Å³,3.3813625981379403 percentage points.
7. **Confidence:** High numerical result; low chemical justification
8. **Limitations:** Morfeus also falls back2.0 for absent B/Fe entries; it is not ground truth. Broad chemical safety of either fallback lacks evidence.

## G14: INCORRECT

1. **StericX claim / audited proposition:** Audit proposition: the filled integration grid exactly reproduces Morfeus for every accepted density. Exact-compatibility audit question, not a quoted universal promise in the source. The documented SUT grid is evaluated against the independent reference grid.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** Float64 linspace and norm≤R in Morfeus versus f32 SUT arithmetic; grid endpoints included.
4. **Implementation:** `src/geometry/buried_volume.rs::integration_grid`
5. **Independent validation:** Four density levels plus odd-n grid and representable-float sphere-boundary perturbations.
6. **Result:** Default n32 counts/membership agree for tested cases; odd grids reveal region-boundary differences, fine grids can differ at sphere lattice boundaries.
7. **Confidence:** High
8. **Limitations:** This rejects unrestricted exact all-density equivalence, not the shared continuum target or measured default-grid agreement. Floating-point sphere-boundary membership and region conventions differ; see convergence and partition records.

## G15: VERIFIED

1. **StericX claim / audited proposition:** Voxel atomic occupancy includes the vdW boundary. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html ; frozen Morfeus cKDTree query_ball_point implementation
3. **Exact expected convention:** distance²≤scaled_radius², corresponding to a closed ball.
4. **Implementation:** `src/geometry/buried_volume.rs::occupied_volumes`
5. **Independent validation:** Frozen one-point tests at1f32 and nextafter inside/outside, independent exact-real comparison.
6. **Result:** Inside/on occupied; next representable outside unoccupied.
7. **Confidence:** High for tested boundary
8. **Limitations:** Open/closed spheres have equal continuum volume; this establishes a lattice convention, not a uniquely physical inequality.

## G16: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Sphere boundary handling includes points on R. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** Retain ||point||≤R.
4. **Implementation:** `src/geometry/buried_volume.rs::integration_grid`
5. **Independent validation:** Private-grid observation at n5 with R and adjacent representable radii.
6. **Result:** Boundary points are present; outside lattice points rejected.
7. **Confidence:** High
8. **Limitations:** f32 rounding can move mathematically marginal points. Atom centers outside R must still contribute when their sphere intersects R, confirmed by boundary-atom witnesses.

## G17: INCORRECT

1. **StericX claim / audited proposition:** Audit proposition: quadrants use exactly the same finite-grid regional convention as Morfeus at every accepted density. Reference-convention compatibility question, not an assertion that an open or closed continuum boundary is physically preferable.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** Canonical xy signs (++),(-+),(--),(+-); Morfeus excludes zero planes.
4. **Implementation:** `src/geometry/buried_volume.rs::quadrant_index/occupied_volumes`
5. **Independent validation:** Every individual quadrant in every orientation, independent/Morfeus comparison.
6. **Result:** Default even-grid base values agree within f32 rounding. n5 qvbur_min differs9.621127681084193 Å³; odd-grid assignment differs structurally.
7. **Confidence:** High
8. **Limitations:** Default even-grid agreement remains verified within the reported numerical limits. At odd grids SUT assigns zeros positive and directly normalizes quadrants; Morfeus excludes zero planes and sums octants. Both approximate the same continuum regions without exact numerical convention equivalence.

## G18: INCORRECT

1. **StericX claim / audited proposition:** Audit proposition: octants use exactly the same finite-grid regional convention as Morfeus at every accepted density. Reference-convention compatibility question; it does not identify either implementation as physical truth.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** Four xy quadrants for positive and negative z; explicit Morfeus ID remapping.
4. **Implementation:** `src/geometry/buried_volume.rs::octant_index/occupied_volumes`
5. **Independent validation:** Every individual octant per orientation compared, canonical signs preserved.
6. **Result:** Even-grid agreement is strong; zero-plane regions at odd n differ.
7. **Confidence:** High
8. **Limitations:** Default even-grid numerical agreement remains favorable. Odd-grid zero-plane assignment differs. Negative-z numeric IDs also differ, but the harness correctly remaps them; this classification concerns the actual boundary convention, not a mistaken ID comparison.

## G19: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Near and far describe donor-facing and distal hemispheres. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html ; Kraken official SI page12
3. **Exact expected convention:** Center→donor is negative z; near z<0, far z≥0 in SUT.
4. **Implementation:** `src/geometry/buried_volume.rs::coordinate_basis/occupied_volumes`
5. **Independent validation:** Independent coordinate frame plus per-octant sums.
6. **Result:** Hemisphere naming/orientation agrees under the stated coordinate convention.
7. **Confidence:** High
8. **Limitations:** This is not the outside-sphere distal molecular volume. Odd-grid normalization creates partition inconsistency.

## G20: INCORRECT

1. **StericX claim / audited proposition:** Audit proposition: each finite-grid regional partition adds exactly to the independently estimated total for all valid configurations. Numerical conservation audit question, distinct from exact additivity of the continuum integral; not a quoted all-grid guarantee.
2. **Primary/reference source:** Continuum additivity of disjoint integrals; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** Regional volume sum should converge to total, with quantified finite-grid discrepancy.
4. **Implementation:** `src/geometry/buried_volume.rs::occupied_volumes`
5. **Independent validation:** Independent conservation check across densities.
6. **Result:** Triphenylphosphine coarse total53.1568947 versus near+far59.7332687 Å³ (+12.37%); default even-grid error≤rounding.
7. **Confidence:** High
8. **Limitations:** Default even-grid additivity holds within rounding in the campaign. The odd-grid discrepancy is a finite-estimator limitation, not a failure of continuum volume additivity. At density 0.0001 near+far remains 0.5790825 Å³ above total. Morfeus also has a coarse-grid gap of +3.4082601 Å³.

## G21: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Maximum adjacent-quadrant difference is reduced correctly. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html ; official Kraken DFT code
3. **Exact expected convention:** max over three plane orientations and four cyclic adjacent pairs of absolute difference.
4. **Implementation:** `src/geometry/buried_volume.rs::compute_from_center`
5. **Independent validation:** Independent explicit twelve-quadrant reduction.
6. **Result:** Default matched cases agree within f32 subtraction rounding.
7. **Confidence:** High
8. **Limitations:** Frame/center selection, zero-plane treatment and symmetry rejection are separate issues.

## G22: INCORRECT

1. **StericX claim / audited proposition:** A positive buried volume with zero quadrant asymmetry is unphysical/degenerate. 
2. **Primary/reference source:** Sphere intersection symmetry; independent analytic lens; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** An axisymmetric nonzero occupied set has equal quadrants and valid zero asymmetry.
4. **Implementation:** `src/geometry/buried_volume.rs::compute_from_center final rejection`
5. **Independent validation:** PH3 with three chemically bonded H; full-cover toy; small sphere wholly inside ordinary donor radius.
6. **Result:** PH3 correctly has31.937214517315244 Å³ lattice volume,17.782969885773625%,zero asymmetry; SUT errors despite valid frame. Full-cover100% also errors.
7. **Confidence:** High, constructive counterexamples
8. **Limitations:** PH3 CLI separately rejects lack of heavy reference atom; public BV API accepts H as reference then hits invalid symmetry rejection. Small-sphere case also reproduces through CLI with config.

## G23: UNCERTAIN

1. **StericX claim / audited proposition:** Audit proposition: default buried-volume precision supports scientific distinctions finer than its discretization error. Scientific-resolution audit question; no explicit application-specific resolution guarantee was identified in the cited kernel. The measured numerical floor does not itself establish a useful chemical effect size.
2. **Primary/reference source:** Analytic two-sphere intersection formula
3. **Exact expected convention:** Compare lattice to continuum, not only another lattice implementation.
4. **Implementation:** `src/geometry/buried_volume.rs::integration_grid/occupied_volumes`
5. **Independent validation:** Exact PH3 one-atom lens and four densities.
6. **Result:** Analytic31.669141713568127 Å³; default error+0.2680721841369511 Å³ or0.149265 percentage points; density.001 error−0.0284346976; .0001 error+0.0006218515 Å³.
7. **Confidence:** High for the measured lens/convergence errors; insufficient evidence for scientific resolution below them
8. **Limitations:** The analytic lens has default error 0.2680722 Å³, and four real/toy molecules differ from very fine totals by up to 0.1971054 Å³. Error cancellation for a particular descriptor difference, experimental relevance, and a universal resolution guarantee are not established. Agreement with another finite grid cannot supply that missing evidence.

## G24: SUPPORTED

1. **StericX claim / audited proposition:** Buried-volume values are rigid-transformation invariant over the tested moderate-coordinate campaign. Scoped numerical invariance proposition: six prespecified representatives and 618 transforms, not a universal finite-precision guarantee.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact expected convention:** Co-rotating molecular frame and center; finite-grid sensitivity quantified.
4. **Implementation:** `src/geometry/buried_volume.rs::coordinate_basis/occupied_volumes`
5. **Independent validation:** 618 transformed cases across six representatives, all public fields and raw region values.
6. **Result:** All nine public BV fields had exactly0 measured deviation on this campaign.
7. **Confidence:** High for campaign; limited generality
8. **Limitations:** Grid transitions can occur for other geometries; huge translations lose f32 geometric information and are not covered by this favorable result.

## G25: SUPPORTED

1. **StericX claim / audited proposition:** Three covalently bonded substituents include donor-bound hydrogen. 
2. **Primary/reference source:** https://doi.org/10.1039/B801115J ; explicit molecular graphs
3. **Exact expected convention:** SMILES graph bonds define chemically intended neighbors; H included.
4. **Implementation:** `src/geometry/buried_volume.rs::bonded_neighbors/donor_neighbor_indices`
5. **Independent validation:** 12 independent RDKit graph molecules plus primary/secondary phosphine contact witnesses.
6. **Result:** All12 standard graphs match; both H-bearing historical witnesses recover chemically bonded neighbors.
7. **Confidence:** High for test molecules
8. **Limitations:** Geometry-only inference is not universal chemistry; charge/bond order/coordination identity are ignored.

## G26: UNCERTAIN

1. **StericX claim / audited proposition:** Cordero covalent radii justify a universal1.3× bond cutoff. 
2. **Primary/reference source:** https://doi.org/10.1039/B801115J
3. **Exact expected convention:** Cordero provides empirical atomic radii, not universal binary bonding truth at1.3× sums.
4. **Implementation:** `src/geometry/buried_volume.rs::BOND_TOLERANCE_FACTOR/bonded_neighbors`
5. **Independent validation:** Primary source interpretation, threshold nextafter witnesses, independent graphs.
6. **Result:** Threshold behavior is deterministic and matches stated heuristic; scientific universality unsupported. Coincident atoms count as bonded before later geometry validation.
7. **Confidence:** High implementation; insufficient chemical universality
8. **Limitations:** No threshold retuning was performed. Stretching a synthetic bond beyond cutoff is not itself proof the chemistry is misclassified.

## G27: SUPPORTED

1. **StericX claim / audited proposition:** Automatic donor selection identifies the expected donors in the tested standard P/N graphs. Scoped test of automatic-selection behavior against independently supplied molecular graphs; general chemical donor availability remains unestablished.
2. **Primary/reference source:** https://doi.org/10.1039/B801115J ; CLI contract
3. **Exact expected convention:** Sole chosen element, or explicit index, plus exactly3 inferred neighbors and≥1heavy neighbor.
4. **Implementation:** `src/descriptors.rs::detect_donor`
5. **Independent validation:** Frozen CLI runs with phosphines/amines, multiple P, planar/coincident cases and PH3.
6. **Result:** Standard tested P/N donors work; multiple-element ambiguity is rejected; PH3 rejected for lack of heavy reference; chemistry beyond geometry unestablished.
7. **Confidence:** Moderate
8. **Limitations:** Explicit index can select any element; trivalent geometry does not establish lone-pair availability, charge state, aromaticity or catalytic coordination.

## G28: SUPPORTED

1. **StericX claim / audited proposition:** Current frame fixes the historical nearest-heavy error. 
2. **Primary/reference source:** Independent known SMILES connectivity; primary/secondary phosphine valence
3. **Exact expected convention:** Three bonded neighbors, including H, rather than three nearest heavy atoms.
4. **Implementation:** `src/geometry/buried_volume.rs::donor_neighbor_indices`
5. **Independent validation:** Reconstruct former selector independently; controls plus H and nearby nonbonded heavy atoms.
6. **Result:** Tertiary control center unchanged; primary/secondary wrong-frame centers displaced2.4302383/2.0563763 Å and BV descriptors differ strongly; current neighborhoods match graphs.
7. **Confidence:** High for reconstructed failure mechanism
8. **Limitations:** The exact six historically zero ligands (575,1485,1487,1490,1491,1495),20 conformers, were then independently reconstructed from primary SDF bond graphs; every old nearest-heavy frame again gives zero while current matched-frame values agree with independent references. Kraken reproduction conventions remain separate.

## G29: INCORRECT

1. **StericX claim / audited proposition:** The geometric virtual-metal placement exactly matches published Kraken. 
2. **Primary/reference source:** Official Kraken SI FigureS2/page10 and frozen DFT source audited in kraken/
3. **Exact expected convention:** Kraken DFT sums raw substituent displacements then normalizes; SUT sums unit bond directions.
4. **Implementation:** `src/geometry/buried_volume.rs::infer_lone_pair_direction`
5. **Independent validation:** Primary code/SI comparison, independent center construction.
6. **Result:** Conventions differ for unequal bond lengths; matching center-distance2.28 alone is insufficient.
7. **Confidence:** High
8. **Limitations:** Separate Kraken audit quantifies published-data consequences. Geometric placement is itself an approximation to actual coordination geometry.

## G30: UNCERTAIN

1. **StericX claim / audited proposition:** Planar fallback is a physically determined lone-pair direction. 
2. **Primary/reference source:** Geometric symmetry; no universal primary rule identified
3. **Exact expected convention:** Planar center has no unique opposite-bond-sum direction; plane normal sign is ambiguous without environment/electronics.
4. **Implementation:** `src/geometry/buried_volume.rs::infer_lone_pair_direction/center_clearance`
5. **Independent validation:** Planar and nearly planar synthetic structures, rotations and permutations.
6. **Result:** Fallback normal/maximum-clearance criterion is a StericX design choice; independent naive normalized-sum centers become unstable near exact planarity.
7. **Confidence:** High description; insufficient physical evidence
8. **Limitations:** Controlled-center kernel comparisons deliberately do not treat SUT-selected center as scientifically established.

## G31: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Pyramidalization P implements the determinant and acute-correction convention used by Morfeus 0.8.0. Runnable-reference implementation proposition. Morfeus attributes the measure to Radhakrishnan and Agranat; original-paper textual equivalence was not established because full-text access was unavailable.
2. **Primary/reference source:** https://doi.org/10.1007/BF00676621 ; https://digital-chemistry-laboratory.github.io/morfeus/api/morfeus.pyramidalization.html
3. **Exact expected convention:** |det(unit donor vectors)| with2−P acute correction.
4. **Implementation:** `src/geometry/pyramidalization.rs::compute`
5. **Independent validation:** Independent determinant calculation, Morfeus, tetrahedral/orthogonal/acute witnesses.
6. **Result:** Ordinary graph cases maxP error1.2938381854787906e-7; tetrahedral0.769800305 versus analytic0.7698003589.
7. **Confidence:** High for runnable reference; moderate original-paper definition
8. **Limitations:** Original paper full text not acquired. Degenerate sentinel behavior is not a valid P measurement.

## G32: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Pyramidalization alpha uses the mean signed plane-normal angle in degrees. 
2. **Primary/reference source:** https://doi.org/10.1007/BF00676621 ; Morfeus0.8 source
3. **Exact expected convention:** Three alpha values, acute bisector sign, degree conversion; ordinary planar90°, orthogonal0°.
4. **Implementation:** `src/geometry/pyramidalization.rs::compute`
5. **Independent validation:** Independent plane geometry and explicit acute/planar/tetrahedral cases.
6. **Result:** Ordinary graph maxerror6.122375019401716e-6°; acute example−82.717° is permitted.
7. **Confidence:** High in nondegenerate range
8. **Limitations:** Near-collinear geometry is assigned0 by SUT despite finite independent29.99992497364023°; sentinel must not be treated scientifically.

## G33: TERMINOLOGY ISSUE

1. **StericX claim / audited proposition:** P=1 corresponds to an ideal tetrahedral apex. 
2. **Primary/reference source:** Analytic tetrahedral unit vectors; https://doi.org/10.1007/BF00676621
3. **Exact expected convention:** Three tetrahedral directions determinant4/(3√3)≈0.76980036; orthogonal triad1.
4. **Implementation:** `src/geometry/pyramidalization.rs::PyramidalizationParams documentation`
5. **Independent validation:** Analytic determinant and frozen directAPI.
6. **Result:** The comment saying1 for ideal tetrahedral-like apex is misleading; tests compute tetrahedral0.7698.
7. **Confidence:** High
8. **Limitations:** Tetrahedral-like is imprecise; no universal physical normalization to1 for tetrahedral geometry is established.

## G34: VERIFIED WITH NUMERICAL LIMIT

1. **StericX claim / audited proposition:** Pyramidalization is neighbor-order invariant. 
2. **Primary/reference source:** https://doi.org/10.1007/BF00676621
3. **Exact expected convention:** Geometric three-vector measure independent of permutation.
4. **Implementation:** `src/geometry/pyramidalization.rs::compute`
5. **Independent validation:** 120 molecular permutations with mapped neighbor identities.
6. **Result:** No scientifically meaningful order variation observed; full residuals in invariance.csv.
7. **Confidence:** High for tested domain
8. **Limitations:** Planar/degenerate numerical branches need explicit handling.

## G35: INCORRECT

1. **StericX claim / audited proposition:** Invalid or nearly degenerate geometry produces an explicit scientific failure. 
2. **Primary/reference source:** Mathematical nondefinition of directions/plane normals
3. **Exact expected convention:** Do not silently interpret sentinel zero as a measured descriptor.
4. **Implementation:** `src/geometry/sterimol.rs::compute; src/geometry/pyramidalization.rs::compute`
5. **Independent validation:** Coincident, short-axis, nearly collinear, NaN/Inf input witnesses.
6. **Result:** Sterimol and pyramidalization direct APIs returnzero sentinels; BV errors. Near-collinear alpha29.999925° becomes0°; short nonzeroaxis givesallzeroSterimol.
7. **Confidence:** High for unrestricted explicit-error claim
8. **Limitations:** Zero sentinel behavior is documented at kernel level. CLI combined-descriptor validation often rejects these cases; callers of public kernels remain exposed.

## G36: INCORRECT

1. **StericX claim / audited proposition:** Finite input coordinates always produce finite scientific outputs. 
2. **Primary/reference source:** IEEE754 finite arithmetic range; analytic norm
3. **Exact expected convention:** Overflow must be detected if finite output is promised.
4. **Implementation:** `src/geometry/sterimol.rs::params_from_projection`
5. **Independent validation:** Frozen finite1e38 Å remote atom witness.
6. **Result:** B5 becomes+Inf for wholly finite positions/radii; observer recordsnull plus explicitnonfinite map.
7. **Confidence:** High
8. **Limitations:** Deliberately extreme and not ordinary chemistry. It identifies an API validation limit, not evidence that f64 is needed for typical coordinates.

## G37: UNCERTAIN

1. **StericX claim / audited proposition:** f32 is sufficiently accurate for ordinary geometry in this audit. 
2. **Primary/reference source:** Independent f64 equations and Morfeus under controlled protocol
3. **Exact expected convention:** Separate rounding error from grid/angle discretization and coordinate convention.
4. **Implementation:** `src/geometry/*.rs`
5. **Independent validation:** Graph molecules, rigid campaign, f32-rounded input re-evaluated in f64.
6. **Result:** Restricted ordinary-case L/B5 errors are mostly 1e-6–1e-5 Å and default BV residuals a few 1e-6 Å³. The targeted real near-Z cases instead have ~0.002 Å alignment errors. Scientific sufficiency requires an application-specific tolerance or effect-size assessment.
7. **Confidence:** High measured errors; insufficient evidence for universal scientific sufficiency
8. **Limitations:** Real near-axis Kraken geometries have larger alignment errors (G39), beyond ordinary rounding. Large translations 1e7–1e20 destroy coordinate detail; nearly degenerate axes amplify rounding. No blanket recommendation to switch f64 is made.

## G38: OUT OF SCOPE

1. **StericX claim / audited proposition:** Geometry descriptor agreement validates predictive/experimental chemistry. 
2. **Primary/reference source:** Scientific identifiability; held-out experimental validation required
3. **Exact expected convention:** Numerical descriptor compatibility and chemistry prediction are different hypotheses.
4. **Implementation:** `README scientific summaries; studies002/004/005`
5. **Independent validation:** Separate scope of independent computation from model/generalization evidence.
6. **Result:** This geometry audit establishes no experimental predictive validity.
7. **Confidence:** High
8. **Limitations:** Requires held-out chemical/experimental evidence evaluated in model audit and prospective experiments.

## G39: INCORRECT

1. **StericX claim / audited proposition:** B5 is the exact radial support about the supplied attachment axis, apart from ordinary f32 arithmetic rounding. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html ; https://docs.rs/glam/0.30.10/src/glam/f32/sse2/quat.rs.html
3. **Exact expected convention:** Projection onto the supplied unit axis must remain invariant under a rigid change of global Cartesian coordinates; a near-Z axis must not be replaced by Z.
4. **Implementation:** `src/geometry/sterimol.rs::compute/compute_with_dummy, Quat::from_rotation_arc in frozen glam 0.30.10`
5. **Independent validation:** Independent float64 dot/norm projections on identical SUT-rounded coordinates/radii/center; two real Kraken conformers minimized to donor plus active atom; exact 90-degree X coordinate permutation.
6. **Result:** Full-case L is 6.859999656677246 versus independent 6.858277336508549 Å (963/48394). Full-case B5 is 7.116985321044922 versus independent 7.114215709601416 Å (1075/49864). Exact coordinate rotation reduces corresponding errors to 1.56e-8 and 1.09e-6 Å; original minimized witnesses retain the failure.
7. **Confidence:** High, preserved full and minimized cases plus targeted invariance experiment
8. **Limitations:** The supplied axes are mathematically well defined and the source geometries are real DFT structures. glam documents approximate near-singular alignment (|dot| > 1−2 f32 EPSILON); this is a StericX accuracy-contract issue, not an undocumented glam defect. It is distinct from B1 angular discretization and does not prove a chemistry prediction is materially affected. See geometry/alignment/ and its rotation/ subdirectory.

