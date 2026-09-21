# Scientific claim inventory

147 claims and explicitly identified audit propositions, each with a separate scientific classification. Frozen system: commit `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7`.

This inventory examines the README, scientific/validation studies, code comments, model documentation and tests as claims to investigate. They are not scientific reference truth. Some entries test a stronger interpretation or a numerical-domain property rather than quote an explicit promise; those distinctions and existing documentation caveats are retained. A high-confidence negative finding is not evidence that every normal workload fails.

Classification uses the requested seven categories. **VERIFIED** requires independent support for the stated scope; **VERIFIED WITH NUMERICAL LIMIT** adds a measured precision/discretization limit; **SUPPORTED** has favorable but incomplete evidence; **UNCERTAIN** has missing or conflicting evidence; **INCORRECT** has a constructive contradiction; **TERMINOLOGY ISSUE** concerns scientific/statistical wording; **OUT OF SCOPE** cannot be established computationally from available evidence.

The counts below count propositions, not independent experiments or a percentage of scientific correctness. Related claims can share evidence without being merged into one PASS.

| Classification | Claims |
|---|---:|
| VERIFIED | 15 |
| VERIFIED WITH NUMERICAL LIMIT | 49 |
| SUPPORTED | 20 |
| UNCERTAIN | 17 |
| INCORRECT | 38 |
| TERMINOLOGY ISSUE | 5 |
| OUT OF SCOPE | 3 |

Evidence reports: [overview](SCIENTIFIC_ACCURACY_AUDIT.md), [geometry](geometry/REPORT.md), [kinetics/conformers](kinetics/KINETICS_CONFORMERS.md), [Kraken](kraken/REPORT.md), [models](models/REPORT.md). Machine-readable records: [claims.json](claims.json).

## Geometry

### G01 — Raw Sterimol L is the maximum axial atomic-sphere extent.

1. **StericX claim / audited proposition:** Raw Sterimol L is the maximum axial atomic-sphere extent. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact convention expected:** max_i[(r_i-o)·u + rho_i], excluding a real dummy. Å.
4. **Implementation location:** `src/geometry/sterimol.rs::compute/params_from_projection`
5. **Independent validation method:** Analytic support equation in independent float64 + Morfeus, 31 finite base cases.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Max error 6.994274599492201e-7 Å for bond-axis base cases.
7. **Confidence:** High for controlled valid inputs
8. **Limitations:** Restricted base cases are not a global accuracy bound. Real Kraken near-axis cases expose a separate 0.00172232 Å L alignment error (G39). Does not validate automatic axis chemistry, unknown radii, large-coordinate overflow, or historical original radius parameterization.

### G02 — The documented one-degree B1 estimate approximates the continuous minimum support radius within the measured errors.

1. **StericX claim / audited proposition:** The documented one-degree B1 estimate approximates the continuous minimum support radius within the measured errors. Scoped numerical audit of the one-degree scan documented in the kernel; not a claim that the source promises exact continuous minimization.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact convention expected:** min over all transverse unit directions of max_i(p_i·v+rho_i).
4. **Implementation location:** `src/geometry/sterimol.rs::params_from_projection`
5. **Independent validation method:** Analytic finite candidate enumeration of support-envelope stationary points/intersections; minimized four-atom witness.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — SUT 1.7261793613433838 Å vs exact 1.7000000000000002 Å; a 0.0261793613433836 Å overestimate from one-degree sampling.
7. **Confidence:** High
8. **Limitations:** Approximation is documented in kernel comments. Exact B1 equality/strict invariance is false; error scales with transverse size. A retained full Kraken case has error 0.0416578313 Å; the minimized witness is not the audit-wide maximum. See G39 for the distinct alignment contribution.

### G03 — B5 is the maximum transverse atomic-sphere extent.

1. **StericX claim / audited proposition:** B5 is the maximum transverse atomic-sphere extent. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact convention expected:** max_i(||p_i||+rho_i), Å.
4. **Implementation location:** `src/geometry/sterimol.rs::params_from_projection`
5. **Independent validation method:** Independent analytic radial support; compare Morfeus sampled maximum separately.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Max bond-axis base error 5.23377122085833e-7 Å; Morfeus B5 itself uses an approximate angular scan.
7. **Confidence:** High
8. **Limitations:** The base maximum is not a global bound: a real near-antiparallel Kraken case has B5 error 0.00276961144 Å (G39). Extreme coordinates also require separate limits.

### G04 — Coordination-axis L applies the historical +0.40 Å convention.

1. **StericX claim / audited proposition:** Coordination-axis L applies the historical +0.40 Å convention. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html ; https://aarontools.readthedocs.io/en/latest/api/substituent.html
3. **Exact convention expected:** Raw API L plus 0.40 Å in CLI coordination mode; bond mode raw.
4. **Implementation location:** `src/descriptors.rs::sterimol_for_conformer; src/cli.rs::STERIMOL_L_CORRECTION`
5. **Independent validation method:** Read frozen source and compare raw API / Morfeus L_value_uncorrected and L_value.
6. **Result:** **SUPPORTED** — Arithmetic/convention placement agrees; raw API and CLI use distinct conventions.
7. **Confidence:** High implementation; moderate historical attribution
8. **Limitations:** Original Verloop chapter not fully accessible. The correction is a historical parameter convention, not a universal physical requirement for arbitrary metal dummies.

### G05 — Explicit hydrogen spheres contribute to Sterimol and real dummy is excluded.

1. **StericX claim / audited proposition:** Explicit hydrogen spheres contribute to Sterimol and real dummy is excluded. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact convention expected:** All non-dummy atoms with their chosen radii; virtual dummy includes real donor.
4. **Implementation location:** `src/geometry/sterimol.rs::compute/compute_with_dummy`
5. **Independent validation method:** Independent explicit-H graph molecules and envelope calculation.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Ordinary geometries agree within reported B1 scan/f32 limits.
7. **Confidence:** High
8. **Limitations:** Kraken DFT code uses H radius1.09 Å; StericX uses1.20 Å, so published reproduction is separate.

### G06 — Default bond-axis descriptors are independent of atom ordering.

1. **StericX claim / audited proposition:** Default bond-axis descriptors are independent of atom ordering. 
2. **Primary/reference source:** Geometric invariance requirement; https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact convention expected:** Same chemical structure and chosen chemical axis should produce the same descriptor.
4. **Implementation location:** `src/descriptors.rs::detect_donor; src/geometry/buried_volume.rs::bonded_neighbors`
5. **Independent validation method:** Frozen actual CLI on all six permutations of three equally distant chemically distinct substituents.
6. **Result:** **INCORRECT** — L changes6.0699997→3.5 Å and B5 changes3.5→6.0699997 Å because the tie selects the first atom index.
7. **Confidence:** High
8. **Limitations:** The underlying explicitly selected-axis kernel is invariant to reordering within f32 limits. This falsifies an unrestricted auto-axis invariance claim, not the mathematics for a specified axis.

### G07 — Sterimol is rigid-rotation/translation invariant.

1. **StericX claim / audited proposition:** Sterimol is rigid-rotation/translation invariant. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html
3. **Exact convention expected:** Physical continuous descriptors invariant for fixed chemical axis.
4. **Implementation location:** `src/geometry/sterimol.rs::compute/compute_with_dummy`
5. **Independent validation method:** 6 representatives×100 seeded random rotations plus X/Y/Z rotations and random translations.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — B1 max deviation0.019944190979003906 Å; L≤1.0013580322265625e-5 Å, B5≤4.76837158203125e-6 Å.
7. **Confidence:** High for measured range
8. **Limitations:** Rotations include modest translations ≤50 Å; no universal bound. B1 scan phase depends on alignment. A targeted exact 90-degree rotation of real near-axis Kraken cases changes L by 0.00172233582 Å and B5 by 0.00276851654 Å (G39), despite exact independent invariance.

### G08 — Buried volume integrates the union of atomic spheres inside a sphere.

1. **StericX claim / audited proposition:** Buried volume integrates the union of atomic spheres inside a sphere. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html ; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** V=integral of union membership inside R; %Vbur=100V/(4πR³/3).
4. **Implementation location:** `src/geometry/buried_volume.rs::occupied_volumes`
5. **Independent validation method:** Independent float64 union-of-balls integration and Morfeus under matched center/radii/grid convention.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — 23 finite base values max total-volume difference3.6215270000639066e-6 Å³; default lattice includes15,408 points.
7. **Confidence:** High for compatible default cases
8. **Limitations:** Numerical quadrature error is larger than implementation residual; invalid symmetry rejection and unknown radii are distinct failures.

### G09 — The default sphere radius3.5 Å follows a published descriptor convention.

1. **StericX claim / audited proposition:** The default sphere radius3.5 Å follows a published descriptor convention. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html
3. **Exact convention expected:** Sphere centered on selected real/virtual coordination point, R3.5 Å.
4. **Implementation location:** `src/geometry/buried_volume.rs::BuriedVolumeConfig; src/cli.rs`
5. **Independent validation method:** Official SambVca manual and Morfeus documentation.
6. **Result:** **VERIFIED** — Convention agrees.
7. **Confidence:** High
8. **Limitations:** This empirically useful convention does not prove a universal physically optimal sphere radius.

### G10 — Bondi radius scaling1.17 follows the reference convention.

1. **StericX claim / audited proposition:** Bondi radius scaling1.17 follows the reference convention. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html ; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** Every included vdW radius multiplied by1.17 before intersection.
4. **Implementation location:** `src/geometry/buried_volume.rs::aligned_atoms`
5. **Independent validation method:** Official manual plus independent scaled-radius calculations.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Scaling agrees on supported atoms and finite ordinary values.
7. **Confidence:** High
8. **Limitations:** Input radii table and overflow behavior must be validated separately.

### G11 — Hydrogens are excluded from buried-volume occupancy by default.

1. **StericX claim / audited proposition:** Hydrogens are excluded from buried-volume occupancy by default. 
2. **Primary/reference source:** https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html ; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** H may define connectivity/frame but contributes no sphere unless include_hydrogens=true.
4. **Implementation location:** `src/geometry/buried_volume.rs::aligned_atoms`
5. **Independent validation method:** Explicit-H methylphosphine/trimethylamine protocol pairs and PH3 analytic lens witness.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Exclusion/inclusion matches controlled Morfeus reference.
7. **Confidence:** High
8. **Limitations:** Removing donor-bound H from the frame would be a different, incorrect chemical operation.

### G12 — All supported radii are exactly Morfeus Bondi radii.

1. **StericX claim / audited proposition:** All supported radii are exactly Morfeus Bondi radii. 
2. **Primary/reference source:** https://doi.org/10.1021/j100785a001 ; frozen/reference_sources/morfeus/data.py
3. **Exact convention expected:** Declared table equality, not only the label Bondi-style.
4. **Implementation location:** `src/geometry/xyz.rs::van_der_waals_radius`
5. **Independent validation method:** Read independent table; per-element synthetic witnesses.
6. **Result:** **INCORRECT** — Br is1.85 Å in SUT vs1.83 in Morfeus; B1.92 in SUT vs Morfeus fallback2.0.
7. **Confidence:** High table difference
8. **Limitations:** The original Bondi table was not available for direct inspection, so the Br discrepancy was not adjudicated against that table; B has no Morfeus Bondi entry; Mantina et al., DOI 10.1021/jp8111556, explicitly gives the later 1.92 Å extension. Difference alone does not identify physical truth.

### G13 — Unknown-element radius fallbacks are chemically conservative.

1. **StericX claim / audited proposition:** Unknown-element radius fallbacks are chemically conservative. 
2. **Primary/reference source:** https://doi.org/10.1021/j100785a001
3. **Exact convention expected:** An arbitrary fallback must not be equated with an element-specific physical radius.
4. **Implementation location:** `src/geometry/xyz.rs::van_der_waals_radius/covalent_radius`
5. **Independent validation method:** He/Xe/Na/Se/Fe toy witnesses; official/reference tables.
6. **Result:** **UNCERTAIN** — SUT1.8 Å can be both too large(He) and too small(Na/Xe); max observed BV difference6.072736800747158 Å³,3.3813625981379403 percentage points.
7. **Confidence:** High numerical result; low chemical justification
8. **Limitations:** Morfeus also falls back2.0 for absent B/Fe entries; it is not ground truth. Broad chemical safety of either fallback lacks evidence.

### G14 — Audit proposition: the filled integration grid exactly reproduces Morfeus for every accepted density.

1. **StericX claim / audited proposition:** Audit proposition: the filled integration grid exactly reproduces Morfeus for every accepted density. Exact-compatibility audit question, not a quoted universal promise in the source. The documented SUT grid is evaluated against the independent reference grid.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** Float64 linspace and norm≤R in Morfeus versus f32 SUT arithmetic; grid endpoints included.
4. **Implementation location:** `src/geometry/buried_volume.rs::integration_grid`
5. **Independent validation method:** Four density levels plus odd-n grid and representable-float sphere-boundary perturbations.
6. **Result:** **INCORRECT** — Default n32 counts/membership agree for tested cases; odd grids reveal region-boundary differences, fine grids can differ at sphere lattice boundaries.
7. **Confidence:** High
8. **Limitations:** This rejects unrestricted exact all-density equivalence, not the shared continuum target or measured default-grid agreement. Floating-point sphere-boundary membership and region conventions differ; see convergence and partition records.

### G15 — Voxel atomic occupancy includes the vdW boundary.

1. **StericX claim / audited proposition:** Voxel atomic occupancy includes the vdW boundary. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html ; frozen Morfeus cKDTree query_ball_point implementation
3. **Exact convention expected:** distance²≤scaled_radius², corresponding to a closed ball.
4. **Implementation location:** `src/geometry/buried_volume.rs::occupied_volumes`
5. **Independent validation method:** Frozen one-point tests at1f32 and nextafter inside/outside, independent exact-real comparison.
6. **Result:** **VERIFIED** — Inside/on occupied; next representable outside unoccupied.
7. **Confidence:** High for tested boundary
8. **Limitations:** Open/closed spheres have equal continuum volume; this establishes a lattice convention, not a uniquely physical inequality.

### G16 — Sphere boundary handling includes points on R.

1. **StericX claim / audited proposition:** Sphere boundary handling includes points on R. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** Retain ||point||≤R.
4. **Implementation location:** `src/geometry/buried_volume.rs::integration_grid`
5. **Independent validation method:** Private-grid observation at n5 with R and adjacent representable radii.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Boundary points are present; outside lattice points rejected.
7. **Confidence:** High
8. **Limitations:** f32 rounding can move mathematically marginal points. Atom centers outside R must still contribute when their sphere intersects R, confirmed by boundary-atom witnesses.

### G17 — Audit proposition: quadrants use exactly the same finite-grid regional convention as Morfeus at every accepted density.

1. **StericX claim / audited proposition:** Audit proposition: quadrants use exactly the same finite-grid regional convention as Morfeus at every accepted density. Reference-convention compatibility question, not an assertion that an open or closed continuum boundary is physically preferable.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** Canonical xy signs (++),(-+),(--),(+-); Morfeus excludes zero planes.
4. **Implementation location:** `src/geometry/buried_volume.rs::quadrant_index/occupied_volumes`
5. **Independent validation method:** Every individual quadrant in every orientation, independent/Morfeus comparison.
6. **Result:** **INCORRECT** — Default even-grid base values agree within f32 rounding. n5 qvbur_min differs9.621127681084193 Å³; odd-grid assignment differs structurally.
7. **Confidence:** High
8. **Limitations:** Default even-grid agreement remains verified within the reported numerical limits. At odd grids SUT assigns zeros positive and directly normalizes quadrants; Morfeus excludes zero planes and sums octants. Both approximate the same continuum regions without exact numerical convention equivalence.

### G18 — Audit proposition: octants use exactly the same finite-grid regional convention as Morfeus at every accepted density.

1. **StericX claim / audited proposition:** Audit proposition: octants use exactly the same finite-grid regional convention as Morfeus at every accepted density. Reference-convention compatibility question; it does not identify either implementation as physical truth.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** Four xy quadrants for positive and negative z; explicit Morfeus ID remapping.
4. **Implementation location:** `src/geometry/buried_volume.rs::octant_index/occupied_volumes`
5. **Independent validation method:** Every individual octant per orientation compared, canonical signs preserved.
6. **Result:** **INCORRECT** — Even-grid agreement is strong; zero-plane regions at odd n differ.
7. **Confidence:** High
8. **Limitations:** Default even-grid numerical agreement remains favorable. Odd-grid zero-plane assignment differs. Negative-z numeric IDs also differ, but the harness correctly remaps them; this classification concerns the actual boundary convention, not a mistaken ID comparison.

### G19 — Near and far describe donor-facing and distal hemispheres.

1. **StericX claim / audited proposition:** Near and far describe donor-facing and distal hemispheres. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html ; Kraken official SI page12
3. **Exact convention expected:** Center→donor is negative z; near z<0, far z≥0 in SUT.
4. **Implementation location:** `src/geometry/buried_volume.rs::coordinate_basis/occupied_volumes`
5. **Independent validation method:** Independent coordinate frame plus per-octant sums.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Hemisphere naming/orientation agrees under the stated coordinate convention.
7. **Confidence:** High
8. **Limitations:** This is not the outside-sphere distal molecular volume. Odd-grid normalization creates partition inconsistency.

### G20 — Audit proposition: each finite-grid regional partition adds exactly to the independently estimated total for all valid configurations.

1. **StericX claim / audited proposition:** Audit proposition: each finite-grid regional partition adds exactly to the independently estimated total for all valid configurations. Numerical conservation audit question, distinct from exact additivity of the continuum integral; not a quoted all-grid guarantee.
2. **Primary/reference source:** Continuum additivity of disjoint integrals; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** Regional volume sum should converge to total, with quantified finite-grid discrepancy.
4. **Implementation location:** `src/geometry/buried_volume.rs::occupied_volumes`
5. **Independent validation method:** Independent conservation check across densities.
6. **Result:** **INCORRECT** — Triphenylphosphine coarse total53.1568947 versus near+far59.7332687 Å³ (+12.37%); default even-grid error≤rounding.
7. **Confidence:** High
8. **Limitations:** Default even-grid additivity holds within rounding in the campaign. The odd-grid discrepancy is a finite-estimator limitation, not a failure of continuum volume additivity. At density 0.0001 near+far remains 0.5790825 Å³ above total. Morfeus also has a coarse-grid gap of +3.4082601 Å³.

### G21 — Maximum adjacent-quadrant difference is reduced correctly.

1. **StericX claim / audited proposition:** Maximum adjacent-quadrant difference is reduced correctly. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html ; official Kraken DFT code
3. **Exact convention expected:** max over three plane orientations and four cyclic adjacent pairs of absolute difference.
4. **Implementation location:** `src/geometry/buried_volume.rs::compute_from_center`
5. **Independent validation method:** Independent explicit twelve-quadrant reduction.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Default matched cases agree within f32 subtraction rounding.
7. **Confidence:** High
8. **Limitations:** Frame/center selection, zero-plane treatment and symmetry rejection are separate issues.

### G22 — A positive buried volume with zero quadrant asymmetry is unphysical/degenerate.

1. **StericX claim / audited proposition:** A positive buried volume with zero quadrant asymmetry is unphysical/degenerate. 
2. **Primary/reference source:** Sphere intersection symmetry; independent analytic lens; https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** An axisymmetric nonzero occupied set has equal quadrants and valid zero asymmetry.
4. **Implementation location:** `src/geometry/buried_volume.rs::compute_from_center final rejection`
5. **Independent validation method:** PH3 with three chemically bonded H; full-cover toy; small sphere wholly inside ordinary donor radius.
6. **Result:** **INCORRECT** — PH3 correctly has31.937214517315244 Å³ lattice volume,17.782969885773625%,zero asymmetry; SUT errors despite valid frame. Full-cover100% also errors.
7. **Confidence:** High, constructive counterexamples
8. **Limitations:** PH3 CLI separately rejects lack of heavy reference atom; public BV API accepts H as reference then hits invalid symmetry rejection. Small-sphere case also reproduces through CLI with config.

### G23 — Audit proposition: default buried-volume precision supports scientific distinctions finer than its discretization error.

1. **StericX claim / audited proposition:** Audit proposition: default buried-volume precision supports scientific distinctions finer than its discretization error. Scientific-resolution audit question; no explicit application-specific resolution guarantee was identified in the cited kernel. The measured numerical floor does not itself establish a useful chemical effect size.
2. **Primary/reference source:** Analytic two-sphere intersection formula
3. **Exact convention expected:** Compare lattice to continuum, not only another lattice implementation.
4. **Implementation location:** `src/geometry/buried_volume.rs::integration_grid/occupied_volumes`
5. **Independent validation method:** Exact PH3 one-atom lens and four densities.
6. **Result:** **UNCERTAIN** — Analytic31.669141713568127 Å³; default error+0.2680721841369511 Å³ or0.149265 percentage points; density.001 error−0.0284346976; .0001 error+0.0006218515 Å³.
7. **Confidence:** High for the measured lens/convergence errors; insufficient evidence for scientific resolution below them
8. **Limitations:** The analytic lens has default error 0.2680722 Å³, and four real/toy molecules differ from very fine totals by up to 0.1971054 Å³. Error cancellation for a particular descriptor difference, experimental relevance, and a universal resolution guarantee are not established. Agreement with another finite grid cannot supply that missing evidence.

### G24 — Buried-volume values are rigid-transformation invariant over the tested moderate-coordinate campaign.

1. **StericX claim / audited proposition:** Buried-volume values are rigid-transformation invariant over the tested moderate-coordinate campaign. Scoped numerical invariance proposition: six prespecified representatives and 618 transforms, not a universal finite-precision guarantee.
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html
3. **Exact convention expected:** Co-rotating molecular frame and center; finite-grid sensitivity quantified.
4. **Implementation location:** `src/geometry/buried_volume.rs::coordinate_basis/occupied_volumes`
5. **Independent validation method:** 618 transformed cases across six representatives, all public fields and raw region values.
6. **Result:** **SUPPORTED** — All nine public BV fields had exactly0 measured deviation on this campaign.
7. **Confidence:** High for campaign; limited generality
8. **Limitations:** Grid transitions can occur for other geometries; huge translations lose f32 geometric information and are not covered by this favorable result.

### G25 — Three covalently bonded substituents include donor-bound hydrogen.

1. **StericX claim / audited proposition:** Three covalently bonded substituents include donor-bound hydrogen. 
2. **Primary/reference source:** https://doi.org/10.1039/B801115J ; explicit molecular graphs
3. **Exact convention expected:** SMILES graph bonds define chemically intended neighbors; H included.
4. **Implementation location:** `src/geometry/buried_volume.rs::bonded_neighbors/donor_neighbor_indices`
5. **Independent validation method:** 12 independent RDKit graph molecules plus primary/secondary phosphine contact witnesses.
6. **Result:** **SUPPORTED** — All12 standard graphs match; both H-bearing historical witnesses recover chemically bonded neighbors.
7. **Confidence:** High for test molecules
8. **Limitations:** Geometry-only inference is not universal chemistry; charge/bond order/coordination identity are ignored.

### G26 — Cordero covalent radii justify a universal1.3× bond cutoff.

1. **StericX claim / audited proposition:** Cordero covalent radii justify a universal1.3× bond cutoff. 
2. **Primary/reference source:** https://doi.org/10.1039/B801115J
3. **Exact convention expected:** Cordero provides empirical atomic radii, not universal binary bonding truth at1.3× sums.
4. **Implementation location:** `src/geometry/buried_volume.rs::BOND_TOLERANCE_FACTOR/bonded_neighbors`
5. **Independent validation method:** Primary source interpretation, threshold nextafter witnesses, independent graphs.
6. **Result:** **UNCERTAIN** — Threshold behavior is deterministic and matches stated heuristic; scientific universality unsupported. Coincident atoms count as bonded before later geometry validation.
7. **Confidence:** High implementation; insufficient chemical universality
8. **Limitations:** No threshold retuning was performed. Stretching a synthetic bond beyond cutoff is not itself proof the chemistry is misclassified.

### G27 — Automatic donor selection identifies the expected donors in the tested standard P/N graphs.

1. **StericX claim / audited proposition:** Automatic donor selection identifies the expected donors in the tested standard P/N graphs. Scoped test of automatic-selection behavior against independently supplied molecular graphs; general chemical donor availability remains unestablished.
2. **Primary/reference source:** https://doi.org/10.1039/B801115J ; CLI contract
3. **Exact convention expected:** Sole chosen element, or explicit index, plus exactly3 inferred neighbors and≥1heavy neighbor.
4. **Implementation location:** `src/descriptors.rs::detect_donor`
5. **Independent validation method:** Frozen CLI runs with phosphines/amines, multiple P, planar/coincident cases and PH3.
6. **Result:** **SUPPORTED** — Standard tested P/N donors work; multiple-element ambiguity is rejected; PH3 rejected for lack of heavy reference; chemistry beyond geometry unestablished.
7. **Confidence:** Moderate
8. **Limitations:** Explicit index can select any element; trivalent geometry does not establish lone-pair availability, charge state, aromaticity or catalytic coordination.

### G28 — Current frame fixes the historical nearest-heavy error.

1. **StericX claim / audited proposition:** Current frame fixes the historical nearest-heavy error. 
2. **Primary/reference source:** Independent known SMILES connectivity; primary/secondary phosphine valence
3. **Exact convention expected:** Three bonded neighbors, including H, rather than three nearest heavy atoms.
4. **Implementation location:** `src/geometry/buried_volume.rs::donor_neighbor_indices`
5. **Independent validation method:** Reconstruct former selector independently; controls plus H and nearby nonbonded heavy atoms.
6. **Result:** **SUPPORTED** — Tertiary control center unchanged; primary/secondary wrong-frame centers displaced2.4302383/2.0563763 Å and BV descriptors differ strongly; current neighborhoods match graphs.
7. **Confidence:** High for reconstructed failure mechanism
8. **Limitations:** The exact six historically zero ligands (575,1485,1487,1490,1491,1495),20 conformers, were then independently reconstructed from primary SDF bond graphs; every old nearest-heavy frame again gives zero while current matched-frame values agree with independent references. Kraken reproduction conventions remain separate.

### G29 — The geometric virtual-metal placement exactly matches published Kraken.

1. **StericX claim / audited proposition:** The geometric virtual-metal placement exactly matches published Kraken. 
2. **Primary/reference source:** Official Kraken SI FigureS2/page10 and frozen DFT source audited in kraken/
3. **Exact convention expected:** Kraken DFT sums raw substituent displacements then normalizes; SUT sums unit bond directions.
4. **Implementation location:** `src/geometry/buried_volume.rs::infer_lone_pair_direction`
5. **Independent validation method:** Primary code/SI comparison, independent center construction.
6. **Result:** **INCORRECT** — Conventions differ for unequal bond lengths; matching center-distance2.28 alone is insufficient.
7. **Confidence:** High
8. **Limitations:** Separate Kraken audit quantifies published-data consequences. Geometric placement is itself an approximation to actual coordination geometry.

### G30 — Planar fallback is a physically determined lone-pair direction.

1. **StericX claim / audited proposition:** Planar fallback is a physically determined lone-pair direction. 
2. **Primary/reference source:** Geometric symmetry; no universal primary rule identified
3. **Exact convention expected:** Planar center has no unique opposite-bond-sum direction; plane normal sign is ambiguous without environment/electronics.
4. **Implementation location:** `src/geometry/buried_volume.rs::infer_lone_pair_direction/center_clearance`
5. **Independent validation method:** Planar and nearly planar synthetic structures, rotations and permutations.
6. **Result:** **UNCERTAIN** — Fallback normal/maximum-clearance criterion is a StericX design choice; independent naive normalized-sum centers become unstable near exact planarity.
7. **Confidence:** High description; insufficient physical evidence
8. **Limitations:** Controlled-center kernel comparisons deliberately do not treat SUT-selected center as scientifically established.

### G31 — Pyramidalization P implements the determinant and acute-correction convention used by Morfeus 0.8.0.

1. **StericX claim / audited proposition:** Pyramidalization P implements the determinant and acute-correction convention used by Morfeus 0.8.0. Runnable-reference implementation proposition. Morfeus attributes the measure to Radhakrishnan and Agranat; original-paper textual equivalence was not established because full-text access was unavailable.
2. **Primary/reference source:** https://doi.org/10.1007/BF00676621 ; https://digital-chemistry-laboratory.github.io/morfeus/api/morfeus.pyramidalization.html
3. **Exact convention expected:** |det(unit donor vectors)| with2−P acute correction.
4. **Implementation location:** `src/geometry/pyramidalization.rs::compute`
5. **Independent validation method:** Independent determinant calculation, Morfeus, tetrahedral/orthogonal/acute witnesses.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Ordinary graph cases maxP error1.2938381854787906e-7; tetrahedral0.769800305 versus analytic0.7698003589.
7. **Confidence:** High for runnable reference; moderate original-paper definition
8. **Limitations:** Original paper full text not acquired. Degenerate sentinel behavior is not a valid P measurement.

### G32 — Pyramidalization alpha uses the mean signed plane-normal angle in degrees.

1. **StericX claim / audited proposition:** Pyramidalization alpha uses the mean signed plane-normal angle in degrees. 
2. **Primary/reference source:** https://doi.org/10.1007/BF00676621 ; Morfeus0.8 source
3. **Exact convention expected:** Three alpha values, acute bisector sign, degree conversion; ordinary planar90°, orthogonal0°.
4. **Implementation location:** `src/geometry/pyramidalization.rs::compute`
5. **Independent validation method:** Independent plane geometry and explicit acute/planar/tetrahedral cases.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Ordinary graph maxerror6.122375019401716e-6°; acute example−82.717° is permitted.
7. **Confidence:** High in nondegenerate range
8. **Limitations:** Near-collinear geometry is assigned0 by SUT despite finite independent29.99992497364023°; sentinel must not be treated scientifically.

### G33 — P=1 corresponds to an ideal tetrahedral apex.

1. **StericX claim / audited proposition:** P=1 corresponds to an ideal tetrahedral apex. 
2. **Primary/reference source:** Analytic tetrahedral unit vectors; https://doi.org/10.1007/BF00676621
3. **Exact convention expected:** Three tetrahedral directions determinant4/(3√3)≈0.76980036; orthogonal triad1.
4. **Implementation location:** `src/geometry/pyramidalization.rs::PyramidalizationParams documentation`
5. **Independent validation method:** Analytic determinant and frozen directAPI.
6. **Result:** **TERMINOLOGY ISSUE** — The comment saying1 for ideal tetrahedral-like apex is misleading; tests compute tetrahedral0.7698.
7. **Confidence:** High
8. **Limitations:** Tetrahedral-like is imprecise; no universal physical normalization to1 for tetrahedral geometry is established.

### G34 — Pyramidalization is neighbor-order invariant.

1. **StericX claim / audited proposition:** Pyramidalization is neighbor-order invariant. 
2. **Primary/reference source:** https://doi.org/10.1007/BF00676621
3. **Exact convention expected:** Geometric three-vector measure independent of permutation.
4. **Implementation location:** `src/geometry/pyramidalization.rs::compute`
5. **Independent validation method:** 120 molecular permutations with mapped neighbor identities.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — No scientifically meaningful order variation observed; full residuals in invariance.csv.
7. **Confidence:** High for tested domain
8. **Limitations:** Planar/degenerate numerical branches need explicit handling.

### G35 — Invalid or nearly degenerate geometry produces an explicit scientific failure.

1. **StericX claim / audited proposition:** Invalid or nearly degenerate geometry produces an explicit scientific failure. 
2. **Primary/reference source:** Mathematical nondefinition of directions/plane normals
3. **Exact convention expected:** Do not silently interpret sentinel zero as a measured descriptor.
4. **Implementation location:** `src/geometry/sterimol.rs::compute; src/geometry/pyramidalization.rs::compute`
5. **Independent validation method:** Coincident, short-axis, nearly collinear, NaN/Inf input witnesses.
6. **Result:** **INCORRECT** — Sterimol and pyramidalization direct APIs returnzero sentinels; BV errors. Near-collinear alpha29.999925° becomes0°; short nonzeroaxis givesallzeroSterimol.
7. **Confidence:** High for unrestricted explicit-error claim
8. **Limitations:** Zero sentinel behavior is documented at kernel level. CLI combined-descriptor validation often rejects these cases; callers of public kernels remain exposed.

### G36 — Finite input coordinates always produce finite scientific outputs.

1. **StericX claim / audited proposition:** Finite input coordinates always produce finite scientific outputs. 
2. **Primary/reference source:** IEEE754 finite arithmetic range; analytic norm
3. **Exact convention expected:** Overflow must be detected if finite output is promised.
4. **Implementation location:** `src/geometry/sterimol.rs::params_from_projection`
5. **Independent validation method:** Frozen finite1e38 Å remote atom witness.
6. **Result:** **INCORRECT** — B5 becomes+Inf for wholly finite positions/radii; observer recordsnull plus explicitnonfinite map.
7. **Confidence:** High
8. **Limitations:** Deliberately extreme and not ordinary chemistry. It identifies an API validation limit, not evidence that f64 is needed for typical coordinates.

### G37 — f32 is sufficiently accurate for ordinary geometry in this audit.

1. **StericX claim / audited proposition:** f32 is sufficiently accurate for ordinary geometry in this audit. 
2. **Primary/reference source:** Independent f64 equations and Morfeus under controlled protocol
3. **Exact convention expected:** Separate rounding error from grid/angle discretization and coordinate convention.
4. **Implementation location:** `src/geometry/*.rs`
5. **Independent validation method:** Graph molecules, rigid campaign, f32-rounded input re-evaluated in f64.
6. **Result:** **UNCERTAIN** — Restricted ordinary-case L/B5 errors are mostly 1e-6–1e-5 Å and default BV residuals a few 1e-6 Å³. The targeted real near-Z cases instead have ~0.002 Å alignment errors. Scientific sufficiency requires an application-specific tolerance or effect-size assessment.
7. **Confidence:** High measured errors; insufficient evidence for universal scientific sufficiency
8. **Limitations:** Real near-axis Kraken geometries have larger alignment errors (G39), beyond ordinary rounding. Large translations 1e7–1e20 destroy coordinate detail; nearly degenerate axes amplify rounding. No blanket recommendation to switch f64 is made.

### G38 — Geometry descriptor agreement validates predictive/experimental chemistry.

1. **StericX claim / audited proposition:** Geometry descriptor agreement validates predictive/experimental chemistry. 
2. **Primary/reference source:** Scientific identifiability; held-out experimental validation required
3. **Exact convention expected:** Numerical descriptor compatibility and chemistry prediction are different hypotheses.
4. **Implementation location:** `README scientific summaries; studies002/004/005`
5. **Independent validation method:** Separate scope of independent computation from model/generalization evidence.
6. **Result:** **OUT OF SCOPE** — This geometry audit establishes no experimental predictive validity.
7. **Confidence:** High
8. **Limitations:** Requires held-out chemical/experimental evidence evaluated in model audit and prospective experiments.

### G39 — B5 is the exact radial support about the supplied attachment axis, apart from ordinary f32 arithmetic rounding.

1. **StericX claim / audited proposition:** B5 is the exact radial support about the supplied attachment axis, apart from ordinary f32 arithmetic rounding. 
2. **Primary/reference source:** https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html ; https://docs.rs/glam/0.30.10/src/glam/f32/sse2/quat.rs.html
3. **Exact convention expected:** Projection onto the supplied unit axis must remain invariant under a rigid change of global Cartesian coordinates; a near-Z axis must not be replaced by Z.
4. **Implementation location:** `src/geometry/sterimol.rs::compute/compute_with_dummy, Quat::from_rotation_arc in frozen glam 0.30.10`
5. **Independent validation method:** Independent float64 dot/norm projections on identical SUT-rounded coordinates/radii/center; two real Kraken conformers minimized to donor plus active atom; exact 90-degree X coordinate permutation.
6. **Result:** **INCORRECT** — Full-case L is 6.859999656677246 versus independent 6.858277336508549 Å (963/48394). Full-case B5 is 7.116985321044922 versus independent 7.114215709601416 Å (1075/49864). Exact coordinate rotation reduces corresponding errors to 1.56e-8 and 1.09e-6 Å; original minimized witnesses retain the failure.
7. **Confidence:** High, preserved full and minimized cases plus targeted invariance experiment
8. **Limitations:** The supplied axes are mathematically well defined and the source geometries are real DFT structures. glam documents approximate near-singular alignment (|dot| > 1−2 f32 EPSILON); this is a StericX accuracy-contract issue, not an undocumented glam defect. It is distinct from B1 angular discretization and does not prove a chemistry prediction is materially affected. See geometry/alignment/ and its rotation/ subdirectory.

## Conformers and kinetics

### C01 — Eyring rate arithmetic implements the activation-free-energy equation.

1. **StericX claim / audited proposition:** Eyring rate arithmetic implements the activation-free-energy equation. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** k=κkBT/h exp(-ΔG‡/RT); κ=1; absolute ΔG in kcal/mol, T inK; first-order units.
4. **Implementation location:** `src/kinetics/eyring.rs::calculate_rate_constant`
5. **Independent validation method:** 75-digit Decimal equation on frozen f32 inputs over temperature/barrier grid.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Normal-rate subset maximum relative error 6.9231e-8; subnormal and extreme failures separately retained.
7. **Confidence:** High
8. **Limitations:** Equation correctness does not identify a physical barrier or establish transition-state assumptions.

### C02 — Physical constants and kcal/Hartree conversion are accurate.

1. **StericX claim / audited proposition:** Physical constants and kcal/Hartree conversion are accurate. 
2. **Primary/reference source:** NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** kB/h/NA exact SI; R=kBNA/4184; Hartree energy fromCODATA.
4. **Implementation location:** `src/kinetics/eyring.rs constants; scripts/stericx_quantum.py constants`
5. **Independent validation method:** Independent Decimal arithmetic from primary constants.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — R truncation relative difference ≈3.22e-10; Hartree conversion difference ≈2.03e-10 kcal mol^-1 Hartree^-1.
7. **Confidence:** High
8. **Limitations:** Numerical effects generally smaller than f32 rounding; physical energy-model errors not tested by constants.

### C03 — Zero ΔΔG gives equal R:S and reversing sign exchanges preference.

1. **StericX claim / audited proposition:** Zero ΔΔG gives equal R:S and reversing sign exchanges preference. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** ΔΔG=G‡S−G‡R, common reactant and equal pathway prefactors; zero→50:50.
4. **Implementation location:** `src/kinetics/eyring.rs::product_ratio`
5. **Independent validation method:** Independent two-pathway rates with common energy shift, signed grid and signed zero.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — R/S max error 3.7596e-6 percentage points; zero and sign symmetry observed.
7. **Confidence:** High
8. **Limitations:** This algebraic convention cannot infer the absolute enantiomer from unsigned experimental ee.

### C04 — Enantiomeric excess arithmetic follows the two-pathway populations.

1. **StericX claim / audited proposition:** Enantiomeric excess arithmetic follows the two-pathway populations. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Absolute ee=100|pR−pS|; signed ee, if requested, must retain configuration convention.
4. **Implementation location:** `src/kinetics/eyring.rs::calculate_enantiomeric_excess`
5. **Independent validation method:** Independent pathway probabilities, zero and small/large differences.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Max absolute error 7.5191e-6 percentage points; tiny ee and minority populations can round to zero.
7. **Confidence:** High
8. **Limitations:** Finite precision saturation is not experimentally complete selectivity.

### C05 — simulate can output an absolute rate from the same ΔΔG used for selectivity.

1. **StericX claim / audited proposition:** simulate can output an absolute rate from the same ΔΔG used for selectivity. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** A difference fixes a rate ratio, not either absolute activation barrier.
4. **Implementation location:** `src/commands/simulate.rs; README simulate example`
5. **Independent validation method:** Two barrier pairs(10,11) and(20,21) have same difference but different rates; frozen CLI observation.
6. **Result:** **INCORRECT** — Same ratio 5.4076065654; R rates 290543.5244 versus 0.01358815005 s^-1.
7. **Confidence:** High
8. **Limitations:** The standalone rate function is correct when supplied an actual absolute barrier; the interface is the problem.

### C06 — Rate and selectivity calculations have a defined finite numerical domain.

1. **StericX claim / audited proposition:** Rate and selectivity calculations have a defined finite numerical domain. Audit of the accepted numerical domain, not an assertion that documentation promises finite outputs for every finite input.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; IEEE754 f32 range
3. **Exact convention expected:** Finite input can map outside representable output range; this must be reported.
4. **Implementation location:** `src/kinetics/eyring.rs; src/commands/simulate.rs`
5. **Independent validation method:** Large positive/negative barriers, temperature extremes, NaN/infinity and zero-temperature inputs.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Wider grid includes 16 rate underflows and 15 nonfinite rates; CLI checks inputs only.
7. **Confidence:** High
8. **Limitations:** This is a bounded precision/validation claim, not evidence that ordinary barriers have chemically material errors.

### C07 — MMFF ensemble weights use the normalized Boltzmann equation.

1. **StericX claim / audited proposition:** MMFF ensemble weights use the normalized Boltzmann equation. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** wi∝exp(-(Ei−Emin)/(RT)); kcal/mol; positive temperature; unit degeneracy.
4. **Implementation location:** `scripts/prepare_data.py::embed_and_optimize`
5. **Independent validation method:** Controlled energy/status injection into frozen routine; independent Decimal weights for single/equal/two/three/offset/extreme/temperature/window cases.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum population error 1.9429e-16 in isolated weighting tests.
7. **Confidence:** High for arithmetic
8. **Limitations:** Mocked optimizer results do not validate search completeness or MMFF energies as thermodynamic free energies.

### C08 — Calling the retained MMFF ensemble fully optimized implies that its conformers converged.

1. **StericX claim / audited proposition:** Calling the retained MMFF ensemble fully optimized implies that its conformers converged. Audit of a possible stronger reading of “optimized”; the frozen code retains optimizer status and does not explicitly promise convergence of every retained state.
2. **Primary/reference source:** RDKit installed2026.3.4 MMFFOptimizeMoleculeConfs docstring and official https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html
3. **Exact convention expected:** Return status 0 means converged; finite nonzero status may still be unconverged.
4. **Implementation location:** `scripts/prepare_data.py::embed_and_optimize`
5. **Independent validation method:** Inject status −1/0/1 and nonfinite energies; inspect retained IDs/statuses.
6. **Result:** **INCORRECT** — Status 1 is retained and weighted; negative status and nonfinite energy are excluded.
7. **Confidence:** High
8. **Limitations:** Do not confuse retaining an unconverged state with a proven material descriptor error; magnitude depends on that geometry.

### C09 — CREST populations correspond to configured QuantumConfig.temperature_k.

1. **StericX claim / audited proposition:** CREST populations correspond to configured QuantumConfig.temperature_k. 
2. **Primary/reference source:** CREST official --temp documentation and frozen v2.12 cregen.f90: https://crest-lab.github.io/crest-docs/page/documentation/keywords.html; IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Population temperature must equal requested T; external table temperature cannot be assumed.
4. **Implementation location:** `scripts/stericx_quantum.py::crest_ensemble and _conformer_thermodynamics`
5. **Independent validation method:** Frozen realistic two-state summary at 298.15 K, requested 500 K; independently evaluate weights.
6. **Result:** **INCORRECT** — Population error 0.1116127958 at 500 K; command omits --temp and parsed table bypasses recomputation.
7. **Confidence:** High
8. **Limitations:** Direct frozen helper/command audit; no actual external CREST execution is claimed.

### C10 — CREST missing-table fallback correctly weights supplied electronic energies.

1. **StericX claim / audited proposition:** CREST missing-table fallback correctly weights supplied electronic energies. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** Convert Hartree differences to kcal/mol, subtract minimum, normalize atconfigured T.
4. **Implementation location:** `scripts/stericx_quantum.py::_conformer_thermodynamics`
5. **Independent validation method:** Frozen two-frame energy input at 298.15/500 K, independent Hartree conversion.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum population error≈1.03e-12, including coordinate-comment energy rounding.
7. **Confidence:** High for tested arithmetic
8. **Limitations:** Degeneracy defaults to one and electronic energies are not state free energies; missing/incorrect degeneracies change expected populations.

### C11 — Parsed CREST weights are valid nonnegative probabilities.

1. **StericX claim / audited proposition:** Parsed CREST weights are valid nonnegative probabilities. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Each population ≥0; positive finite normalization total.
4. **Implementation location:** `scripts/stericx_quantum.py::_conformer_thermodynamics`
5. **Independent validation method:** Synthetic malformed population table with weights −0.5,1.5.
6. **Result:** **INCORRECT** — Helper accepts negative population because only finiteness and sum>0 checked.
7. **Confidence:** High
8. **Limitations:** Malformed-table validation test; does not assert real CREST normally emits negative populations.

### C12 — Native supplied conformer weights correctly average Sterimol descriptors.

1. **StericX claim / audited proposition:** Native supplied conformer weights correctly average Sterimol descriptors. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Normalize supplied nonnegative weights then Σwi di; Å; preserve ensemble extrema/count.
4. **Implementation location:** `src/reaction.rs::{conformer_weights,record_from_ensemble}; src/commands/parse.rs`
5. **Independent validation method:** Native frozen CLI on an analytic XYZ pair; independent packed-record decoding.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Accepted cases max descriptor error 2.6703e-7 Å; 1:3 supplied weights give L=5.2, B1=1.7, B5=4.97.
7. **Confidence:** High for specified cases
8. **Limitations:** Does not establish the supplied weights are thermodynamically justified.

### C13 — Missing weights with supplied conformer energies produce Boltzmann populations.

1. **StericX claim / audited proposition:** Missing weights with supplied conformer energies produce Boltzmann populations. Stronger Boltzmann interpretation tested by the audit; the function itself explicitly implements a uniform fallback when weights are absent.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Energy-derived populations require an explicit equation and T; uniform averaging follows equal state weights.
4. **Implementation location:** `src/reaction.rs::conformer_weights`
5. **Independent validation method:** Native parse CSV with unequal energies 0;1 and missing weights.
6. **Result:** **INCORRECT** — Native path uses a uniform average, L=4.7 Å, irrespective of energies; differs from a Boltzmann expectation.
7. **Confidence:** High
8. **Limitations:** Uniform fallback is explicit implementation behavior; documentation must not describe all native averages as energy-derived.

### C14 — Missing or invalid native conformer geometry does not silently produce a partial ensemble.

1. **StericX claim / audited proposition:** Missing or invalid native conformer geometry does not silently produce a partial ensemble. Validation behavior tested by the audit, not a quoted promise of partial-ensemble renormalization.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; nativeingestionspecifiedsemantics
3. **Exact convention expected:** Reject an incomplete ensemble or explicitly track excluded population; avoid silent partial results.
4. **Implementation location:** `src/commands/parse.rs`
5. **Independent validation method:** A missing file and a malformed XYZ as the second conformer.
6. **Result:** **SUPPORTED** — Both return errors and no packed record; the native path does not silently discard the bad conformer.
7. **Confidence:** High for cases
8. **Limitations:** Python preparation intentionally filters some failed energies; these are distinct stage semantics.

### C15 — Buried-volume ensemble weightedmeans implement normalized sums.

1. **StericX claim / audited proposition:** Buried-volume ensemble weightedmeans implement normalized sums. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Use the same populations for each scalar; treat extrema and conformer-identity reductions separately.
4. **Implementation location:** `src/geometry/buried_volume.rs::aggregate`
5. **Independent validation method:** One/two/three/nonunit/permuted analytic descriptors versus Decimal.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum weighted-mean error 6.3578e-7 Å³.
7. **Confidence:** High for finite validatedinputs
8. **Limitations:** Tiny/huge weight normalization and NaN parameters are separate failures; no complete chemical ensemble claim.

### C16 — Positive finite weights can be normalized independent of common scale.

1. **StericX claim / audited proposition:** Positive finite weights can be normalized independent of common scale. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Multiplying all weightsby any positive finite factor should not change valid probabilities when representable.
4. **Implementation location:** `src/reaction.rs::conformer_weights; src/geometry/buried_volume.rs::aggregate`
5. **Independent validation method:** Compare weights 1:3 with 1e-10:3e-10; also test sum overflow.
6. **Result:** **INCORRECT** — A valid tiny-weight ratio is rejected as zero total; a large finite sum can overflow.
7. **Confidence:** High
8. **Limitations:** Ordinary normalized weights are unaffected; a correction should rescale before summation.

### C17 — Public buried-volume aggregation has a scientifically defined result for accepted inputs.

1. **StericX claim / audited proposition:** Public buried-volume aggregation has a scientifically defined result for accepted inputs. Public-API domain/validation test; no evidence that the ordinary validated CLI supplies NaN descriptors.
2. **Primary/reference source:** Weightedmean mathematicaldomain and publicAPI
3. **Exact convention expected:** Nonfinite components require an error or an explicit missing-data meaning.
4. **Implementation location:** `src/geometry/buried_volume.rs::aggregate`
5. **Independent validation method:** Direct frozen public API: NaN descriptor with weight 1.
6. **Result:** **INCORRECT** — Returns Ok with NaN vbur_boltz; exact nonfinite output preserved.
7. **Confidence:** High
8. **Limitations:** Direct API finding; normal CLI geometry validation may prevent this input.

### C18 — Stored conformer energy_span is maximum minus minimum.

1. **StericX claim / audited proposition:** Stored conformer energy_span is maximum minus minimum. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Energy span=max(E)−min(E); alternatively require relative energies to have minimum zero.
4. **Implementation location:** `src/reaction.rs::conformer_energy_span; PackedReactionRecord metadata`
5. **Independent validation method:** Accepted native CSV energies [1,2], independently decode the packed record.
6. **Result:** **INCORRECT** — Stored span is 2 instead of 1; code uses max and does not require minimum zero.
7. **Confidence:** High
8. **Limitations:** When relative energies are correctly referenced to the minimum, max equals max-minus-min.

### C19 — ee-to-ΔΔG conversion reproduces the standard magnitude equation.

1. **StericX claim / audited proposition:** ee-to-ΔΔG conversion reproduces the standard magnitude equation. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** |ΔΔG|=RTln((100+|ee|)/(100−|ee|)); |ee|<100.
4. **Implementation location:** `scripts/prepare_data.py::ee_to_ddg`
5. **Independent validation method:** Independent equation for validpositive/negative ee and two temperatures.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Valid interior inputs agree within floating-point precision.
7. **Confidence:** High
8. **Limitations:** Requires kinetic selectivity assumptions; absolute value intentionally discards enantiomer identity.

### C20 — 100% and out-of-range ee can be converted to a finite exact barrier.

1. **StericX claim / audited proposition:** 100% and out-of-range ee can be converted to a finite exact barrier. Mathematical interpretation of an implemented clipping behavior; no measured detection limit was supplied.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** 100% is an infinite mathematical limit; >100% is invalid; censored measurements require bounds.
4. **Implementation location:** `scripts/prepare_data.py::ee_to_ddg`
5. **Independent validation method:** Cases ±100/±101 versus the analytic limit/domain.
6. **Result:** **INCORRECT** — All are clipped to 99.999%, giving finite 7.231911... kcal/mol at 298.15 K without censoring metadata.
7. **Confidence:** High
8. **Limitations:** If clipping represents a detection limit, report the assumption and lower bound instead of an exact observation.

### C21 — Ni-hDA stored targets use their recorded 298.15 K temperature.

1. **StericX claim / audited proposition:** Ni-hDA stored targets use their recorded 298.15 K temperature. 
2. **Primary/reference source:** Corrected primary Ni-hDA SI Table S3; primary author CSV; NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** Experimental 80 °C, use 353.15 K (allow original rounded R).
4. **Implementation location:** `data/reactions_raw.csv; scripts/prepare_data.py::normalize_public_sigman; Study 011 temperature description`
5. **Independent validation method:** Independent reverse conversion of all 11 labeled source pairs.
6. **Result:** **INCORRECT** — Ten pairs imply 353.113701 K instead of 298.15 K; maximum room-temperature target difference ≈0.333685 kcal/mol.
7. **Confidence:** High
8. **Limitations:** Do not silently recompute targets; preserve published values and correct provenance separately.

### C22 — Ni-hDA ID2064 target and ee are mutually consistent.

1. **StericX claim / audited proposition:** Ni-hDA ID2064 target and ee are mutually consistent. 
2. **Primary/reference source:** Corrected primary Ni-hDA Table S3 and frozen author CSV
3. **Exact convention expected:** Use the same temperature and selectivity equation for each row.
4. **Implementation location:** `data/official/ni_hda_kraken.csv row 2064`
5. **Independent validation method:** Independently convert 3% ee at 353.15 K and infer T from the published target.
6. **Result:** **UNCERTAIN** — Source 3% ee and target ≈0.028072 correspond to ≈2% ee at 353 K; implied T=235.37 K.
7. **Confidence:** High confidence in the inconsistency; uncertain which datum is correct
8. **Limitations:** The disagreement exists upstream too; resolution requires author/experimental records.

### C23 — Population-weighted steric descriptors establish experimental conformer populations.

1. **StericX claim / audited proposition:** Population-weighted steric descriptors establish experimental conformer populations. Experimental interpretation tested for evidentiary support, not a claim that the README explicitly asserts measured conformer populations.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** True equilibrium populations depend on appropriate free energies, state counts, solvent and temperature.
4. **Implementation location:** `scripts/prepare_data.py; scripts/stericx_quantum.py; README ensembleclaims`
5. **Independent validation method:** Assess energy sources, degeneracy assumptions and failed-conformer handling against statistical thermodynamics.
6. **Result:** **OUT OF SCOPE** — Numerical weights are verified conditionally; no independent experimental population measurements were supplied.
7. **Confidence:** High about scope
8. **Limitations:** MMFF/xTB electronic-energy weighting and finite conformer search remain model assumptions.

## Models, statistics, uncertainty and screening

### M01 — Reaction feature vector is reproducibly constructed

1. **StericX claim / audited proposition:** Reaction feature vector is reproducibly constructed Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** [1,L,B1,B5,q,B1*q,B5*q,IR]; interactions in record precision.
4. **Implementation location:** `src/model/features.rs:15`
5. **Independent validation method:** Independent NumPy feature construction and frozen CLI fits.
6. **Result:** **VERIFIED** — Selected columns and native predictions match the declared algebra.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M02 — The fixed feature vocabulary makes the model mechanistically constrained

1. **StericX claim / audited proposition:** The fixed feature vocabulary makes the model mechanistically constrained Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Mechanistic constraints would require a justified relation to reaction mechanism; feature naming alone is insufficient.
4. **Implementation location:** `src/model/fit.rs:288`
5. **Independent validation method:** Inspect imposed restrictions and compare stated model name.
6. **Result:** **TERMINOLOGY ISSUE** — Actual constraints are a fixed vocabulary, term cap, pair-correlation cutoff and forward BIC; no enforced mechanism, sign or kinetic constraint.
7. **Confidence:** high
8. **Limitations:** Interpretable interactions can be useful hypotheses; chemical appropriateness is reaction-dependent.

### M03 — OLS coefficients minimize residual squares

1. **StericX claim / audited proposition:** OLS coefficients minimize residual squares Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm
3. **Exact convention expected:** Least-squares regression with intercept; unpenalized slopes.
4. **Implementation location:** `src/model/fit.rs:383`
5. **Independent validation method:** Independent NumPy SVD fits on 32 synthetic observations and 10 Ni-hDA observations.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum coefficient error 1.164e-7 (native) and 3.032e-8 (synthetic), consistent with stored f32 coefficients. The nominal OLS solve adds a 1e-10 ridge floor.
7. **Confidence:** high
8. **Limitations:** Strictly not exact unregularized OLS; floor negligible in tested well-conditioned designs, not necessarily singular designs.

### M04 — R² is calculated correctly

1. **StericX claim / audited proposition:** R² is calculated correctly Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/model_evaluation.html
3. **Exact convention expected:** 1−RSS/Σ(y−mean(y))²
4. **Implementation location:** `src/model/fit.rs:765; src/model/evaluation.rs`
5. **Independent validation method:** Independent NumPy and per-row residuals; results/*_math.json.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Native and synthetic values agree to coefficient rounding.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M05 — MAE is calculated correctly

1. **StericX claim / audited proposition:** MAE is calculated correctly Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/model_evaluation.html
3. **Exact convention expected:** mean(abs(y−prediction))
4. **Implementation location:** `src/model/fit.rs:765; src/model/evaluation.rs`
5. **Independent validation method:** Independent NumPy and per-row residuals; results/*_math.json.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Native and synthetic values agree to coefficient rounding.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M06 — RMSE is calculated correctly

1. **StericX claim / audited proposition:** RMSE is calculated correctly Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/model_evaluation.html
3. **Exact convention expected:** sqrt(mean((y−prediction)²))
4. **Implementation location:** `src/model/fit.rs:765; src/model/evaluation.rs`
5. **Independent validation method:** Independent NumPy and per-row residuals; results/*_math.json.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Native and synthetic values agree to coefficient rounding.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M07 — Fixed-feature LOO Q² measures held-out predictive residuals

1. **StericX claim / audited proposition:** Fixed-feature LOO Q² measures held-out predictive residuals Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/model_evaluation.html
3. **Exact convention expected:** Q²=1−PRESS/TSS; refit slopes/intercept after removing each observation.
4. **Implementation location:** `src/model/fit.rs:479`
5. **Independent validation method:** Independent SVD leave-one-out fits.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Independent Ni-hDA Q² = 0.00204170720964314; StericX = 0.002041716065248389.
7. **Confidence:** high
8. **Limitations:** Valid as a conditional fixed-feature diagnostic only, because selection was outside folds.

### M08 — Audit question: do fixed-feature diagnostics validate the complete selection workflow?

1. **StericX claim / audited proposition:** Audit question: do fixed-feature diagnostics validate the complete selection workflow? Audit question; not an attributed StericX quotation.
2. **Primary/reference source:** https://www.jmlr.org/papers/v11/cawley10a.html
3. **Exact convention expected:** Feature selection must be learned anew without each held-out response.
4. **Implementation location:** `src/model/fit.rs:203–243`
5. **Independent validation method:** Rerun the frozen CLI with each outer observation withheld; independently refit selection using SVD.
6. **Result:** **INCORRECT** — The actual CLI fails 7 of 10 outer folds because no descriptor improves the intercept model. No complete CLI outer-LOO Q² exists. An explicitly different reference workflow using training-mean fallback gives Q² = -1.5369635574.
7. **Confidence:** high
8. **Limitations:** Current fields explicitly say fixed_feature_loo. The negative result rejects a broader interpretation; the reference fallback is not StericX behavior.

### M09 — Train/test partition prevents held-out rows affecting scaling or fit

1. **StericX claim / audited proposition:** Train/test partition prevents held-out rows affecting scaling or fit Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.jmlr.org/papers/v11/cawley10a.html
3. **Exact convention expected:** Only train rows enter selection, means, variances, coefficients and resampling.
4. **Implementation location:** `src/model/training.rs:51; src/model/fit.rs:167`
5. **Independent validation method:** Perturb the held-out row descriptors, temperature and target to extreme values.
6. **Result:** **VERIFIED** — The full training report remains identical (holdout_partition.json).
7. **Confidence:** high
8. **Limitations:** This does not establish untouched experimental data or remove feature-selection leakage inside reported CV.

### M10 — Standardization is recomputed within regularization CV

1. **StericX claim / audited proposition:** Standardization is recomputed within regularization CV Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html
3. **Exact convention expected:** Use each inner-training fold mean and population SD.
4. **Implementation location:** `src/model/fit.rs:479,555,611`
5. **Independent validation method:** Independent nested LOO with sklearn regressors and freshly estimated scaling.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Nested-alpha ridge/LASSO errors match within f32 rounding.
7. **Confidence:** high
8. **Limitations:** Descriptor set still selected outside every outer fold.

### M11 — Ridge uses the stated regularized objective

1. **StericX claim / audited proposition:** Ridge uses the stated regularized objective Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html
3. **Exact convention expected:** RSS+αΣβ² in standardized coordinates; intercept unpenalized.
4. **Implementation location:** `src/model/fit.rs:383`
5. **Independent validation method:** sklearn Ridge SVD; same declared alpha grid.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Independent nested-alpha ridge RMSE = 0.732179928258; StericX = 0.732179919787.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M12 — LASSO uses the stated sparse objective

1. **StericX claim / audited proposition:** LASSO uses the stated sparse objective Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html
3. **Exact convention expected:** RSS/(2n)+αΣabs(β), standardized predictors, unpenalized intercept.
4. **Implementation location:** `src/model/fit.rs:423`
5. **Independent validation method:** sklearn Lasso with tight convergence, independently nested alpha.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Independent nested-alpha LASSO RMSE = 0.800454806613; StericX = 0.800454796181.
7. **Confidence:** high
8. **Limitations:** The production solver stops after 2000 cycles without a dual-gap report. Difficult correlated inputs beyond the tested designs remain unverified.

### M13 — The regularized nested_loo label fully describes what is nested

1. **StericX claim / audited proposition:** The regularized nested_loo label fully describes what is nested Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.jmlr.org/papers/v11/cawley10a.html
3. **Exact convention expected:** Both feature selection and alpha choice must be inside each outer fold.
4. **Implementation location:** `src/model/fit.rs:555`
5. **Independent validation method:** Trace selected argument; independently recompute full selection folds.
6. **Result:** **TERMINOLOGY ISSUE** — Alpha choice and scaling are nested. Feature selection is fixed using all outer training labels.
7. **Confidence:** high
8. **Limitations:** Interpret this as fixed-feature nested-alpha LOO, not full-pipeline predictive validation.

### M14 — BIC feature selection follows a declared rule

1. **StericX claim / audited proposition:** BIC feature selection follows a declared rule Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://doi.org/10.1214/aos/1176344136
3. **Exact convention expected:** BIC = n ln(RSS/n) + p ln(n); the additional forward-step rule requires a decrease greater than 2.
4. **Implementation location:** `src/model/fit.rs:337,794`
5. **Independent validation method:** Independent matrix-SVD forward search; compare chosen feature sets.
6. **Result:** **SUPPORTED** — Both synthetic and native selected sets agree.
7. **Confidence:** high
8. **Limitations:** The improvement cutoff, correlation cutoff 0.95, and n/3 term cap are design choices. Forward search need not find the global BIC minimum.

### M15 — VIF diagnoses selected-feature collinearity

1. **StericX claim / audited proposition:** VIF diagnoses selected-feature collinearity Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.statsmodels.org/stable/generated/statsmodels.stats.outliers_influence.variance_inflation_factor.html
3. **Exact convention expected:** VIF_j = 1/(1-R_j²), regressing predictor j on the other selected predictors.
4. **Implementation location:** `src/model/fit.rs:839`
5. **Independent validation method:** Inverse reference correlation matrix on two-feature synthetic data.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Independent inverse-correlation calculation gives 1.002144098 for the synthetic predictors; the native one-predictor case gives 1.
7. **Confidence:** high
8. **Limitations:** The denominator is capped below at 1e-9, so singular VIF is represented as 1e9 rather than infinity.

### M16 — Leave-group-out removes all rows with a group label

1. **StericX claim / audited proposition:** Leave-group-out removes all rows with a group label Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.jmlr.org/papers/v11/cawley10a.html
3. **Exact convention expected:** Whole specified group withheld; no per-row residue from that group in the fold.
4. **Implementation location:** `src/model/fit.rs:508`
5. **Independent validation method:** Independent group partitions and SVD refits; paired synthetic groups and native scaffold groups.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Independent group Q² = 0.001315520959; StericX = 0.001315540471.
7. **Confidence:** high
8. **Limitations:** Groups are supplied by the caller; missing labels use row identity. Feature selection still sees all outer labels.

### M17 — Audit question: is historical ligand 723 a scaffold-disjoint holdout?

1. **StericX claim / audited proposition:** Audit question: is historical ligand 723 a scaffold-disjoint holdout? Audit question; not an attributed StericX quotation.
2. **Primary/reference source:** https://github.com/SigmanGroup/Ni-Catalyzed-hDA/blob/main/Enantioselectivity_Model.ipynb
3. **Exact convention expected:** A scaffold holdout must share no chosen scaffold group with training.
4. **Implementation location:** `data/reactions_raw.csv; docs/study_001`
5. **Independent validation method:** RDKit canonical SMILES and Murcko scaffolds from frozen inputs.
6. **Result:** **INCORRECT** — Ligand 723 shares a Murcko scaffold with training ligands 1057 and 1058. No exact canonical isomeric SMILES duplicates were found among the eleven records.
7. **Confidence:** high
8. **Limitations:** Study 001 calls this a historical holdout; it does not claim scaffold-disjoint validation. Scaffold definition is the recorded RDKit Murcko convention.

### M18 — Bootstrap replicates measure fixed-feature coefficient uncertainty

1. **StericX claim / audited proposition:** Bootstrap replicates measure fixed-feature coefficient uncertainty Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://doi.org/10.1214/aos/1176344552
3. **Exact convention expected:** Resample paired training rows, refit chosen features, retain whole coefficient vectors.
4. **Implementation location:** `src/model/fit.rs:876`
5. **Independent validation method:** Inspect sampling unit; independently propagate stored vectors and rederive quantiles.
6. **Result:** **SUPPORTED** — Recomputed coefficient percentiles differ by at most 4.44e-16. Stored whole-vector propagation is consistent with the stated paired-row resampling design.
7. **Confidence:** high
8. **Limitations:** Replicate fits use ridge 1e-8, compared with a 1e-10 floor for nominal OLS. This checks the procedure and propagation, not small-sample coverage. Selection is fixed and ligand groups are not block-resampled.

### M19 — Audit question: does the fixed-feature permutation p-value test the full feature search?

1. **StericX claim / audited proposition:** Audit question: does the fixed-feature permutation p-value test the full feature search? Audit question; not an attributed StericX quotation.
2. **Primary/reference source:** https://www.jmlr.org/papers/v11/cawley10a.html
3. **Exact convention expected:** Repeat all response-driven selection on every shuffled target vector.
4. **Implementation location:** `src/model/fit.rs:946`
5. **Independent validation method:** 2000 paired independent target permutations, with fixed features versus repeated BIC selection.
6. **Result:** **INCORRECT** — Fixed-feature p = 0.06746626687; repeated-selection p = 0.25137431284.
7. **Confidence:** high
8. **Limitations:** Production describes fixed-feature permutations. Independent RNG draws differ from StericX; these are paired diagnostics of conditioning, not an exact RNG replay.

### M20 — Add-one response-permutation p-value is computed as stated

1. **StericX claim / audited proposition:** Add-one response-permutation p-value is computed as stated Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** p = (number at least as extreme + 1)/(B + 1), using training R² and fixed selected columns.
4. **Implementation location:** `src/model/fit.rs:946`
5. **Independent validation method:** Inspect counter and evaluate independent paired permutation distribution.
6. **Result:** **SUPPORTED** — Formula and conditioning match the declaration.
7. **Confidence:** high
8. **Limitations:** Counter arithmetic and design are supported by source inspection; the production RNG was not independently proven uniform.

### M21 — The 95% Student-t multiplier is numerically accurate on the audited degrees of freedom

1. **StericX claim / audited proposition:** The 95% Student-t multiplier is numerically accurate on the audited degrees of freedom Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html
3. **Exact convention expected:** Inverse two-sided tail at alpha = 0.05 and positive residual degrees of freedom.
4. **Implementation location:** `src/model/domain.rs:662`
5. **Independent validation method:** SciPy t.isf for df = 1, 2, 3, 8, 29, 100 and 1e6.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum absolute 95% quantile error = 1.720721343e-10.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M22 — Audit question: does the public Student-t inverse handle all accepted tail inputs?

1. **StericX claim / audited proposition:** Audit question: does the public Student-t inverse handle all accepted tail inputs? Audit question; not an attributed StericX quotation.
2. **Primary/reference source:** https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html
3. **Exact convention expected:** Return the finite quantile for positive df and 0 < alpha < 1, or explicitly report a numerical limit.
4. **Implementation location:** `src/model/domain.rs:662`
5. **Independent validation method:** SciPy t.isf and the analytic Cauchy tail, with alpha = 1e-8 and df = 1.
6. **Result:** **INCORRECT** — StericX returns 1048576; the reference is 63661977.23675813 (98.35% relative error).
7. **Confidence:** high
8. **Limitations:** The silently capped bracket does not affect the fixed 95% intervals at the audited positive integer degrees of freedom.

### M23 — OLS prediction interval includes future residual scatter

1. **StericX claim / audited proposition:** OLS prediction interval includes future residual scatter Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm
3. **Exact convention expected:** ŷ±t_(.975,n−p)s sqrt(1+h), s²=RSS/(n−p).
4. **Implementation location:** `src/model/domain.rs:368`
5. **Independent validation method:** Independent inverse design matrix and SciPy critical value.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Synthetic screening prediction-interval endpoint error <= 5.33e-15.
7. **Confidence:** high
8. **Limitations:** Conditional on a fixed correct linear mean and IID homoscedastic normal errors. Does not include descriptor, selection, mechanism, or laboratory-shift uncertainty.

### M24 — OLS mean-response confidence interval differs from prediction interval

1. **StericX claim / audited proposition:** OLS mean-response confidence interval differs from prediction interval Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm
3. **Exact convention expected:** ŷ±t s sqrt(h), omitting the new-observation residual term.
4. **Implementation location:** `src/model/domain.rs:490`
5. **Independent validation method:** Independent design inverse and SciPy critical values on training, interpolation and extrapolation points.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Confidence intervals use sqrt(h), whereas prediction intervals use sqrt(1+h); observed values match.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M25 — Screening distinguishes bootstrap mean-response uncertainty from a prediction interval

1. **StericX claim / audited proposition:** Screening distinguishes bootstrap mean-response uncertainty from a prediction interval Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://doi.org/10.1214/aos/1176344552
3. **Exact convention expected:** Joint coefficient-replicate propagation estimates fitted-mean uncertainty; a future observation also has residual scatter.
4. **Implementation location:** `src/commands/screen.rs:1349; README Prediction uncertainty`
5. **Independent validation method:** Recompute replicate predictions and their empirical 2.5% and 97.5% quantiles.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — The method is correctly named percentile_bootstrap_mean_response. Endpoint differences <= 1.11e-16 on the tested library.
7. **Confidence:** high
8. **Limitations:** The current README explicitly distinguishes these intervals. Neither bootstrap coverage nor uncertainty from selection or chemical shift is established.

### M26 — The marginal coefficient band is necessarily conservative

1. **StericX claim / audited proposition:** The marginal coefficient band is necessarily conservative Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://doi.org/10.1214/aos/1176344552
3. **Exact convention expected:** Marginal 95% coefficient intervals do not generally form a joint 95% confidence region.
4. **Implementation location:** `src/commands/screen.rs:1421; README:611`
5. **Independent validation method:** Create 1000 joint coefficient replicates with disjoint 2.4% positive and negative tails.
6. **Result:** **INCORRECT** — The marginal band collapses to one value and covers 90.4% of joint predictions; the joint 95% interval spans approximately 200 response units.
7. **Confidence:** high
8. **Limitations:** A mathematical serialized-model counterexample disproves a universal guarantee. It is not a measured chemistry error rate.

### M27 — Audit question: are nominal 95% intervals calibrated on the available Ni-hDA ranking panels?

1. **StericX claim / audited proposition:** Audit question: are nominal 95% intervals calibrated on the available Ni-hDA ranking panels? Audit question; not an attributed StericX quotation.
2. **Primary/reference source:** https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm
3. **Exact convention expected:** Approximately 95% coverage over the stated prediction population, subject to model assumptions.
4. **Implementation location:** `docs/study_011/rankings.csv`
5. **Independent validation method:** Recalculate coverage from every frozen Study 011 prediction and target.
6. **Result:** **INCORRECT** — 92/141 = 65.248% coverage; interpolation subset 70/108 = 64.815%.
7. **Confidence:** high
8. **Limitations:** Overlapping panels reuse eleven ligands. These are descriptive coverage fractions, not 141 independent trials. Aggregation mismatch complicates causal attribution.

### M28 — Applicability standardizes selected descriptors correctly

1. **StericX claim / audited proposition:** Applicability standardizes selected descriptors correctly Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html
3. **Exact convention expected:** z_j = (x_j - training mean_j)/training SD_j.
4. **Implementation location:** `src/model/domain.rs:301`
5. **Independent validation method:** Known anisotropic coordinates and independent affine transformation.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Training, nearby and extreme points match, subject to the f32 input representation.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M29 — Nearest-training distance is Euclidean in standardized space

1. **StericX claim / audited proposition:** Nearest-training distance is Euclidean in standardized space Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html
3. **Exact convention expected:** Minimum Euclidean distance to any standardized training point.
4. **Implementation location:** `src/model/domain.rs:389`
5. **Independent validation method:** SciPy cdist versus the public API observer and frozen screening CLI.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — CLI distance error is zero on the tested library; an exact training point has distance zero.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M30 — Training spacing thresholds are computed as documented

1. **StericX claim / audited proposition:** Training spacing thresholds are computed as documented Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Use each training point’s nearest other point, then max, mean + SD or mean + 2 SD as the chosen threshold.
4. **Implementation location:** `src/model/domain.rs:84,158`
5. **Independent validation method:** Inspect calibration formulas and recompute training nearest-neighbor distances independently.
6. **Result:** **SUPPORTED** — The default maximum threshold matches independent distances; alternate formulas match their declared algebra.
7. **Confidence:** high
8. **Limitations:** These are policies, not calibrated probabilities. Duplicates can make the threshold zero; an isolated training point can widen the maximum rule.

### M31 — Feature-range checks detect extremal extrapolation

1. **StericX claim / audited proposition:** Feature-range checks detect extremal extrapolation Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Flag values below a selected-feature minimum or above its maximum; equality is inside.
4. **Implementation location:** `src/model/domain.rs:688`
5. **Independent validation method:** Synthetic boundary, in-range and extreme inputs, including a point outside the convex hull but inside every range.
6. **Result:** **VERIFIED** — Range flags behave as declared. A bounding-box interior point outside the convex hull is still labeled interpolation.
7. **Confidence:** high
8. **Limitations:** A feature box and nearest-neighbor rule do not define chemical applicability or the training convex hull.

### M32 — Leverage matches the regression design

1. **StericX claim / audited proposition:** Leverage matches the regression design Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm
3. **Exact convention expected:** h = xᵀ(XᵀX)^(-1)x, with intercept and the recorded training standardization.
4. **Implementation location:** `src/model/domain.rs:324; src/model/fit.rs:661`
5. **Independent validation method:** NumPy inverse and SVD on known and fitted design matrices.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — CLI leverage error <= 7.91e-14 in the tested library.
7. **Confidence:** high
8. **Limitations:** The fit’s small ridge floor perturbs the unregularized equation. Leverage alone does not calibrate prediction error.

### M33 — Mahalanobis distance is recovered correctly in nonsingular designs

1. **StericX claim / audited proposition:** Mahalanobis distance is recovered correctly in nonsingular designs Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html
3. **Exact convention expected:** Mahalanobis distance = sqrt((n-1)(h-1/n)) for an ordinary full-rank linear design with intercept.
4. **Implementation location:** `src/model/domain.rs:442`
5. **Independent validation method:** Independent sample-covariance inverse on known two-dimensional geometry.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — CLI error <= 1.43e-13 on the nonsingular test library.
7. **Confidence:** high
8. **Limitations:** The identity requires full rank and the correct covariance convention.

### M34 — Audit question: does the public stored-geometry API identify singular covariance before reporting Mahalanobis distance?

1. **StericX claim / audited proposition:** Audit question: does the public stored-geometry API identify singular covariance before reporting Mahalanobis distance? Audit question; not an attributed StericX quotation.
2. **Primary/reference source:** https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html
3. **Exact convention expected:** Ordinary Mahalanobis distance needs an invertible covariance; a regularized surrogate should be identified.
4. **Implementation location:** `src/model/domain.rs:442`
5. **Independent validation method:** Inject a rank-one covariance stabilized with the same 1e-10 slope floor into TrainingGeometry.
6. **Result:** **INCORRECT** — Reports 316227.75293364684 and no unavailability reason for an off-manifold query.
7. **Confidence:** high
8. **Limitations:** This is an accepted serialized-geometry/API witness. Native forward selection did not generate this specific singular fitted model.

### M35 — The emitted trust label reliable denotes calibrated chemical reliability

1. **StericX claim / audited proposition:** The emitted trust label reliable denotes calibrated chemical reliability Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Reliable prediction would require held-out calibration beyond a descriptor range and nearest-neighbor rule.
4. **Implementation location:** `src/commands/screen.rs:1492; src/model/domain.rs:788`
5. **Independent validation method:** Trace the label decision and independently recalculate Study 011 errors by domain.
6. **Result:** **TERMINOLOGY ISSUE** — Interpolation MAE = 0.9332371 kcal/mol with 64.815% nominal 95% interval coverage.
7. **Confidence:** high
8. **Limitations:** README already states that descriptor distance is not calibrated reliability; the emitted label is stronger than the evidence.

### M36 — Ascending objective ranking

1. **StericX claim / audited proposition:** Ascending objective ranking Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Increasing predicted response
4. **Implementation location:** `src/commands/screen.rs:282,1805`
5. **Independent validation method:** Compare independent sorts with frozen synthetic-library outputs and metadata variants.
6. **Result:** **VERIFIED** — Observed ranking matches the declared objective.
7. **Confidence:** high
8. **Limitations:** Ranking arithmetic does not validate the chemical predictions.

### M37 — Descending objective ranking

1. **StericX claim / audited proposition:** Descending objective ranking Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Decreasing predicted response
4. **Implementation location:** `src/commands/screen.rs:282,1805`
5. **Independent validation method:** Compare independent sorts with frozen synthetic-library outputs and metadata variants.
6. **Result:** **VERIFIED** — Observed ranking matches the declared objective.
7. **Confidence:** high
8. **Limitations:** Ranking arithmetic does not validate the chemical predictions.

### M38 — Magnitude objective preserves negative predictions

1. **StericX claim / audited proposition:** Magnitude objective preserves negative predictions Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Decreasing absolute predicted response while retaining the original sign
4. **Implementation location:** `src/commands/screen.rs:282,1805`
5. **Independent validation method:** Compare independent sorts with frozen synthetic-library outputs and metadata variants.
6. **Result:** **VERIFIED** — Observed ranking matches the declared objective.
7. **Confidence:** high
8. **Limitations:** Ranking arithmetic does not validate the chemical predictions.

### M39 — Prediction ties use the declared deterministic order

1. **StericX claim / audited proposition:** Prediction ties use the declared deterministic order Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Resolve ties by identifier, then original position
4. **Implementation location:** `src/commands/screen.rs:282,1805`
5. **Independent validation method:** Compare independent sorts with frozen synthetic-library outputs and metadata variants.
6. **Result:** **VERIFIED** — Observed ranking matches the declared objective.
7. **Confidence:** high
8. **Limitations:** Ranking arithmetic does not validate the chemical predictions.

### M40 — Missing required descriptors are excluded explicitly

1. **StericX claim / audited proposition:** Missing required descriptors are excluded explicitly Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Exclude and explain candidates missing a descriptor used by the selected model.
4. **Implementation location:** `src/commands/screen.rs:1644`
5. **Independent validation method:** A frozen CSV includes one row without required Sterimol L.
6. **Result:** **VERIFIED** — The row is excluded with sterimol_l recorded as missing; the other nine candidates remain.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M41 — Previously tested exclusions are traceable

1. **StericX claim / audited proposition:** Previously tested exclusions are traceable Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Exclude matching stable identifiers and report identifiers that were not found.
4. **Implementation location:** `src/commands/screen.rs:978`
5. **Independent validation method:** An exclusion CSV contains high twice and not_present once.
6. **Result:** **VERIFIED** — The matching candidate is excluded; the unresolved identifier is reported.
7. **Confidence:** high
8. **Limitations:** Identifier matching is not chemical-identity matching.

### M42 — Diversity selection preserves prediction values

1. **StericX claim / audited proposition:** Diversity selection preserves prediction values Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Diversity changes selection or ordering, while preserving predictions; weight zero recovers ordinary ranking.
4. **Implementation location:** `src/commands/screen.rs:390,1839`
5. **Independent validation method:** Compare per-identifier predictions for default, diversity weight 0.65 and weight 0.
6. **Result:** **VERIFIED** — Maximum prediction change is zero; weight-zero ordering matches ordinary ranking.
7. **Confidence:** high
8. **Limitations:** Finite audited inputs; not a proof for every floating point input.

### M43 — OOD candidates are visible unless explicitly filtered

1. **StericX claim / audited proposition:** OOD candidates are visible unless explicitly filtered Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Include OOD candidates by default; an explicit domain-only filter records removals.
4. **Implementation location:** `src/commands/screen.rs:1760`
5. **Independent validation method:** Extreme descriptor candidates and the in-domain-only flag.
6. **Result:** **VERIFIED** — Default output retains extrapolations; the requested filter leaves a traceable exclusion record.
7. **Confidence:** high
8. **Limitations:** A warning does not establish model reliability.

### M44 — Published Ni-hDA univariate coefficients and metrics are reproduced

1. **StericX claim / audited proposition:** Published Ni-hDA univariate coefficients and metrics are reproduced Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Ni-Catalyzed-hDA/blob/main/Enantioselectivity_Model.ipynb
3. **Exact convention expected:** The official ten training IDs and vbur_max_delta_qvbur_min, with univariate OLS.
4. **Implementation location:** `studies/study_001_ni_hda.py`
5. **Independent validation method:** Fresh primary CSV and official notebook; independent SVD and fixed-feature LOO.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Slope = 0.0892042662475, intercept = -0.197934470338; R² = 0.8193060694 and conditional Q² = 0.7521012497.
7. **Confidence:** high
8. **Limitations:** The original search considered 191 descriptors on the same ten labels. Ligand 723 is historical, not prospective, and ID 2064 has conflicting primary ee/target values (M67).

### M45 — The native compact-feature ablation has weak conditional Ni-hDA validation

1. **StericX claim / audited proposition:** The native compact-feature ablation has weak conditional Ni-hDA validation Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Ni-Catalyzed-hDA/blob/main/Enantioselectivity_Model.ipynb
3. **Exact convention expected:** Reproduce the declared native fit without substituting the published buried-volume predictor.
4. **Implementation location:** `docs/study_001; src/model/fit.rs`
5. **Independent validation method:** Frozen CLI fits, independent SVD diagnostics and full outer-selection CLI attempts.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Native fixed-feature Q² = 0.002041716; seven of ten full-selection outer folds fail to fit.
7. **Confidence:** high
8. **Limitations:** Study 001 already reports this negative ablation. The mean-fallback reference Q² is not a completed StericX pipeline score.

### M46 — Study 002 Ni-hDA model arithmetic reproduces its frozen inputs

1. **StericX claim / audited proposition:** Study 002 Ni-hDA model arithmetic reproduces its frozen inputs Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Ni-Catalyzed-hDA/blob/main/Enantioselectivity_Model.ipynb
3. **Exact convention expected:** The same ten training IDs and the study’s declared per-ligand buried-volume predictor.
4. **Implementation location:** `studies/study_002`
5. **Independent validation method:** Independent SVD refits from frozen per-ligand inputs, with every fitted and LOO value saved.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Training, fixed-feature LOO and historical ligand-723 predictions reproduce within saved CSV precision.
7. **Confidence:** high
8. **Limitations:** This verifies model arithmetic, not the independently audited geometry-generation process or new experimental performance.

### M47 — Study 003 Ni-hDA model arithmetic reproduces its frozen inputs

1. **StericX claim / audited proposition:** Study 003 Ni-hDA model arithmetic reproduces its frozen inputs Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Ni-Catalyzed-hDA/blob/main/Enantioselectivity_Model.ipynb
3. **Exact convention expected:** The same ten training IDs and the study’s declared per-ligand buried-volume predictor.
4. **Implementation location:** `studies/study_003`
5. **Independent validation method:** Independent SVD refits from frozen per-ligand inputs, with every fitted and LOO value saved.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Training, fixed-feature LOO and historical ligand-723 predictions reproduce within saved CSV precision.
7. **Confidence:** high
8. **Limitations:** This verifies model arithmetic, not the independently audited geometry-generation process or new experimental performance.

### M48 — Audit question: has the prospective deck established experimental predictive validity?

1. **StericX claim / audited proposition:** Audit question: has the prospective deck established experimental predictive validity? Audit question; not an attributed StericX quotation.
2. **Primary/reference source:** https://github.com/SigmanGroup/Ni-Catalyzed-hDA/blob/main/Enantioselectivity_Model.ipynb
3. **Exact convention expected:** New measured outcomes obtained after predictions were fixed, without refitting.
4. **Implementation location:** `docs/study_003/PREREGISTRATION.md`
5. **Independent validation method:** Inspect the frozen protocol, prediction deck and available measurement status.
6. **Result:** **OUT OF SCOPE** — No prospective experimental outcome dataset was available for this audit.
7. **Confidence:** high
8. **Limitations:** The separate kinetics audit identifies a target-temperature mismatch; computational reproduction cannot supply missing experiments.

### M49 — Study 011 ranking performance is accurately reported

1. **StericX claim / audited proposition:** Study 011 ranking performance is accurately reported Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://scikit-learn.org/stable/modules/model_evaluation.html
3. **Exact convention expected:** Retain all 57 planned panels, including ten failures, and the locked failure penalties.
4. **Implementation location:** `docs/study_011; studies/study_011_ligand_ranking.py`
5. **Independent validation method:** Recompute 141 prediction rows, 47 successful panels, and all-57 summary denominators.
6. **Result:** **VERIFIED** — Top-1 = 0.1578947 versus random 1/3; top-2 set overlap = 0.4649123 versus random 2/3; mean successful-panel Spearman rho = -0.3404255.
7. **Confidence:** high
8. **Limitations:** Panels overlap; these are not 57 independent chemistry experiments.

### M50 — Screening uses the same descriptor aggregation as fitting

1. **StericX claim / audited proposition:** Screening uses the same descriptor aggregation as fitting Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** StericX-specific declared design; no external scientific validation is implied.
3. **Exact convention expected:** Candidate and training predictors must use the same descriptor and conformer-aggregation definitions.
4. **Implementation location:** `src/commands/screen.rs:1120; src/descriptors.rs; src/reaction.rs`
5. **Independent validation method:** Trace geometry input paths and independently recalculate frozen Study 011 prediction differences.
6. **Result:** **INCORRECT** — Largest screen-versus-fit disagreement = 1.20976108497 kcal/mol; representative geometry versus weighted ensemble is the documented mismatch.
7. **Confidence:** high
8. **Limitations:** The README and Study 011 disclose it. A correctly constructed precomputed descriptor CSV can avoid this specific mismatch.

### M51 — Cross-coupling reaction I reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction I reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_007_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **SUPPORTED** — N=34; current threshold 32.5285568237 (Left), MCC 0.624294473581; reference-descriptor threshold 32.4179439998 (Left), MCC 0.624294473581. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M52 — Cross-coupling reaction II reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction II reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_007_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **UNCERTAIN** — N=89; current threshold 32.9017391205 (Left), MCC 0.541895556035; reference-descriptor threshold 32.7365006523 (Left), MCC 0.541895556035. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M53 — Cross-coupling reaction III reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction III reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_007_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **SUPPORTED** — N=89; current threshold 31.756231308 (Left), MCC 0.470566123689; reference-descriptor threshold 31.5539973915 (Left), MCC 0.494451386058. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M54 — Cross-coupling reaction IV reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction IV reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_007_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **UNCERTAIN** — N=89; current threshold 31.9379549026 (Left), MCC 0.457947350849; reference-descriptor threshold 31.8923546413 (Left), MCC 0.457947350849. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M55 — Cross-coupling reaction V reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction V reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_007_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **SUPPORTED** — N=89; current threshold 50.8177585602 (Left), MCC 0.364079120849; reference-descriptor threshold 51.5269754011 (Left), MCC 0.364079120849. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M56 — Cross-coupling reaction RS1 reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction RS1 reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_007_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **UNCERTAIN** — N=89; current threshold 31.9379549026 (Left), MCC 0.546614927337; reference-descriptor threshold 31.8923546413 (Left), MCC 0.546614927337. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M57 — Cross-coupling reaction VII reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction VII reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_009_pd_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **UNCERTAIN** — N=55; current threshold 25.1200685501 (Right), MCC 0.368706994901; reference-descriptor threshold 25.1752935155 (Right), MCC 0.403437192056. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M58 — Cross-coupling reaction VIII reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction VIII reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_009_pd_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **UNCERTAIN** — N=53; current threshold 28.9362668991 (Right), MCC 0.73997002516; reference-descriptor threshold 28.8672311774 (Right), MCC 0.705275956512. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M59 — Cross-coupling reaction IX reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction IX reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_009_pd_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **UNCERTAIN** — N=30; current threshold 30.2667427063 (Right), MCC 0.782585580871; reference-descriptor threshold 30.5333931222 (Right), MCC 0.721687836487. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M60 — Cross-coupling reaction X reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction X reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_009_pd_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **SUPPORTED** — N=28; current threshold 30.5782699585 (Right), MCC 0.864063848822; reference-descriptor threshold 30.5333931222 (Right), MCC 0.864063848822. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M61 — Cross-coupling reaction XI reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction XI reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_009_pd_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **SUPPORTED** — N=30; current threshold 29.1050109863 (Right), MCC 0.782585580871; reference-descriptor threshold 28.8194767989 (Right), MCC 0.792527080644. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M62 — Cross-coupling reaction XII reproduces the published classifier

1. **StericX claim / audited proposition:** Cross-coupling reaction XII reproduces the published classifier Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Use the primary active convention yield >= cutoff, declared class weights and a depth-one weighted-Gini classifier.
4. **Implementation location:** `studies/study_009_pd_crosscoupling.py`
5. **Independent validation method:** Complete primary SI extraction; independent exhaustive cutpoint search versus sklearn, on reference and fresh StericX descriptors.
6. **Result:** **SUPPORTED** — N=71; current threshold 29.6112422943 (Right), MCC 0.376270264631; reference-descriptor threshold 29.5815085729 (Right), MCC 0.376270264631. All per-row reference/sklearn predictions agree.
7. **Confidence:** high for arithmetic; limited for exact published reproduction
8. **Limitations:** The complete open preprint SI is primary evidence; final Science SI retrieval returned 403. Printed yields and source versions do not exactly reproduce every published summary metric. All twelve model fits are resubstitution, not unseen-ligand validation.

### M63 — Cross-coupling active/inactive target convention matches reference

1. **StericX claim / audited proposition:** Cross-coupling active/inactive target convention matches reference Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://github.com/SigmanGroup/Threshold/blob/main/Sigman_single_parameter_threshold_analysis.ipynb
3. **Exact convention expected:** Official notebook: inactive if yield < cutoff, otherwise active (>=).
4. **Implementation location:** `studies/study_007_crosscoupling.py; studies/study_009_pd_crosscoupling.py`
5. **Independent validation method:** Extract and preserve every primary-table row exactly at a cutoff.
6. **Result:** **INCORRECT** — Study scripts use strict >. Seven row labels differ across I, II, IV, V, VII and VIII.
7. **Confidence:** high
8. **Limitations:** Published rounded yields can also differ from original unrounded measurements; that is separate from the definite code-inequality mismatch.

### M64 — The cross-coupling bootstrap intervals describe apparent-fit variability

1. **StericX claim / audited proposition:** The cross-coupling bootstrap intervals describe apparent-fit variability Paraphrase of documented behavior or implementation contract; see location.
2. **Primary/reference source:** https://doi.org/10.1214/aos/1176344552
3. **Exact convention expected:** Fitting and scoring the same resample evaluates apparent fit, not held-out prediction.
4. **Implementation location:** `studies/study_007_crosscoupling.py:bootstrap_ci; study_009 equivalent`
5. **Independent validation method:** Inspect the resampling unit and scoring set in both frozen study scripts.
6. **Result:** **SUPPORTED** — Each bootstrap sample is used for fitting and scoring. The reported intervals do not establish unseen-reaction or unseen-ligand calibration.
7. **Confidence:** high
8. **Limitations:** The interval may summarize resampled in-sample fit; no external predictive coverage follows from it.

### M65 — Study 007 withholds the reaction in its leave-one-reaction-out test

1. **StericX claim / audited proposition:** Study 007 withholds the reaction in its leave-one-reaction-out test Explicit study claim or audit question as stated.
2. **Primary/reference source:** https://www.jmlr.org/papers/v11/cawley10a.html
3. **Exact convention expected:** No target from the held reaction may enter its fit. Ligand identity overlap must be reported separately.
4. **Implementation location:** `studies/study_007_crosscoupling.py:333`
5. **Independent validation method:** Independent six-fold threshold fits and ID-set intersections.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Reaction targets are withheld. All 89 test ligands recur in training for five folds; 20/34 recur for reaction I.
7. **Confidence:** high within the stated scope
8. **Limitations:** This is reaction-level transfer among known ligand families, not scaffold-disjoint or generally unseen-ligand prediction.

### M66 — The Boltzmann descriptor is a scientifically justified reference for Study 007’s single-geometry comparison

1. **StericX claim / audited proposition:** The Boltzmann descriptor is a scientifically justified reference for Study 007’s single-geometry comparison Explicit study claim or audit question as stated.
2. **Primary/reference source:** Statistical-mechanical ensemble definition; see kinetics/conformer audit.
3. **Exact convention expected:** Justify the comparison target from the geometry and ensemble populations; mathematical equality is not required for a sensitivity comparison.
4. **Implementation location:** `docs/study_007/STUDY_007.md:59`
5. **Independent validation method:** Inspect the actual fair-reference wording and reported imperfect agreement; distinguish that methodology from a hypothetical exact-equivalence claim.
6. **Result:** **UNCERTAIN** — The study calls the Boltzmann descriptor the fair published reference and reports R² = 0.9735 and MAE = 0.446 percentage points. It does not assert mathematical equality. This audit has not independently established whether that reference choice is appropriate for all eighteen cases.
7. **Confidence:** high for source attribution; limited for methodological suitability
8. **Limitations:** No new eighteen-geometry comparison was run here. A separate audit hypothesis of universal single-geometry/ensemble equality is mathematically false, but is not a claim attributed to Study 007.

### M67 — The retained Ni-hDA target and ee columns are mutually consistent

1. **StericX claim / audited proposition:** The retained Ni-hDA target and ee columns are mutually consistent Explicit study claim or audit question as stated.
2. **Primary/reference source:** https://acs.figshare.com/articles/journal_contribution/30933433
3. **Exact convention expected:** At 353.15 K, absolute DDG = RT ln((100+ee)/(100-ee)).
4. **Implementation location:** `data/reactions_raw.csv; primary reaction_data.csv; corrected SI Table S3`
5. **Independent validation method:** Inspect rendered corrected SI page S22; independently convert every retained ee.
6. **Result:** **UNCERTAIN** — ID 2064 has ee=3% but target 0.028072105, corresponding to approximately 2%; independent 3% conversion gives 0.0421195099234 kcal/mol.
7. **Confidence:** high within the stated scope
8. **Limitations:** The primary sources conflict internally. Do not change StericX targets until authors or original measurements establish the intended value.

### M68 — The corrected Ni-hDA SI retains the model’s eleven ligand IDs and printed DDG targets

1. **StericX claim / audited proposition:** The corrected Ni-hDA SI retains the model’s eleven ligand IDs and printed DDG targets Explicit study claim or audit question as stated.
2. **Primary/reference source:** https://acs.figshare.com/articles/journal_contribution/30933433
3. **Exact convention expected:** Compare the actual corrected table, without assuming the correction changes or validates every value.
4. **Implementation location:** `models/sources/correction_si.pdf, Table S3, page S22`
5. **Independent validation method:** Row-by-row checked transcription against the primary CSV.
6. **Result:** **SUPPORTED** — The retained IDs and printed DDG column match; ee-to-DDG inconsistency remains for 2064.
7. **Confidence:** high within the stated scope
8. **Limitations:** The correction letter was inaccessible (403). Only the table and relevant modeling section were reviewed; no chromatogram reassessment or full correction scope is claimed.

### M69 — Frozen evaluation rejects predictions inconsistent with the supplied finite model

1. **StericX claim / audited proposition:** Frozen evaluation rejects predictions inconsistent with the supplied finite model Explicit study claim or audit question as stated.
2. **Primary/reference source:** Declared integrity contract in src/model/evaluation.rs module documentation.
3. **Exact convention expected:** Reject nonfinite predictions before comparing them to a finite model response or reporting accuracy.
4. **Implementation location:** `src/model/evaluation.rs:81; src/commands/evaluate.rs:51`
5. **Independent validation method:** Frozen CLI evaluation of three copied prediction CSVs: NaN, infinity and a finite mismatch.
6. **Result:** **INCORRECT** — NaN is accepted with exit 0; MAE and RMSE become NaN in stdout and null in JSON. Infinity and the finite mismatch are rejected with exit 2.
7. **Confidence:** high within the stated scope
8. **Limitations:** A concrete CLI validation failure. The single-row R² being unavailable is expected independently and is not itself the failure.

### M70 — Study 007’s formula guard rejects isomeric geometries before comparison

1. **StericX claim / audited proposition:** Study 007’s formula guard rejects isomeric geometries before comparison Explicit study claim or audit question as stated.
2. **Primary/reference source:** Distinct explicit molecular graphs and independent RDKit CalcMolFormula; source claim in Study 007 section 4.
3. **Exact convention expected:** Chemical identity requires connectivity and, when relevant, stereochemistry; equal element counts alone do not distinguish constitutional isomers.
4. **Implementation location:** `docs/study_007/STUDY_007.md:57; studies/study_007_crosscoupling.py:503–507`
5. **Independent validation method:** Freeze explicit CCCP and CCPC graphs, XYZ and SDF files; independently verify C3H9P and different P neighbor sets; execute the frozen formula helper definitions.
6. **Result:** **INCORRECT** — Primary n-propylphosphine and secondary methylethylphosphine both return ((C,3),(H,9),(P,1)). The cross-isomer pairing passes the formula guard, despite one versus two heavy neighbors at phosphorus.
7. **Confidence:** high within the stated scope
8. **Limitations:** Synthetic helper-level falsification of the universal guard claim only. No claim that any of the eighteen historical geometry mappings was actually wrong; other lookup and presence conditions were not rerun.

## Kraken reproduction

### K01 — Study 004/README: 1,541 validated ligands and 31,611 DFT conformers from 1,566 reference IDs.

1. **StericX claim / audited proposition:** Study 004/README: 1,541 validated ligands and 31,611 DFT conformers from 1,566 reference IDs. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Every available conformer is attempted; failed/native-incomplete ensembles are reported, never reduced from successful subsets.
4. **Implementation location:** `studies/study_004_scaled.py; data/ligand_db/kraken_phosphines.manifest.json`
5. **Independent validation method:** Fresh acquisition and independent SDF topology inventory, full native observation and complete-ensemble counting.
6. **Result:** **SUPPORTED** — Current public snapshot provides 1,546 ligands/31,721 conformers; 20 IDs lack DFT geometries. Full successful ensembles: BV 1,543, Sterimol 1,544, pyramidalization 1,546. Historical exact coverage count is a historical artifact, not the current independently eligible universe.
7. **Confidence:** High for frozen coverage.
8. **Limitations:** SDF connectivity comes from independent OpenBabel inference; available public snapshot can differ from original calculation archive.

### K02 — Minimum of max adjacent-quadrant differences agrees with Kraken at approximately R² = 0.9852.

1. **StericX claim / audited proposition:** Minimum of max adjacent-quadrant differences agrees with Kraken at approximately R² = 0.9852. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Published primary 3.5 Å sphere, 2.28 Å raw bond vector center, Bondi 1.17, H excluded, density 0.001, three orientations and complete-ensemble minimum.
4. **Implementation location:** `src/geometry/buried_volume.rs; studies/study_004_scaled.py`
5. **Independent validation method:** Fresh native outputs against published values; all conventional error statistics, residuals, full 20 worst ensembles under four untuned Morfeus variants.
6. **Result:** **SUPPORTED** — N 1,543, MAE 0.270886692 Å³, RMSE 0.490777865 Å³, max 4.559474133 Å³, R² = 0.985128630. Original-convention Morfeus reduces all 20 headline residuals to ≤ 0.006275429 Å³. Approximate agreement is reproduced; exact convention equivalence is not.
7. **Confidence:** High for measured results.
8. **Limitations:** Default density and center differ; aggregate accuracy does not establish physical or predictive validity.

### K03 — Study 004 describes the entire Kraken vbur family as reproduced.

1. **StericX claim / audited proposition:** Study 004 describes the entire Kraken vbur family as reproduced. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Every claimed field/reduction must be checked separately; original source includes distal/total volumes, ratios and Boltzmann reductions beyond eight geometric fields.
4. **Implementation location:** `studies/study_004_vbur_family.py::FAMILY; primary source boltzproperties/mmproperties`
5. **Independent validation method:** Enumerate primary properties and compare every implemented matching unweighted reduction: 56 metrics across 14 fields × 4 reductions.
6. **Result:** **TERMINOLOGY ISSUE** — Eight absolute BV fields are a subset; even 56 current comparisons do not include source-energy Boltzmann or unimplemented distal/total families.
7. **Confidence:** High.
8. **Limitations:** Unimplemented source properties are OUT OF SCOPE for numerical SUT reproduction; missing weights remain UNCERTAIN rather than verified.

### K04 — Study 004 reports six Sterimol extrema under matched Kraken coordination-axis convention.

1. **StericX claim / audited proposition:** Study 004 reports six Sterimol extrema under matched Kraken coordination-axis convention. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Original virtual Pd→P axis with raw bond vector center, H radius 1.09 Å, 3,600 directions, +0.40 Å L.
4. **Implementation location:** `src/geometry/sterimol.rs; src/descriptors.rs; studies/study_004_sterimol.py`
5. **Independent validation method:** All 31,721 conformers in independent primary-convention Morfeus; explicit same center/radii and radius-only variants; full minima/maxima/other reductions.
6. **Result:** **INCORRECT** — Strict matched-convention wording is false: native H radius 1.20 Å, unit bond vector center and 360 directions differ. Six native extrema MAEs 0.0963–0.1124 Å; max L error 0.98638 Å against published values.
7. **Confidence:** High for conventions and measurements.
8. **Limitations:** A different documented descriptor convention can be legitimate; choosing one requires scientific intent, not tuning against these outputs.

### K05 — Native P/alpha reproduce the intended Morfeus construction for ordinary three-neighbor Kraken geometries.

1. **StericX claim / audited proposition:** Native P/alpha reproduce the intended Morfeus construction for ordinary three-neighbor Kraken geometries. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Same donor, coordinates and three selected neighbors; P dimensionless and alpha degrees.
4. **Implementation location:** `src/geometry/pyramidalization.rs`
5. **Independent validation method:** Independent contemporary Morfeus on all 31,721 exported conformers; nearest-three and SDF-topology variants independently agree in this corpus.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Max identical-geometry |ΔP| = 3.2010947048632943e−7; |Δalpha| = 4.8856123654239525e−5 degrees.
7. **Confidence:** High within tested corpus.
8. **Limitations:** Does not establish undocumented degeneracies or every automatic neighbor policy; see geometry audit. Published-value extrema are a separate claim K09.

### K06 — STUDY 004_RESIDUAL attributes the Kraken DFT center to an xTB localized orbital and says it coincides with the geometric center for tertiary phosphines.

1. **StericX claim / audited proposition:** STUDY 004_RESIDUAL attributes the Kraken DFT center to an xTB localized orbital and says it coincides with the geometric center for tertiary phosphines. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Actual published DFT code sums raw P − neighbor vectors then normalizes; native sums individually normalized vectors. Neither is an electronic orbital center in this DFT path.
4. **Implementation location:** `docs/study_004/STUDY_004_RESIDUAL.md; src/geometry/buried_volume.rs::lone_pair_direction`
5. **Independent validation method:** Immutable source inspection plus original center/native center explicit comparisons across 590 diagnostic conformers, including all 20 headline ensembles.
6. **Result:** **INCORRECT** — Largest headline case ID 1796 has zero P–H bonds, center angle≈ 9.715 degrees and 4.55947 Å³ discrepancy; original center/density resolves it to 0.00418 Å³. Many major cases have no P–H bonds.
7. **Confidence:** High.
8. **Limitations:** Finding is specific to this DFT path; other xTB/LMO workflows are not inferred identical.

### K07 — Matched Kraken convention implies identical total, near/far and finite-grid outputs.

1. **StericX claim / audited proposition:** Matched Kraken convention implies identical total, near/far and finite-grid outputs. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Original source retains last of three frame evaluations for total/near/far; native retains first. Quadrant/octant extrema combine all three.
4. **Implementation location:** `src/geometry/buried_volume.rs::compute_from_center; primary morfeus_properties`
5. **Independent validation method:** Direct reference variants explicitly match first orientation before native comparison, then change only density or original center.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Expanded matched-default comparison covers 2,372 cases: total BV max error 0.01166136723496436 Å³; percent BV max error 0.006494365624803322 percentage points; adjacent-quadrant max error 0.011656810711858867 Å³. These are approximately one grid cell. The earlier 581-case sample had a smaller total error; expansion is retained.
7. **Confidence:** High for measured finite-grid effects.
8. **Limitations:** Matched frame/density is necessary; first-versus-last convention is scientifically equivalent only in the continuum, not bitwise on finite grids.

### K08 — Committed table contains native descriptors and declares bond-axis Sterimol, rather than containing the published coordination-axis values.

1. **StericX claim / audited proposition:** Committed table contains native descriptors and declares bond-axis Sterimol, rather than containing the published coordination-axis values. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Artifact provenance/declared axis must be distinguished from scientific numerical accuracy.
4. **Implementation location:** `data/ligand_db/kraken_phosphines.manifest.json; src/commands/db.rs`
5. **Independent validation method:** Inspect frozen manifest, source table and consumers separately from all fresh scientific comparisons.
6. **Result:** **VERIFIED** — Frozen manifest explicitly declares computed native values and sterimol_axis=bond. This establishes identity/declared convention only.
7. **Confidence:** High for artifact metadata.
8. **Limitations:** Does not independently verify every historical packed record or authorize exchanging bond-axis and coordination-axis model descriptors.

### K09 — Study 005 interprets four extrema residuals as essentially export-coordinate/Rust-precision error.

1. **StericX claim / audited proposition:** Study 005 interprets four extrema residuals as essentially export-coordinate/Rust-precision error. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Any attribution must explain molecule-specific residual magnitude with independent matched-geometry calculations.
4. **Implementation location:** `docs/study_005; studies/study_005_pyramidalization.py`
5. **Independent validation method:** Full corpus Morfeus and native; 1,028,096 rounding corner calculations on 251 major outlier conformers.
6. **Result:** **UNCERTAIN** — Native/published max alpha discrepancy 1.155986692 degrees is reproduced independently (1.156001260 degrees). Same geometry native error≤ 4.89e−5 degrees; sampled rounding deviation≤ 0.010351413 degrees overall and 0.008927286 degrees for worst ID 1290. Simple precision attribution is unsupported. A conservative continuous rounding-box bound gives maximum |ΔP| of 0.000623445508668 across all 251 examined conformers, including both ID 1290 cases; its published P discrepancy is about 0.007218735, assuming unchanged ensemble and atom identity.
7. **Confidence:** High that tested numerical effects are insufficient; provenance cause uncertain.
8. **Limitations:** Corner sensitivity is not a rigorous interior bound. Original full precision logs/ensembles and historical software needed; do not simply blame published dataset.

### K10 — Comments describe B5 as exact maximum radial support and imply an exact alignment for all valid finite axes.

1. **StericX claim / audited proposition:** Comments describe B5 as exact maximum radial support and imply an exact alignment for all valid finite axes. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** L=max axial support+0.40; B5=max perpendicular distance+radius on intended axis.
4. **Implementation location:** `src/geometry/sterimol.rs::compute_with_dummy; glam0.30.10::Quat::from_rotation_arc`
5. **Independent validation method:** Independent dot/perpendicular formulas on 31,713 successful native axes, corroborated with direct Morfeus on both worst real conformers.
6. **Result:** **INCORRECT** — L error 0.0036931285331647246 Å at KRAKEN:16:31923; B5 error 0.004176904379527002 Å at KRAKEN:1407:54385. Near parallel/antiparallel approximation discards a small real tilt.
7. **Confidence:** High; exact original failing inputs and independent results retained.
8. **Limitations:** Most errors≈ 5e−7 Å; magnitude is small relative to many ligand differences but nonzero and larger than claimed exactness. No production fix made.

### K11 — Geometric proximity reliably supplies three bonded substituents for intended library phosphines.

1. **StericX claim / audited proposition:** Geometric proximity reliably supplies three bonded substituents for intended library phosphines. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Known molecular topology should not gain nonbonded nearby P/Si atoms solely because of generous radius cutoff.
4. **Implementation location:** `src/geometry/buried_volume.rs::donor_neighbor_indices; studies/study_004_scaled.py prefilter`
5. **Independent validation method:** Fresh independent SDF graph and original Kraken connectivity compared with native; distances and element thresholds preserved for failures.
6. **Result:** **INCORRECT** — Eight conformers from IDs 1281/1907 get a fourth neighbor (P or Si); SDF and original primary connectivity have three. Both BV and coordination Sterimol error.
7. **Confidence:** High for explicit disagreement with two independent connectivity conventions.
8. **Limitations:** Chemical bond assignment still inferred; no experimental bond-order truth claimed. The historical prefilter mirrored native cutoff and hid this failure mode.

### K12 — Native BV rejects nonzero occupancy with global max_delta_qvbur=0 as degenerate.

1. **StericX claim / audited proposition:** Native BV rejects nonzero occupancy with global max_delta_qvbur=0 as degenerate. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** A symmetrically occupied valid sphere can have zero anisotropy; finite grid zero does not prove an invalid frame.
4. **Implementation location:** `src/geometry/buried_volume.rs::compute_from_center`
5. **Independent validation method:** Frozen actual PH3-like ligand 1299 conformer 54318; independent geometry/frame and Morfeus occupancy comparison.
6. **Result:** **INCORRECT** — Native raises symmetric-zero error for valid three-neighbor input; independent calculation succeeds.
7. **Confidence:** High with actual input and reference output.
8. **Limitations:** Do not infer every zero is valid; validate geometric degeneracy independently. Synthetic symmetric cases examined in geometry audit.

### K13 — Full reproduction could extend published ligand-level Boltzmann values from available geometries.

1. **StericX claim / audited proposition:** Full reproduction could extend published ligand-level Boltzmann values from available geometries. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Need per-conformer source corrected free energies/populations, temperature/degeneracy, selection and identity mapping. Published aggregates are comparison outputs only.
4. **Implementation location:** `studies/study_004_*; original read_ligand weighting`
5. **Independent validation method:** Current/v2 API schemas/endpoints on five conformers, alternate exports, original repo/release tree, full SI archive, updated Sigman repo checked and frozen.
6. **Result:** **UNCERTAIN** — Historical weights not recovered. API returns data:null; SI contains aggregate tables; updated repo 33 new ensembles aredifferentcalculations. Full corpus Boltzmann not reproduced.
7. **Confidence:** High in checked-source limitation, not universal absence.
8. **Limitations:** Obtain original source energies; synthetic Boltzmann verification elsewhere does not validate historical populations.

### K14 — Matching primaryconventions suffices to exactly reproduce all published Sterimol values fromcurrentAPI geometries.

1. **StericX claim / audited proposition:** Matching primaryconventions suffices to exactly reproduce all published Sterimol values fromcurrentAPI geometries. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Same original geometry/ensemble, axis, radii and historical dependency implementation are required.
4. **Implementation location:** `Primary API and original Kraken/SI artifacts; a provenance claim, not exclusively a native-code claim`
5. **Independent validation method:** Contemporary Morfeus primary convention on all 31,721 geometries, complete ensembles, all outliers retained.
6. **Result:** **UNCERTAIN** — Remaining largest errors: L max 0.955032324 Å (ID 1290), B5 min 0.649228687 Å (ID 1805), B1 max 0.274270968 Å (ID 560). Bothnativeandindependentreferencedisagreewithsomepublished extrema.
7. **Confidence:** High in measured disagreement; attribution uncertain.
8. **Limitations:** Need original calculation provenance before declaring either the public table or current implementation wrong.

### K15 — Study 005 says the native kernel agrees with Morfeus to 4.4e−16 for P and 2.8e−14 degrees for alpha.

1. **StericX claim / audited proposition:** Study 005 says the native kernel agrees with Morfeus to 4.4e−16 for P and 2.8e−14 degrees for alpha. 
2. **Primary/reference source:** Gensch et al., JACS (2022), [DOI 10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718); primary DFT code at immutable commit `4eaad505c1343e6083032b4a3fda47e004e19734`; frozen MolSSI API/SI inputs. See [Kraken report](kraken/REPORT.md) for sources and hashes.
3. **Exact convention expected:** Those numbers must describe native Rust f32 output if attributed to the native kernel, not a Python float64 identity check.
4. **Implementation location:** `docs/study_005/STUDY_005.md:21; src/geometry/pyramidalization.rs`
5. **Independent validation method:** Actual native values and IEEE bits frozen for all 31,721 conformers, compared independently with Morfeus on the identical exported geometry.
6. **Result:** **INCORRECT** — Measured maxima are 3.2010947048632943e−7 for P and 4.8856123654239525e−5 degrees for alpha. A Python algebraic identity at float64 precision cannot support the stated native-output precision.
7. **Confidence:** High, with full corpus native and independent records.
8. **Limitations:** The measured native errors are small in ordinary geometries; this finding corrects numerical-precision attribution, not every aspect of the underlying definition.
