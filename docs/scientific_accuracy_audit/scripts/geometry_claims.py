# ruff: noqa: E501, RUF001
"""Render separately classified geometry claims; no combined PASS labels."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S = "https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html"
B = "https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html"
P = "https://doi.org/10.1007/BF00676621"
C = "https://doi.org/10.1039/B801115J"
V = "https://doi.org/10.1021/j100785a001"
M = "https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html"
claims = []


def add(
    id,
    claim,
    source,
    convention,
    location,
    method,
    result,
    status,
    confidence,
    limitations,
):
    claims.append(
        dict(
            id=id,
            claim=claim,
            primary_or_reference_source=source,
            exact_convention_expected=convention,
            implementation_location=location,
            independent_validation_method=method,
            result=result,
            status=status,
            confidence=confidence,
            limitations=limitations,
        )
    )


add(
    "G01",
    "Raw Sterimol L is the maximum axial atomic-sphere extent.",
    S,
    "max_i[(r_i-o)·u + rho_i], excluding a real dummy. Å.",
    "src/geometry/sterimol.rs::compute/params_from_projection",
    "Analytic support equation in independent float64 + Morfeus, 31 finite base cases.",
    "Max error 6.994274599492201e-7 Å for bond-axis base cases.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High for controlled valid inputs",
    "Restricted base cases are not a global accuracy bound. Real Kraken near-axis cases expose a separate 0.00172232 Å L alignment error (G39). Does not validate automatic axis chemistry, unknown radii, large-coordinate overflow, or historical original radius parameterization.",
)
add(
    "G02",
    "B1 is the continuous minimum transverse support radius.",
    S,
    "min over all transverse unit directions of max_i(p_i·v+rho_i).",
    "src/geometry/sterimol.rs::params_from_projection",
    "Analytic finite candidate enumeration of support-envelope stationary points/intersections; minimized four-atom witness.",
    "SUT 1.7261793613433838 Å vs exact 1.7000000000000002 Å; a 0.0261793613433836 Å overestimate from one-degree sampling.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "Approximation is documented in kernel comments. Exact B1 equality/strict invariance is false; error scales with transverse size. A retained full Kraken case has error 0.0416578313 Å; the minimized witness is not the audit-wide maximum. See G39 for the distinct alignment contribution.",
)
add(
    "G03",
    "B5 is the maximum transverse atomic-sphere extent.",
    S,
    "max_i(||p_i||+rho_i), Å.",
    "src/geometry/sterimol.rs::params_from_projection",
    "Independent analytic radial support; compare Morfeus sampled maximum separately.",
    "Max bond-axis base error 5.23377122085833e-7 Å; Morfeus B5 itself uses an approximate angular scan.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "The base maximum is not a global bound: a real near-antiparallel Kraken case has B5 error 0.00276961144 Å (G39). Extreme coordinates also require separate limits.",
)
add(
    "G04",
    "Coordination-axis L applies the historical +0.40 Å convention.",
    S + " ; https://aarontools.readthedocs.io/en/latest/api/substituent.html",
    "Raw API L plus 0.40 Å in CLI coordination mode; bond mode raw.",
    "src/descriptors.rs::sterimol_for_conformer; src/cli.rs::STERIMOL_L_CORRECTION",
    "Read frozen source and compare raw API / Morfeus L_value_uncorrected and L_value.",
    "Arithmetic/convention placement agrees; raw API and CLI use distinct conventions.",
    "SUPPORTED",
    "High implementation; moderate historical attribution",
    "Original Verloop chapter not fully accessible. The correction is a historical parameter convention, not a universal physical requirement for arbitrary metal dummies.",
)
add(
    "G05",
    "Explicit hydrogen spheres contribute to Sterimol and real dummy is excluded.",
    S,
    "All non-dummy atoms with their chosen radii; virtual dummy includes real donor.",
    "src/geometry/sterimol.rs::compute/compute_with_dummy",
    "Independent explicit-H graph molecules and envelope calculation.",
    "Ordinary geometries agree within reported B1 scan/f32 limits.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "Kraken DFT code uses H radius1.09 Å; StericX uses1.20 Å, so published reproduction is separate.",
)
add(
    "G06",
    "Default bond-axis descriptors are independent of atom ordering.",
    "Geometric invariance requirement; " + S,
    "Same chemical structure and chosen chemical axis should produce the same descriptor.",
    "src/descriptors.rs::detect_donor; src/geometry/buried_volume.rs::bonded_neighbors",
    "Frozen actual CLI on all six permutations of three equally distant chemically distinct substituents.",
    "L changes6.0699997→3.5 Å and B5 changes3.5→6.0699997 Å because the tie selects the first atom index.",
    "INCORRECT",
    "High",
    "The underlying explicitly selected-axis kernel is invariant to reordering within f32 limits. This falsifies an unrestricted auto-axis invariance claim, not the mathematics for a specified axis.",
)
add(
    "G07",
    "Sterimol is rigid-rotation/translation invariant.",
    S,
    "Physical continuous descriptors invariant for fixed chemical axis.",
    "src/geometry/sterimol.rs::compute/compute_with_dummy",
    "6 representatives×100 seeded random rotations plus X/Y/Z rotations and random translations.",
    "B1 max deviation0.019944190979003906 Å; L≤1.0013580322265625e-5 Å, B5≤4.76837158203125e-6 Å.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High for measured range",
    "Rotations include modest translations ≤50 Å; no universal bound. B1 scan phase depends on alignment. A targeted exact 90-degree rotation of real near-axis Kraken cases changes L by 0.00172233582 Å and B5 by 0.00276851654 Å (G39), despite exact independent invariance.",
)
add(
    "G08",
    "Buried volume integrates the union of atomic spheres inside a sphere.",
    M + " ; " + B,
    "V=integral of union membership inside R; %Vbur=100V/(4πR³/3).",
    "src/geometry/buried_volume.rs::occupied_volumes",
    "Independent float64 union-of-balls integration and Morfeus under matched center/radii/grid convention.",
    "23 finite base values max total-volume difference3.6215270000639066e-6 Å³; default lattice includes15,408 points.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High for compatible default cases",
    "Numerical quadrature error is larger than implementation residual; invalid symmetry rejection and unknown radii are distinct failures.",
)
add(
    "G09",
    "The default sphere radius3.5 Å follows a published descriptor convention.",
    M,
    "Sphere centered on selected real/virtual coordination point, R3.5 Å.",
    "src/geometry/buried_volume.rs::BuriedVolumeConfig; src/cli.rs",
    "Official SambVca manual and Morfeus documentation.",
    "Convention agrees.",
    "VERIFIED",
    "High",
    "This empirically useful convention does not prove a universal physically optimal sphere radius.",
)
add(
    "G10",
    "Bondi radius scaling1.17 follows the reference convention.",
    M + " ; " + B,
    "Every included vdW radius multiplied by1.17 before intersection.",
    "src/geometry/buried_volume.rs::aligned_atoms",
    "Official manual plus independent scaled-radius calculations.",
    "Scaling agrees on supported atoms and finite ordinary values.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "Input radii table and overflow behavior must be validated separately.",
)
add(
    "G11",
    "Hydrogens are excluded from buried-volume occupancy by default.",
    M + " ; " + B,
    "H may define connectivity/frame but contributes no sphere unless include_hydrogens=true.",
    "src/geometry/buried_volume.rs::aligned_atoms",
    "Explicit-H methylphosphine/trimethylamine protocol pairs and PH3 analytic lens witness.",
    "Exclusion/inclusion matches controlled Morfeus reference.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "Removing donor-bound H from the frame would be a different, incorrect chemical operation.",
)
add(
    "G12",
    "All supported radii are exactly Morfeus Bondi radii.",
    V + " ; frozen/reference_sources/morfeus/data.py",
    "Declared table equality, not only the label Bondi-style.",
    "src/geometry/xyz.rs::van_der_waals_radius",
    "Read independent table; per-element synthetic witnesses.",
    "Br is1.85 Å in SUT vs1.83 in Morfeus; B1.92 in SUT vs Morfeus fallback2.0.",
    "INCORRECT",
    "High table difference",
    "The original Bondi table was not available for direct inspection, so the Br discrepancy was not adjudicated against that table; B has no Morfeus Bondi entry; Mantina et al., DOI 10.1021/jp8111556, explicitly gives the later 1.92 Å extension. Difference alone does not identify physical truth.",
)
add(
    "G13",
    "Unknown-element radius fallbacks are chemically conservative.",
    V,
    "An arbitrary fallback must not be equated with an element-specific physical radius.",
    "src/geometry/xyz.rs::van_der_waals_radius/covalent_radius",
    "He/Xe/Na/Se/Fe toy witnesses; official/reference tables.",
    "SUT1.8 Å can be both too large(He) and too small(Na/Xe); max observed BV difference6.072736800747158 Å³,3.3813625981379403 percentage points.",
    "UNCERTAIN",
    "High numerical result; low chemical justification",
    "Morfeus also falls back2.0 for absent B/Fe entries; it is not ground truth. Broad chemical safety of either fallback lacks evidence.",
)
add(
    "G14",
    "The filled integration grid reproduces Morfeus for every density.",
    B,
    "Float64 linspace and norm≤R in Morfeus versus f32 SUT arithmetic; grid endpoints included.",
    "src/geometry/buried_volume.rs::integration_grid",
    "Four density levels plus odd-n grid and representable-float sphere-boundary perturbations.",
    "Default n32 counts/membership agree for tested cases; odd grids reveal region-boundary differences, fine grids can differ at sphere lattice boundaries.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "Not exact all-density compatibility. See convergence and partition records; sphere-boundary choices have measure zero only in continuum.",
)
add(
    "G15",
    "Voxel atomic occupancy includes the vdW boundary.",
    B + " ; frozen Morfeus cKDTree query_ball_point implementation",
    "distance²≤scaled_radius², corresponding to a closed ball.",
    "src/geometry/buried_volume.rs::occupied_volumes",
    "Frozen one-point tests at1f32 and nextafter inside/outside, independent exact-real comparison.",
    "Inside/on occupied; next representable outside unoccupied.",
    "VERIFIED",
    "High for tested boundary",
    "Open/closed spheres have equal continuum volume; this establishes a lattice convention, not a uniquely physical inequality.",
)
add(
    "G16",
    "Sphere boundary handling includes points on R.",
    B,
    "Retain ||point||≤R.",
    "src/geometry/buried_volume.rs::integration_grid",
    "Private-grid observation at n5 with R and adjacent representable radii.",
    "Boundary points are present; outside lattice points rejected.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "f32 rounding can move mathematically marginal points. Atom centers outside R must still contribute when their sphere intersects R, confirmed by boundary-atom witnesses.",
)
add(
    "G17",
    "Quadrant volumes reproduce the reference regional convention.",
    B,
    "Canonical xy signs (++),(-+),(--),(+-); Morfeus excludes zero planes.",
    "src/geometry/buried_volume.rs::quadrant_index/occupied_volumes",
    "Every individual quadrant in every orientation, independent/Morfeus comparison.",
    "Default even-grid base values agree within f32 rounding. n5 qvbur_min differs9.621127681084193 Å³; odd-grid assignment differs structurally.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "SUT assigns zeros positive and normalizes each unequal-size region; no unrestricted Morfeus equivalence claim is justified.",
)
add(
    "G18",
    "Octant volumes reproduce the reference regional convention.",
    B,
    "Four xy quadrants for positive and negative z; explicit Morfeus ID remapping.",
    "src/geometry/buried_volume.rs::octant_index/occupied_volumes",
    "Every individual octant per orientation compared, canonical signs preserved.",
    "Even-grid agreement is strong; zero-plane regions at odd n differ.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "Numeric octant IDs differ from Morfeus for negative z, so direct numeric-ID comparisons would be wrong; mapping is explicit in the harness.",
)
add(
    "G19",
    "Near and far describe donor-facing and distal hemispheres.",
    B + " ; Kraken official SI page12",
    "Center→donor is negative z; near z<0, far z≥0 in SUT.",
    "src/geometry/buried_volume.rs::coordinate_basis/occupied_volumes",
    "Independent coordinate frame plus per-octant sums.",
    "Hemisphere naming/orientation agrees under the stated coordinate convention.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "This is not the outside-sphere distal molecular volume. Odd-grid normalization creates partition inconsistency.",
)
add(
    "G20",
    "Near+far and regional partitions consistently represent the same total volume.",
    "Continuum additivity of disjoint integrals; " + B,
    "Regional volume sum should converge to total, with quantified finite-grid discrepancy.",
    "src/geometry/buried_volume.rs::occupied_volumes",
    "Independent conservation check across densities.",
    "Triphenylphosphine coarse total53.1568947 versus near+far59.7332687 Å³ (+12.37%); default even-grid error≤rounding.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "An exact-additivity claim for all accepted configurations is incorrect. At very fine odd grid near+far remains 0.5790825 Å³ above total; refinement parity matters. Morfeus also has a coarse-grid partition gap of +3.4082601 Å³, so this is not evidence that its convention is physically true.",
)
add(
    "G21",
    "Maximum adjacent-quadrant difference is reduced correctly.",
    B + " ; official Kraken DFT code",
    "max over three plane orientations and four cyclic adjacent pairs of absolute difference.",
    "src/geometry/buried_volume.rs::compute_from_center",
    "Independent explicit twelve-quadrant reduction.",
    "Default matched cases agree within f32 subtraction rounding.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High",
    "Frame/center selection, zero-plane treatment and symmetry rejection are separate issues.",
)
add(
    "G22",
    "A positive buried volume with zero quadrant asymmetry is unphysical/degenerate.",
    "Sphere intersection symmetry; independent analytic lens; " + B,
    "An axisymmetric nonzero occupied set has equal quadrants and valid zero asymmetry.",
    "src/geometry/buried_volume.rs::compute_from_center final rejection",
    "PH3 with three chemically bonded H; full-cover toy; small sphere wholly inside ordinary donor radius.",
    "PH3 correctly has31.937214517315244 Å³ lattice volume,17.782969885773625%,zero asymmetry; SUT errors despite valid frame. Full-cover100% also errors.",
    "INCORRECT",
    "High, constructive counterexamples",
    "PH3 CLI separately rejects lack of heavy reference atom; public BV API accepts H as reference then hits invalid symmetry rejection. Small-sphere case also reproduces through CLI with config.",
)
add(
    "G23",
    "Buried-volume default precision supports finer-than-discretization scientific distinctions.",
    "Analytic two-sphere intersection formula",
    "Compare lattice to continuum, not only another lattice implementation.",
    "src/geometry/buried_volume.rs::integration_grid/occupied_volumes",
    "Exact PH3 one-atom lens and four densities.",
    "Analytic31.669141713568127 Å³; default error+0.2680721841369511 Å³ or0.149265 percentage points; density.001 error−0.0284346976; .0001 error+0.0006218515 Å³.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High for lens witness",
    "Four real/toy molecules default-versus-very-fine total differences up to0.1971054 Å³; no universal numerical error guarantee follows from five examples.",
)
add(
    "G24",
    "Buried-volume values are rigid-transformation invariant.",
    B,
    "Co-rotating molecular frame and center; finite-grid sensitivity quantified.",
    "src/geometry/buried_volume.rs::coordinate_basis/occupied_volumes",
    "618 transformed cases across six representatives, all public fields and raw region values.",
    "All nine public BV fields had exactly0 measured deviation on this campaign.",
    "SUPPORTED",
    "High for campaign; limited generality",
    "Grid transitions can occur for other geometries; huge translations lose f32 geometric information and are not covered by this favorable result.",
)
add(
    "G25",
    "Three covalently bonded substituents include donor-bound hydrogen.",
    C + " ; explicit molecular graphs",
    "SMILES graph bonds define chemically intended neighbors; H included.",
    "src/geometry/buried_volume.rs::bonded_neighbors/donor_neighbor_indices",
    "12 independent RDKit graph molecules plus primary/secondary phosphine contact witnesses.",
    "All12 standard graphs match; both H-bearing historical witnesses recover chemically bonded neighbors.",
    "SUPPORTED",
    "High for test molecules",
    "Geometry-only inference is not universal chemistry; charge/bond order/coordination identity are ignored.",
)
add(
    "G26",
    "Cordero covalent radii justify a universal1.3× bond cutoff.",
    C,
    "Cordero provides empirical atomic radii, not universal binary bonding truth at1.3× sums.",
    "src/geometry/buried_volume.rs::BOND_TOLERANCE_FACTOR/bonded_neighbors",
    "Primary source interpretation, threshold nextafter witnesses, independent graphs.",
    "Threshold behavior is deterministic and matches stated heuristic; scientific universality unsupported. Coincident atoms count as bonded before later geometry validation.",
    "UNCERTAIN",
    "High implementation; insufficient chemical universality",
    "No threshold retuning was performed. Stretching a synthetic bond beyond cutoff is not itself proof the chemistry is misclassified.",
)
add(
    "G27",
    "Automatic donor selection identifies chemically appropriate donors generally.",
    C + " ; CLI contract",
    "Sole chosen element, or explicit index, plus exactly3 inferred neighbors and≥1heavy neighbor.",
    "src/descriptors.rs::detect_donor",
    "Frozen CLI runs with phosphines/amines, multiple P, planar/coincident cases and PH3.",
    "Standard tested P/N donors work; multiple-element ambiguity is rejected; PH3 rejected for lack of heavy reference; chemistry beyond geometry unestablished.",
    "SUPPORTED",
    "Moderate",
    "Explicit index can select any element; trivalent geometry does not establish lone-pair availability, charge state, aromaticity or catalytic coordination.",
)
add(
    "G28",
    "Current frame fixes the historical nearest-heavy error.",
    "Independent known SMILES connectivity; primary/secondary phosphine valence",
    "Three bonded neighbors, including H, rather than three nearest heavy atoms.",
    "src/geometry/buried_volume.rs::donor_neighbor_indices",
    "Reconstruct former selector independently; controls plus H and nearby nonbonded heavy atoms.",
    "Tertiary control center unchanged; primary/secondary wrong-frame centers displaced2.4302383/2.0563763 Å and BV descriptors differ strongly; current neighborhoods match graphs.",
    "SUPPORTED",
    "High for reconstructed failure mechanism",
    "The exact six historically zero ligands (575,1485,1487,1490,1491,1495),20 conformers, were then independently reconstructed from primary SDF bond graphs; every old nearest-heavy frame again gives zero while current matched-frame values agree with independent references. Kraken reproduction conventions remain separate.",
)
add(
    "G29",
    "The geometric virtual-metal placement exactly matches published Kraken.",
    "Official Kraken SI FigureS2/page10 and frozen DFT source audited in kraken/",
    "Kraken DFT sums raw substituent displacements then normalizes; SUT sums unit bond directions.",
    "src/geometry/buried_volume.rs::infer_lone_pair_direction",
    "Primary code/SI comparison, independent center construction.",
    "Conventions differ for unequal bond lengths; matching center-distance2.28 alone is insufficient.",
    "INCORRECT",
    "High",
    "Separate Kraken audit quantifies published-data consequences. Geometric placement is itself an approximation to actual coordination geometry.",
)
add(
    "G30",
    "Planar fallback is a physically determined lone-pair direction.",
    "Geometric symmetry; no universal primary rule identified",
    "Planar center has no unique opposite-bond-sum direction; plane normal sign is ambiguous without environment/electronics.",
    "src/geometry/buried_volume.rs::infer_lone_pair_direction/center_clearance",
    "Planar and nearly planar synthetic structures, rotations and permutations.",
    "Fallback normal/maximum-clearance criterion is a StericX design choice; independent naive normalized-sum centers become unstable near exact planarity.",
    "UNCERTAIN",
    "High description; insufficient physical evidence",
    "Controlled-center kernel comparisons deliberately do not treat SUT-selected center as scientifically established.",
)
add(
    "G31",
    "Pyramidalization P implements the Radhakrishnan/Morfeus expression.",
    P
    + " ; https://digital-chemistry-laboratory.github.io/morfeus/api/morfeus.pyramidalization.html",
    "|det(unit donor vectors)| with2−P acute correction.",
    "src/geometry/pyramidalization.rs::compute",
    "Independent determinant calculation, Morfeus, tetrahedral/orthogonal/acute witnesses.",
    "Ordinary graph cases maxP error1.2938381854787906e-7; tetrahedral0.769800305 versus analytic0.7698003589.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High for runnable reference; moderate original-paper definition",
    "Original paper full text not acquired. Degenerate sentinel behavior is not a valid P measurement.",
)
add(
    "G32",
    "Pyramidalization alpha uses the mean signed plane-normal angle in degrees.",
    P + " ; Morfeus0.8 source",
    "Three alpha values, acute bisector sign, degree conversion; ordinary planar90°, orthogonal0°.",
    "src/geometry/pyramidalization.rs::compute",
    "Independent plane geometry and explicit acute/planar/tetrahedral cases.",
    "Ordinary graph maxerror6.122375019401716e-6°; acute example−82.717° is permitted.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High in nondegenerate range",
    "Near-collinear geometry is assigned0 by SUT despite finite independent29.99992497364023°; sentinel must not be treated scientifically.",
)
add(
    "G33",
    "P=1 corresponds to an ideal tetrahedral apex.",
    "Analytic tetrahedral unit vectors; " + P,
    "Three tetrahedral directions determinant4/(3√3)≈0.76980036; orthogonal triad1.",
    "src/geometry/pyramidalization.rs::PyramidalizationParams documentation",
    "Analytic determinant and frozen directAPI.",
    "The comment saying1 for ideal tetrahedral-like apex is misleading; tests compute tetrahedral0.7698.",
    "TERMINOLOGY ISSUE",
    "High",
    "Tetrahedral-like is imprecise; no universal physical normalization to1 for tetrahedral geometry is established.",
)
add(
    "G34",
    "Pyramidalization is neighbor-order invariant.",
    P,
    "Geometric three-vector measure independent of permutation.",
    "src/geometry/pyramidalization.rs::compute",
    "120 molecular permutations with mapped neighbor identities.",
    "No scientifically meaningful order variation observed; full residuals in invariance.csv.",
    "VERIFIED WITH NUMERICAL LIMIT",
    "High for tested domain",
    "Planar/degenerate numerical branches need explicit handling.",
)
add(
    "G35",
    "Invalid or nearly degenerate geometry produces an explicit scientific failure.",
    "Mathematical nondefinition of directions/plane normals",
    "Do not silently interpret sentinel zero as a measured descriptor.",
    "src/geometry/sterimol.rs::compute; src/geometry/pyramidalization.rs::compute",
    "Coincident, short-axis, nearly collinear, NaN/Inf input witnesses.",
    "Sterimol and pyramidalization direct APIs returnzero sentinels; BV errors. Near-collinear alpha29.999925° becomes0°; short nonzeroaxis givesallzeroSterimol.",
    "INCORRECT",
    "High for unrestricted explicit-error claim",
    "Zero sentinel behavior is documented at kernel level. CLI combined-descriptor validation often rejects these cases; callers of public kernels remain exposed.",
)
add(
    "G36",
    "Finite input coordinates always produce finite scientific outputs.",
    "IEEE754 finite arithmetic range; analytic norm",
    "Overflow must be detected if finite output is promised.",
    "src/geometry/sterimol.rs::params_from_projection",
    "Frozen finite1e38 Å remote atom witness.",
    "B5 becomes+Inf for wholly finite positions/radii; observer recordsnull plus explicitnonfinite map.",
    "INCORRECT",
    "High",
    "Deliberately extreme and not ordinary chemistry. It identifies an API validation limit, not evidence that f64 is needed for typical coordinates.",
)
add(
    "G37",
    "f32 is sufficiently accurate for ordinary geometry in this audit.",
    "Independent f64 equations and Morfeus under controlled protocol",
    "Separate rounding error from grid/angle discretization and coordinate convention.",
    "src/geometry/*.rs",
    "Graph molecules, rigid campaign, f32-rounded input re-evaluated in f64.",
    "Restricted ordinary-case L/B5 errors are mostly 1e-6–1e-5 Å and default BV residuals a few 1e-6 Å³. The targeted real near-Z cases instead have ~0.002 Å alignment errors. Scientific sufficiency requires an application-specific tolerance or effect-size assessment.",
    "UNCERTAIN",
    "High measured errors; insufficient evidence for universal scientific sufficiency",
    "Real near-axis Kraken geometries have larger alignment errors (G39), beyond ordinary rounding. Large translations 1e7–1e20 destroy coordinate detail; nearly degenerate axes amplify rounding. No blanket recommendation to switch f64 is made.",
)
add(
    "G38",
    "Geometry descriptor agreement validates predictive/experimental chemistry.",
    "Scientific identifiability; held-out experimental validation required",
    "Numerical descriptor compatibility and chemistry prediction are different hypotheses.",
    "README scientific summaries; studies002/004/005",
    "Separate scope of independent computation from model/generalization evidence.",
    "This geometry audit establishes no experimental predictive validity.",
    "OUT OF SCOPE",
    "High",
    "Requires held-out chemical/experimental evidence evaluated in model audit and prospective experiments.",
)
add(
    "G39",
    "B5 is the exact radial support about the supplied attachment axis, apart from ordinary f32 arithmetic rounding.",
    S + " ; https://docs.rs/glam/0.30.10/src/glam/f32/sse2/quat.rs.html",
    "Projection onto the supplied unit axis must remain invariant under a rigid change of global Cartesian coordinates; a near-Z axis must not be replaced by Z.",
    "src/geometry/sterimol.rs::compute/compute_with_dummy, Quat::from_rotation_arc in frozen glam 0.30.10",
    "Independent float64 dot/norm projections on identical SUT-rounded coordinates/radii/center; two real Kraken conformers minimized to donor plus active atom; exact 90-degree X coordinate permutation.",
    "Full-case L is 6.859999656677246 versus independent 6.858277336508549 Å (963/48394). Full-case B5 is 7.116985321044922 versus independent 7.114215709601416 Å (1075/49864). Exact coordinate rotation reduces corresponding errors to 1.56e-8 and 1.09e-6 Å; original minimized witnesses retain the failure.",
    "INCORRECT",
    "High, preserved full and minimized cases plus targeted invariance experiment",
    "The supplied axes are mathematically well defined and the source geometries are real DFT structures. glam documents approximate near-singular alignment (|dot| > 1−2 f32 EPSILON); this is a StericX accuracy-contract issue, not an undocumented glam defect. It is distinct from B1 angular discretization and does not prove a chemistry prediction is materially affected. See geometry/alignment/ and its rotation/ subdirectory.",
)
# Separate documented formulas/design choices from stronger audit propositions.
# Numerical agreement on a bounded population must not classify an unrestricted
# exact-equivalence or scientific-resolution proposition as verified.
scope_updates = {
    "G02": {
        "claim": "The documented one-degree B1 estimate approximates the continuous minimum support radius within the measured errors.",
        "claim_origin": "Scoped numerical audit of the one-degree scan documented in the kernel; not a claim that the source promises exact continuous minimization.",
    },
    "G14": {
        "claim": "Audit proposition: the filled integration grid exactly reproduces Morfeus for every accepted density.",
        "claim_origin": "Exact-compatibility audit question, not a quoted universal promise in the source. The documented SUT grid is evaluated against the independent reference grid.",
        "status": "INCORRECT",
        "limitations": "This rejects unrestricted exact all-density equivalence, not the shared continuum target or measured default-grid agreement. Floating-point sphere-boundary membership and region conventions differ; see convergence and partition records.",
    },
    "G17": {
        "claim": "Audit proposition: quadrants use exactly the same finite-grid regional convention as Morfeus at every accepted density.",
        "claim_origin": "Reference-convention compatibility question, not an assertion that an open or closed continuum boundary is physically preferable.",
        "status": "INCORRECT",
        "limitations": "Default even-grid agreement remains verified within the reported numerical limits. At odd grids SUT assigns zeros positive and directly normalizes quadrants; Morfeus excludes zero planes and sums octants. Both approximate the same continuum regions without exact numerical convention equivalence.",
    },
    "G18": {
        "claim": "Audit proposition: octants use exactly the same finite-grid regional convention as Morfeus at every accepted density.",
        "claim_origin": "Reference-convention compatibility question; it does not identify either implementation as physical truth.",
        "status": "INCORRECT",
        "limitations": "Default even-grid numerical agreement remains favorable. Odd-grid zero-plane assignment differs. Negative-z numeric IDs also differ, but the harness correctly remaps them; this classification concerns the actual boundary convention, not a mistaken ID comparison.",
    },
    "G20": {
        "claim": "Audit proposition: each finite-grid regional partition adds exactly to the independently estimated total for all valid configurations.",
        "claim_origin": "Numerical conservation audit question, distinct from exact additivity of the continuum integral; not a quoted all-grid guarantee.",
        "status": "INCORRECT",
        "limitations": "Default even-grid additivity holds within rounding in the campaign. The odd-grid discrepancy is a finite-estimator limitation, not a failure of continuum volume additivity. At density 0.0001 near+far remains 0.5790825 Å³ above total. Morfeus also has a coarse-grid gap of +3.4082601 Å³.",
    },
    "G23": {
        "claim": "Audit proposition: default buried-volume precision supports scientific distinctions finer than its discretization error.",
        "claim_origin": "Scientific-resolution audit question; no explicit application-specific resolution guarantee was identified in the cited kernel. The measured numerical floor does not itself establish a useful chemical effect size.",
        "status": "UNCERTAIN",
        "confidence": "High for the measured lens/convergence errors; insufficient evidence for scientific resolution below them",
        "limitations": "The analytic lens has default error 0.2680722 Å³, and four real/toy molecules differ from very fine totals by up to 0.1971054 Å³. Error cancellation for a particular descriptor difference, experimental relevance, and a universal resolution guarantee are not established. Agreement with another finite grid cannot supply that missing evidence.",
    },
    "G24": {
        "claim": "Buried-volume values are rigid-transformation invariant over the tested moderate-coordinate campaign.",
        "claim_origin": "Scoped numerical invariance proposition: six prespecified representatives and 618 transforms, not a universal finite-precision guarantee.",
    },
    "G27": {
        "claim": "Automatic donor selection identifies the expected donors in the tested standard P/N graphs.",
        "claim_origin": "Scoped test of automatic-selection behavior against independently supplied molecular graphs; general chemical donor availability remains unestablished.",
    },
    "G31": {
        "claim": "Pyramidalization P implements the determinant and acute-correction convention used by Morfeus 0.8.0.",
        "claim_origin": "Runnable-reference implementation proposition. Morfeus attributes the measure to Radhakrishnan and Agranat; original-paper textual equivalence was not established because full-text access was unavailable.",
    },
}
for claim in claims:
    claim.update(scope_updates.get(claim["id"], {}))
(ROOT / "geometry/claims.json").write_text(json.dumps(claims, indent=2) + "\n")
lines = [
    "# Geometry claim inventory",
    "",
    "Each result is independently classified. An implementation-level agreement is not evidence of chemical predictive validity. Detailed evidence is in [geometry/REPORT.md](geometry/REPORT.md), with equations and primary-source limitations in [geometry/METHODS.md](geometry/METHODS.md).",
    "",
]
for c in claims:
    lines += [
        f"## {c['id']}: {c['status']}",
        "",
        f"1. **StericX claim / audited proposition:** {c['claim']} "
        + c.get("claim_origin", ""),
        f"2. **Primary/reference source:** {c['primary_or_reference_source']}",
        f"3. **Exact expected convention:** {c['exact_convention_expected']}",
        f"4. **Implementation:** `{c['implementation_location']}`",
        f"5. **Independent validation:** {c['independent_validation_method']}",
        f"6. **Result:** {c['result']}",
        f"7. **Confidence:** {c['confidence']}",
        f"8. **Limitations:** {c['limitations']}",
        "",
    ]
(ROOT / "claims_geometry.md").write_text("\n".join(lines) + "\n")
print(len(claims), "claims")
