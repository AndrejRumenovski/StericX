# Independent geometry audit methods

Initial source references were consulted on 19 September 2026; the final geometry review added a source-access check and dependency inspection on 21 September 2026. Frozen HTTP bodies and access failures are retained in `sources/manifest.json`. Original papers that cannot be accessed are identified as access limitations; an inaccessible DOI is not evidence that the equation was verified from the paper. [Source-access details](sources/SOURCE_ACCESS.md) distinguish full text, abstracts, official software and unsuccessful downloads.

## System and independence

`inputs.json` was generated and hashed before any reference comparison. SMILES connectivity in the RDKit-generated structures supplies an independent graph: the expected donor neighbors do not come from StericX's distance cutoff. Synthetic structures carry explicit intentions; no chemical connectivity claim is made for artificially stretched threshold structures. ETKDG/MMFF is a geometry generator, not an experimental structural reference. The study includes unusual structures precisely to look for failures, not to imply they are stable molecules.

`frozen_buried_volume.rs` contains a byte-identical copy of the frozen source followed by `observer_append.rs`. The added functions expose previously private grids, aligned coordinates and individual region volumes. These observations belong to the **system under investigation**. They are not an independent reference calculation, and they do not establish scientific correctness. The ordinary frozen executable remains the reference for actual CLI behavior. No production source was changed.

`geometry_reference.py` imports NumPy, SciPy and Morfeus, never StericX. Its independent Sterimol B1 algorithm analytically enumerates stationary directions and pairwise intersections of atom support curves, instead of repeating StericX's one-degree scan. Its independent volume calculation uses float64 grid construction and cKDTree union-of-balls membership, with explicit geometric region masks. The Morfeus results are a third implementation; agreement with Morfeus establishes implementation compatibility under stated conventions, not physical truth.

## Sterimol equations and conventions

For unit axis **u** from dummy/origin **o** toward attached atom, each non-dummy atom has axial coordinate `z_i = (r_i-o)·u`, transverse vector `p_i = (r_i-o)-z_i u`, and van der Waals radius `rho_i`. The supported continuous-sphere geometric definitions are

- `L_raw = max_i(z_i + rho_i)`; Morfeus reports `L = L_raw + 0.40 Å`.
- `B5 = max_i(||p_i|| + rho_i)`.
- `B1 = min_{||v||=1,v⊥u} max_i(p_i·v + rho_i)`.

These equations follow the support function of the projected union of atomic spheres. For a finite envelope of sinusoids, a minimum occurs at an individual stationary point or a crossing of two support functions. All such candidates are evaluated, with no angular tuning against the SUT.

Original historical reference: [Verloop, Hoogenstraaten & Tipker, *Drug Design* 7 (1976), 165–207](https://doi.org/10.1016/B978-0-12-060307-7.50010-9). The original chapter could not be read through the DOI endpoint during this audit; exact historical wording and original radius parameterization therefore remain limited. The executable convention is independently available in [official Morfeus Sterimol documentation](https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html), frozen version 0.8.0 source, and [AaronTools' discussion of alternative L conventions](https://aarontools.readthedocs.io/en/latest/api/substituent.html). A +0.40 correction is a convention with historical cap-atom motivation, not a universal physical constant. Raw API L and corrected CLI L must not be conflated.

StericX uses Bondi-style radii, includes explicit hydrogens in Sterimol, and excludes the real attachment atom in bond-axis mode. In virtual-dummy mode every real atom, including the donor, contributes. Whole free-ligand projection along one donor-substituent bond is a StericX-specific axis choice; it is not automatically equivalent to a substituent isolated by cutting a molecular graph. Morfeus defaults to CRC radii, so this audit explicitly supplies Bondi radii in controlled comparisons. Published Kraken hydrogen-radius conventions are examined separately in the Kraken audit.

## Buried-volume equations and conventions

For center **c**, sphere radius `R`, selected atoms **a_i** and scaled radii `s rho_i`, the continuum target is

`V_occ = ∫_{||x-c||≤R} 1[∃i: ||x-a_i||≤s rho_i] dx`,

`%Vbur = 100 V_occ / (4πR³/3)`.

Regional volumes use the same integral restricted to a stated sign region. Boundaries are measure zero in the continuum but have nonzero weight on a finite grid. Thus equality inequalities, zero-plane assignment and per-region normalization must be specified for numerical compatibility.

[Poater et al.'s official SambVca manual](https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html) recommends a 3.5 Å sphere, Bondi radii multiplied by 1.17 and exclusion of hydrogens. These are calibrated descriptor conventions, not uniquely correct physical atomic sizes. The method's original publication is [Poater et al., *Eur. J. Inorg. Chem.* (2009), 1759–1766](https://doi.org/10.1002/ejic.200801160); the later implementation is [Falivene et al., *Organometallics* (2016), 35, 2286–2293](https://doi.org/10.1021/acs.organomet.6b00371). [Morfeus official documentation](https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html) independently documents those defaults. Morfeus and the published Kraken DFT protocol use density 0.001 Å³; the frozen StericX CLI and this audit’s primary comparison use its unchanged default of 0.01 Å³. A matched-default implementation comparison does not reproduce the published integration density.

The Morfeus filled projection grid has `n = round((8R³/density)^(1/3))` equally spaced coordinates including each cube face. Points satisfying `||x||≤R` are retained. The resulting actual sphere point count is not generally `(4πR³/3)/density`. Total volume is the occupied point fraction multiplied by analytic sphere volume. Regional volume is the regional fraction multiplied by analytic sphere volume divided by four/eight. At odd n, strict sign masks omit zero planes in Morfeus; StericX assigns zero planes to positive regions. This is tested rather than presumed negligible.

Canonical region IDs here are quadrants `(++),(-+),(--),(+-)` and octants in that order first for z≥0, then z<0. Morfeus's negative-z numeric IDs differ; the explicit mapping is `[0,1,2,3,7,6,5,4]`. With the center→donor direction aligned to negative z, negative z is the donor-facing **near** hemisphere. This is a geometrical convention, not the separate distal ligand-volume descriptor outside the integration sphere.

For three chemically selected plane atoms the audit compares every region of every orientation and the extrema of all 12 quadrant / 24 octant values. The maximum adjacent-quadrant difference is `max_{orientation,j}|q_j-q_(j-1 mod4)|`. The center and plane are controlled explicitly wherever possible. A geometrically guessed lone-pair direction is an approximation and must be evaluated separately from the integration algorithm.

## Pyramidalization equations

For normalized donor→neighbor vectors `a,b,c`, the ordinary branch is `P=|det[a,b,c]|`. Each alpha is the angle between one vector and the appropriately directed normal of the other two, signed negative when that vector projects positively onto their bisector. `alpha` is their mean in degrees; the acute branch uses `P=2-P` when the mean is negative. Independent calculation uses determinant and three geometric plane angles. The source is [Radhakrishnan & Agranat, “Measures of pyramidalization”, *Struct. Chem.* (1991), 2, 107–115](https://doi.org/10.1007/BF00676621), with exact runnable implementation documented by [Morfeus](https://digital-chemistry-laboratory.github.io/morfeus/api/morfeus.pyramidalization.html). The original article's inaccessible full text limits historical-definition verification. The descriptor is not simply restricted to [0,1] under the acute correction. Ideal tetrahedral three-neighbor directions give `P=4/(3√3)=0.7698003589`, `alpha=35.26438968°`; an orthogonal triad gives P=1. A planar ordinary center has P=0 and alpha=90°. Degenerate planes do not define alpha; reporting a numerical zero there is a sentinel convention, not a valid pyramidalization measurement.

## Radii, connectivity and chemical assumptions

[B​​ondi, *J. Phys. Chem.* (1964), 68, 441–451](https://doi.org/10.1021/j100785a001) provides the historical van der Waals reference; [Cordero et al., *Dalton Trans.* (2008), 2832–2838](https://doi.org/10.1039/B801115J) provides crystallographically estimated covalent radii. A 1.3 multiplier on the sum of covalent radii is a heuristic connectivity criterion chosen by StericX; the Cordero paper does not establish that every pair within that cutoff is bonded. Radii fallback values for unknown elements are implementation choices, not independently validated chemistry. Chemically bound hydrogens belong in donor geometry even when hydrogen spheres are excluded from volume integration.

## Error summaries and numerical campaigns

Descriptor-wise errors use SUT minus independent reference, with MAE, RMSE, maximum and median absolute error, ordinary least-squares slope/intercept (SUT versus reference), and identity-line R² `1-SSE/Σ(reference-mean)²`. Constant-reference cases have undefined slope/R². No outlier is removed. R² never substitutes for an absolute error assessment. Rotations use 100 seeded Haar-distributed random SO(3) rotations plus separate X/Y/Z rotations for each of six molecular representatives, each with a random translation. The campaign reports maximum absolute deviation, mean signed deviation, mean absolute deviation and standard deviation from its original structure. Twenty full atom permutations per representative preserve explicit chemical identity mappings.

Convergence uses unchanged defaults for the primary comparison, then densities 0.1, 0.01, 0.001 and 0.0001 Å³, plus 0.008 to force an odd grid. Finer-grid values are convergence diagnostics, not physical ground truth. Precision differences are investigated separately from changes in conventions or input geometries. Large translations are deliberately ill-conditioned f32 inputs; their failures do not imply ordinary molecular coordinates have those errors.

## Targeted near-axis alignment experiment

Two full Kraken geometries and donor-plus-active-atom reductions retain identical SUT-rounded coordinates, radii and explicit center. Independent direct dot-product and radial projections avoid quaternion alignment entirely. A 90-degree X rotation is implemented as `[x,y,z] → [x,−z,y]` for every atom and the center, an exact sign/permutation transform that introduces no trigonometric input rounding. Inputs and SUT outputs are frozen before reference comparison in `alignment/rotation/`. This targets the near-parallel and near-antiparallel branches of the pinned glam dependency; it is distinct from the prespecified six-molecule random-rotation campaign and from B1 angular discretization.
