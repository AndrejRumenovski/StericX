# Experimental structure and descriptor reproduction

## Inputs and scope

We independently retrieved **135 original conformer exports for eleven ligands** in
Cadge et al.'s Ni-catalyzed homo-Diels–Alder selectivity screen. All original records,
identifiers, atom ordering and unsuccessful requests are retained. The twelve
quantified experimental observations include ligand 906, whose metadata request
returned HTTP 500 and whose row is absent from the authors' descriptor CSV; it is
excluded explicitly. Unquantified L9 is not assigned an ee. The primary outcome is
the SI integer ee at 353.15 K, with author ΔΔG preserved separately. The corrected
SI's Table S3 is unchanged. See [dataset review](DATASET_REVIEW.md).

The [Kraken publication](https://doi.org/10.1021/jacs.1c09718) describes searches of
free and Ni(CO)3-complexed ligands followed by selection and free-ligand DFT
optimization. Our inputs are the MolSSI Kraken conformer exports, not newly
optimized structures. The API did not expose calculation logs or usable energies
for these specific exports: all 135 valid `/conformers/data/{id}` responses contained
`data: null`. The earlier guessed endpoint returned 404 and is retained as a failed
acquisition attempt. Library-level provenance does **not** verify an individual
export's energy, optimization level or thermal population. “Published-library” is
the operative distinction throughout this study.

`prepare_chemistry.py` validates neutral formal charge, one connected fragment,
one trivalent P with exactly three source-identified neighbors, explicit hydrogen
counts, connectivity against source SMILES and identical ordered graphs within each
ligand. No metal atoms are present. The virtual coordination center is a geometric
convention, not an experimentally determined Ni coordinate. Atom mapping is identity
within each retained export; full coordinates, bonds and frame axes are serialized
in native observations. V2000 exports round coordinates to four decimal places;
the engine then stores coordinates/radii in f32. These uncertainties are distinct
from continuous-integration error.

## Stereochemistry findings

Ordinary RDKit 3D assignment exactly reproduced the encoded source labels for 101
exports and left P partly unassigned for 34. This was **not** accepted as proof of
matching stereochemistry. `audit_stereochemistry.py` compares ordered P-neighbor
scalar triple products with ten independently embedded, source-stereo-constrained
geometries per explicitly labeled P center. A virtual fourth-point diagnostic,
with the correct implicit-neighbor parity conversion after removal, agrees.

| Ligand | P comparison against encoded source SMILES | Consequence |
|---|---|---|
| 401 | All 15 opposite; fixed carbon configurations agree | Source/export diastereomer disagreement; cannot silently treat as verified experimental ligand geometry |
| 1057 | All five match | Initial virtual-point interpretation without parity conversion was rejected |
| 1058 | Eight match, six opposite; carbon configuration agrees | Distributed “ensemble” mixes configurations |
| Other eight | Source P unassigned or not stereogenic | No unsupported absolute P assignment claimed |

Relative BINOL-axis inspection must use chemically appropriate aromaticity: default
RDKit aromaticity marks the P-containing ring aromatic and hides its single biaryl
bond. With original SDF bond orders and MDL aromaticity, the indexed axial sector is
constant within 785, 2062, 2064 and 2067; 2063 has eleven in one sector and one in the
other. Spiro ligand 498 has a consistent indexed handedness. These are relative
geometric checks, not experimentally verified absolute atropisomer assignments.
The primary scalar extrema and absolute ee are invariant to reflecting an entire
molecule; flipping only P in a ligand with fixed carbon stereocenters is different.
Oriented regional labels cannot be pooled across reflections without a mapping.

The original eleven-ligand analysis remains the **as-distributed library association**.
We also report exclusion and filtered-structure sensitivities. The first P amendment
was recorded before fitting new models; the later axial/aromaticity sensitivity is
explicitly post-primary. No original measurement or structure is overwritten.

## Authoritative Rust calculation

`observer/src/bin/chemistry.rs` reads each original SDF with `Ensemble::from_file`,
assigns its source ID, calls native descriptor/statistics APIs and records every
success/error. All 135 published conformers and all eleven ensembles succeed.

| Quantity | Frozen convention |
|---|---|
| Sterimol L/B1/B5 | P–virtual-metal coordination axis; all real atoms including H; native Bondi-style radii; default 360 angular samples; +0.40 Å native L correction |
| Virtual metal | 2.28 Å from P, opposite the normalized three-neighbor sum; native f32 frame conventions |
| Buried volume | Integration R=3.5 Å, density=0.01 Å³, radius scale=1.17, H excluded |
| Frames | Three P-neighbor planes, with original plane-atom identity preserved |
| Regions | Quadrants (++),(-+),(--),(+-); positive-Z octants first, negative-Z last; near is negative local Z |
| Summaries | Total volume mean over frames; extrema over frames; maximum adjacent-quadrant difference in Å³, not percent |
| Pyramidalization | Existing native P-distance and angular descriptors; no new empirical correction |
| Ensembles | First numeric conformer ID, minimum, maximum, mean, population SD and range; the first ID is not called lowest-energy |
| Published weighting | None: energies and a common energy model are unavailable |

Regional raw volumes and percentages are separately defined. The finite-grid regional
volumes use their own sampled denominators. New continuous region percentages use
exact geometric region volume (one quarter/eighth of the integration ball), with
outward-rounded denominator bounds. This convention distinction is preserved.

For every published conformer, continuous support-envelope Sterimol and all three
certified regional frames were additionally computed. **405/405** regional requests
met 0.05 Å³ absolute error with a 1,048,576-cell per-frame budget. Certificates apply
to the explicitly serialized local balls after the established f32 frame transform;
they do not certify source coordinates, atomic-radius choices or electronic structure.

## Independent checks and practical precision

Python independently reproduces all eleven-feature arithmetic reductions and all
complete generated-energy weighted means (relative/absolute check 2e-12 or tighter).
Upstream Morfeus and DBSTEP independently evaluate all **135** matched geometries
and radii. B1/B5 comparisons use their mathematically derived finite-angular-sampling
bounds; an external sampled answer is not asserted to be an exact continuous extremum.
Morfeus donor-referenced L and DBSTEP exact radial B5 agree within a scale-derived
binary64 arithmetic margin. DBSTEP's different absolute-axis L is retained separately.
All 135 comparisons pass. Independent SciPy centered-covariance regression reproduces
all sixteen primary univariate outer-fold fits to within 1e-11 kcal/mol.

Across the published conformers, maximum native-versus-continuous discrepancies are:
L 1.092e-6 Å, B1 0.04223 Å, B5 4.756e-7 Å, adjacent-quadrant difference 0.31055 Å³,
and any individual quadrant volume 0.21485 Å³. These are **observed finite-grid/angle
errors on these geometries**, not universal bounds. Continuous regional errors remain
separately enclosed. Using certified regional midpoints in the same models changes
single/ensemble RMSE from 0.84977/0.38941 to 0.84988/0.39548 kcal/mol: the numerical
change is much smaller than the statistical uncertainty of the proposed improvement.

## Reproduction artifacts

[Evidence bundles](data/README.md) preserve original and generated inputs, unsuccessful
preparations, native observations, provenance and analyses. Human-readable derived
CSVs/JSON are in `results/`. Run instructions are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
Raw timing/counters and full calculation records remain under
`.stericx/experimental_validation/`. No generated conformation replaces a published
geometry, and no missing source energy is inferred from conformer IDs or aggregate
reference descriptors.
