# Independent DBSTEP reference

Reference: the official [patonlab/DBSTEP repository](https://github.com/patonlab/DBSTEP),
version **2.2.0**, pinned commit
[`4365353e0a2f24d05336e619759e350bc630a3cb`](https://github.com/patonlab/DBSTEP/tree/4365353e0a2f24d05336e619759e350bc630a3cb).
No DBSTEP code or scientific constants are modified. The environment freezer
requires a clean checkout and compares every installed DBSTEP Python source
file byte for byte with that checkout. It also hashes all 757 tracked reference
files, including the official example geometries, density cubes and tests.

## Reproduction

```bash
git clone https://github.com/patonlab/DBSTEP.git .stericx/paton/reference/DBSTEP
git -C .stericx/paton/reference/DBSTEP checkout 4365353e0a2f24d05336e619759e350bc630a3cb
uv venv --python 3.12.13 .stericx/paton/reference/venv
uv pip install --python .stericx/paton/reference/venv/bin/python \
  -r docs/validation/dbstep/requirements.lock
.stericx/paton/reference/venv/bin/python scripts/freeze_dbstep_reference.py \
  --checkout .stericx/paton/reference/DBSTEP --output NEW_RECEIPT_DIRECTORY
cd .stericx/paton/reference/DBSTEP
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ../venv/bin/python -m pytest -q tests
```

The original setup installed the clean local checkout with `uv pip install`,
then Morfeus 0.8.0 and pytest; the portable lock uses the same official VCS commit.
[environment.json](environment.json) records Python 3.12.13, platform, installed
package versions, direct installation URL, interpreter hash, source hashes and
the successful source-equivalence check. [requirements.lock](requirements.lock)
pins all 14 installed distributions. `freeze.log` and `upstream_tests.log`
preserve execution. Freezing refuses to overwrite an existing receipt.

The unmodified upstream suite completed with **413 passed and 3 skipped**.
The skips are the two optional RDKit molecule tests and the optional PyMOL
module; neither dependency is installed in this reference environment. The
checkout remained clean after testing. All included classic, cube, parser,
volume, ensemble and coordinate tests passed. This is an installation check,
not proof of agreement with StericX.

## Inputs and comparison scope

The existing retained StericX cohort is the ordered 11-ligand/56-conformer set
sealed in [the sensitivity manifest](../../paton_recommendations/conformers/sensitivity_v1/manifest.json).
That manifest records every original XYZ hash, ordered energy metadata and
descriptor setting. Official DBSTEP fixture inputs, including
`tests/cube_files/Ne_medium.cube`, `benzene_coarse.cube` and `CH2CMe3_100.cube`,
are independently hashed in this environment receipt. Exact subsets, transformed
geometries, reference atom mapping, equations and numerical acceptance criteria
must be frozen in the cross-validation protocol **before** running comparisons.
Upstream tests establish that the reference installation operates; they are not
StericX cross-validation results.

Relevant settings to record explicitly are surface (`vdw` or `density`),
measurement (`classic` or `grid`), supplied or named radii, scaling, H/metal
inclusion, atom indices, center, rotation, sphere radius, grid spacing or density,
isodensity threshold and any scan/shell settings. DBSTEP indices are one-based;
StericX public indices are zero-based. A volumetric point spacing in Å is not
the same parameter as StericX's volume-per-point density in Å³.

## Conventions identified from the pinned source

* DBSTEP's Bondi H radius is 1.09 Å. StericX's existing steric radius convention
  uses 1.20 Å. Preserve both defaults and label any comparison with explicitly
  supplied common radii. Do not alter either table to reduce residuals.
* DBSTEP classic L evaluates `max(abs(z_i) + r_i)` in its aligned frame, despite
  a source comment describing positive-axis extent. StericX uses signed
  `max(z_i + r_i)` and its coordination API adds 0.40 Å. These definitions can
  coincide for forward geometries, but need not coincide when atoms project
  behind the reference plane. Preserve counterexamples and classify them.
* DBSTEP classic B1 scans 361 endpoint-inclusive angular positions; the native
  default uses 360 unique positions. Frame phase and floating arithmetic must
  be distinguished from differences in the continuous support function.
* DBSTEP's `occ_vol` denotes the whole occupied molecular grid volume. It is
  **not** the occupied volume inside the evaluation sphere. Its buried-volume
  percentage uses occupied sphere points divided by total sphere points.
  StericX's absolute buried volume multiplies its sampled fraction by the
  analytic sphere volume. Compare common definitions explicitly.
* DBSTEP volume scans change sphere radius at a fixed center. They do not
  translate the center along a bond. Density mode uses a strict `rho > isoval`
  occupation predicate. The pinned default is 0.0016; older/cached documentation
  mentions 0.002. Any validated density pathway must require an explicit
  threshold and preserve the actual grid and units.

These are source-based observations, not results of fitting parameters or
modifying the reference. Density-grid compatibility and physical interpretation
will be established separately before implementation is admitted.
