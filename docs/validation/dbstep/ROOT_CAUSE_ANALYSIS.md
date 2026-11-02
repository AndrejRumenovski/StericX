# Hard-sphere volume discrepancy: root-cause investigation

Investigation starts from StericX commit
`f04b156f7b8b31700e1c6ab978e337b00cc734e0`. This initial record freezes the original
failure and comparison definitions before any production correction. The original
strict gate remains failed. New independent analysis, convergence and adversarial
results will be appended with separate receipts; this document is not yet a
release admission.

## Original failure inventory

The remaining numerical failure is **SIG-NIHDA-401, conformer index 0**, with
sampling radius **2.5 Å**, nominal density **0.01 Å³**, radii scale **1.17**,
hydrogen occupancy excluded, and the recorded explicit center and donor/neighbor
identities. All three donor-plane orientations fail. The same witness reappears
in the radius-profile comparison. It is one geometry/configuration, not several
independent chemical failures.

The [original request](resolution_baseline/original_failure/requests.jsonl) and
[original observer response](resolution_baseline/original_failure/native_stdout.jsonl)
are verbatim copies, including accepted coordinates, atomic radii, aligned atoms,
grid points and individual bins. [Source-conformer metadata](resolution_baseline/original_failure/source_conformer.json)
preserves source geometry and identity. SHA-256 identities, all parameters and
the source-bound executable are in the [baseline manifest](resolution_baseline/manifest.json).

Differences below are absolute; relative differences divide by the absolute
DBSTEP value. The volume comparison uses the DBSTEP sampled fraction multiplied
by the analytic sphere volume, not its whole-molecule occupied-grid volume.

| Result | StericX | DBSTEP derived common normalization | Absolute difference | Relative difference | Original tolerance |
|---|---:|---:|---:|---:|---:|
| Plane atom 15, Å³ | 16.4183273315 | 16.4423074752 | 0.0239801436 | 0.145844% | 0.00024967135 Å³ |
| Plane atom 13, Å³ | 16.4418487549 | 16.4659654716 | 0.0241167167 | 0.146464% | 0.00024967135 Å³ |
| Plane atom 10, Å³ | 16.4300880432 | 16.4541364734 | 0.0240484302 | 0.146154% | 0.00024967135 Å³ |
| Three-plane mean, Å³ | 16.4300880432 | 16.4541364734 | 0.0240484302 | 0.146154% | 0.00024967135 Å³ |
| Mean buried percentage | 25.1033248901 | 25.1400686788 | 0.0367437887 percentage points | 0.146156% | 0.00038146973 percentage points |
| Sphere points, each frame | 5,565 | 5,533 | 32 | 0.578348% | Exact equality |

The original volume tolerance is `32 × 2^-23 × 4πR³/3`; the percentage tolerance
is `32 × 2^-23 × 100`. They were specified before the comparison and remain
unchanged in its historical lane. The other 68 cases / 204 frames passed them.
The [machine-readable failure ledger](resolution_baseline/original_failures.json)
preserves every failed row and denominator diagnostic, including repeated profile
results and first-attempt failures.

Earlier failed execution attempts remain separate:

* `volume_v1` supplied unsupported `radius` instead of official `r` for R=2.5
  and 4.5 Å. DBSTEP ignored it and kept its default; the controller then rejected
  the inconsistent reconstructed counts. This was a confirmed study-controller
  error, corrected in `volume_v2` without changing its tolerance. The original
  outputs and controller remain frozen.
* The first observer build lacked a direct serde dependency. It produced no
  scientific comparison. The corrected build preserved locked scientific crates.
* The first density controller failed to serialize a NumPy boolean after its
  calculations. Its corrected successor passed the original independent criteria.
* Sterimol's six backward-axis convention records, default H-radius offsets and
  finite-angle B1 residuals are retained labelled comparisons, not additional
  unexplained strict volume failures. All 61 conditional Sterimol cases and all
  12 density comparisons passed their separately declared criteria.

## Implementations and evidence frozen

The fresh baseline native executable has SHA-256
`3b9079338386a2532b17ebe844a28d72344a7bdf91d83526d4b3d0fd00438424`, identical to
the executable used in the prior confirmed benchmark. Its new observer hash is
`6b3a5ccf0c8d95275f4c32b3a96207fc4beaecab9eeab054d257161abf24dcef`.
The build snapshots the Rust sources, locks, helper programs, commands and
toolchain. The observer only exposes native state; it is not an independent
reference implementation.

Official DBSTEP **2.2.0**, commit
`4365353e0a2f24d05336e619759e350bc630a3cb`, remains unmodified. Its fresh
[environment receipt](resolution_baseline/dbstep/environment.json) verifies a
clean checkout and byte-identical installed Python sources, with all tracked
reference-file hashes and installed dependency versions. Python is 3.12.13;
Rust/Cargo are 1.97.0. Morfeus 0.8.0 and the scientific dependency versions are
recorded in that receipt and its lock. Original reports and inputs are hashed
before new evaluation, not overwritten.

## Scientific definitions and implementation comparison

The continuum target is the volume of the evaluation sphere intersected with
the union of the selected atomic spheres. Its percentage divides by `4πR³/3`.
Finite grids approximate this integral; boundaries have zero continuum measure
but nonzero weight in a sampled estimator. [SambVca's official methodology](https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html)
describes voxel classification and the occupied/sampled-sphere ratio. Its
physical definition does not require all programs to choose identical nodes.

| Item | StericX | Official pinned DBSTEP | Original comparison |
|---|---|---|---|
| Sphere radius | Explicit f32 R; default 3.5 Å | Explicit `r`; default 3.5 Å | Same accepted R in corrected controller |
| Reference center | Explicit center or geometric donor construction | Selected atom position; virtual `Bq` can supply center | Identical captured center; no center fitting |
| Atomic radii | Existing Bondi-style table with documented extensions | Official Bondi table; H=1.09 Å versus native 1.20 Å | Captured common radii in finite-point lane; native-table differences labelled in public-API lane |
| Scaling | Default 1.17, f32 product | Default 1.0, explicit `scalevdw`, f64 product | Explicit scale matched; product representation retained |
| Hydrogen | Occupancy excluded by default; bonded H remains in frame | Included by default; `noH=True` excludes it | Explicit policy matched; H-included default-table lane has distinct surfaces |
| Donor/metal | Three explicit donor neighbors define native orientations; virtual center is not an atom | No three-neighbor requirement for volume; explicit metal/exclusion controls | Same retained atoms; zero-radius Bq adds no occupation |
| Coordinate transform | Three donor-plane frames; total is their symmetric mean | Volume-only public API translates to center, without native three-frame reduction | Supplied aligned atoms identical in finite-point lane; source-frame public grid is distinct |
| Grid construction | f32 endpoint-inclusive cube; side count from nominal point volume; then sphere mask | f64 molecular bounding box and explicit spacing; coordinates rounded to eight decimals | Captured native nodes supplied to support functions; separately generated public grids retained |
| Sphere membership | Inclusive f32 squared-distance comparison during grid construction | Support path uses cKDTree denominator and summed-square numerator; public fast path uses subtractive denominator | A second DBSTEP boundary test changes the supplied sample set at R=2.5 |
| Atomic membership | Inclusive f32 squared distance and rounded radius square | cKDTree or direct f64 inclusive distance test | Full point decisions must be compared, not only counts |
| Quadrants/octants | Zero planes assigned positive; separately normalized regions and three-frame reductions | No equivalent native three-frame regional summary is supplied by the compared volume API | Morfeus regions require explicit sign mapping and zero-plane convention |
| Precision | f32 geometry/grid/predicates and reported values | f64 grid/predicates, eight-decimal axis snapping in public path | Precision and node selection differ even when decimal inputs match |
| Absolute normalization | Occupied fraction × analytic sphere volume | Buried percentage from counted volumes; `occ_vol` is whole-grid molecular volume | Compare explicitly derived same-normalization quantities only |

The pinned [official source](https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/sterics.py)
and `Dbstep.py` establish these executable conventions. [Morfeus documentation](https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html)
supplies a third documented implementation. Original method references are
[Poater et al. (2009)](https://doi.org/10.1002/ejic.200801160) and
[Falivene et al. (2016)](https://doi.org/10.1021/acs.organomet.6b00371); current
publisher retrieval was unavailable, so source equations and accessible official
manuals provide the directly inspected implementation evidence. No inaccessible
paper is claimed to have been read in full.

## First identified mathematical divergence

The prior frozen exact-distance analysis finds 32 captured native nodes whose
exact squared distances exceed 6.25 Å² by approximately 1.08e-7 to 2.17e-7 Å².
The f32 sum rounds them onto the inclusive sphere boundary. The f64 reclip
removes them. Native occupied counts are 1,396 / 1,398 / 1,397; DBSTEP atom
occupation on the supplied nodes gives the same counts. Six occupied nodes per
frame disappear during the subsequent sphere test, producing 1,390 / 1,392 /
1,391 over the smaller denominator.

Thus the historical comparison described as a matched finite lattice actually
compares two different sample domains in this witness. Its arithmetic-only
tolerance cannot by itself bound that change. This identifies a potential
**reference comparison error**, alongside the demonstrated **numerical
integration difference**; it does not yet authorize changing a release gate.

## Independent reproduction and classification

The new observation adapter calls the unchanged Rust kernel for each node and
separately for the complete grid. Its production source prefix and legacy
responses match the frozen baseline exactly. It is a measurement of the engine,
not the independent reference. The independent research implementation uses
direct union predicates, exact rational boundary witnesses, analytical sphere
intersections and adaptive continuous-volume bounds.

Across **69 cases / 207 frames / 3,545,298 sampled occupancy decisions**, Rust,
unmodified DBSTEP occupancy and the independent reference give **exactly equal
atomic masks**. Counts reconstruct all native aggregate/quadrant/octant bits.
The maximum common-domain volume residual across all frames is
**7.86485133e-6 Å³**, within the unchanged original arithmetic tolerance in every
case. All **138** public mean-volume/percentage comparisons also satisfy their
unchanged arithmetic criteria.

| Original witness | Rust, unchanged Å³ | DBSTEP same-domain derived Å³ | Residual Å³ |
|---|---:|---:|---:|
| Plane 15 | 16.4183273315 | 16.4183263867 | 9.44858e-7 |
| Plane 13 | 16.4418487549 | 16.4418483443 | 4.10624e-7 |
| Plane 10 | 16.4300880432 | 16.4300873655 | 6.77741e-7 |
| Three-frame mean | 16.4300880432 | 16.4300873655 | 6.77741e-7 |

The percentage comparison is 25.1033248901% native versus 25.1033243486% on the
same domain, a 5.41529e-7 percentage-point residual. These values use unchanged
geometries/radii and normalization, without a fitted integration parameter.

**Classification of every original strict failure:** reference comparison error
(the purported common finite domain was clipped again), together with a
finite-precision numerical integration difference (binary32 versus binary64
sphere membership). The first divergent operation is the second sphere test,
not atom mapping, coordinate transformation or atomic occupancy. The original
three frame failures and two public-summary failures, and their repeated profile
witness, remain reproduced in the historical lane. No production Rust error is
established by these results. The smallest justified correction is the separate
common-domain comparison implemented in `analyze_dbstep_point_domains.py`.

Full evidence and reconstruction commands are in
[POINT_DOMAIN_RESULTS.md](resolution_analysis/POINT_DOMAIN_RESULTS.md). The
[molecular adversary study](resolution_analysis/MOLECULAR_ADVERSARIES.md) adds
49 requests, explicit hydrogens, primary/secondary phosphines, permutations,
rotations, translations and nearly degenerate frames. All 48 mathematically
valid requests succeed; the exactly collinear frame returns a structured error.
An initially incorrect test expectation for a small nonzero plane was corrected
using its exact nonzero cross product, with the unsuccessful attempt retained.

## Separate DBSTEP reference implementation defect

Official DBSTEP's volume-only fast path counts sphere denominator points using
successive subtraction of coordinate squares, while its occupied numerator uses
a summed-square comparison. Those floating-point predicates can disagree.
Independent exact rational witnesses and analytical sphere intersections confirm
that this can produce more occupied points than denominator points; its public
API then clamps the percentage to 100%.

For the frozen Bq-center/carbon test at R=1 Å, carbon center (0.8,0,0), Bondi
radius 1.7 Å and spacing 0.2 Å, the unmodified public API counts 510 occupied
nodes over 491 denominator nodes. Its unclamped ratio is 103.86965%, clamped to
100%, despite five unoccupied nodes in its own 515-point summed-square sphere.
The analytical continuum result is 98.42890625%. The exact rational finite-grid
ratio is a different quantity, 490/491, or 99.79633401%.

This is a **confirmed DBSTEP internal predicate inconsistency**. It is separate
from the original ligand-401 failure, which used the tree path. The source is
unmodified; this defective public fraction is not used as an exact finite-set
oracle. [Source review and all 12 public probes](resolution_analysis/DBSTEP_SOURCE_REVIEW.md)
preserve the settings, counts, alternative finite definitions and continuum truth.

## Release-policy status

The [separate-lane policy](ACCEPTANCE_POLICY.md) was explicitly
[approved by the user](POLICY_APPROVAL.md) before final release-gate evaluation.
It preserves the historical failed equality lane and original tolerances, and
requires independent boundary and continuum analysis plus complete scientific,
engineering and performance evidence. Approval of the policy is not a release
pass. Convergence results and the final source-bound gate decision are reported
separately when complete.
