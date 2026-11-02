# Independent source review of buried-volume definitions

This review concerns official, unmodified DBSTEP 2.2.0 commit
`4365353e0a2f24d05336e619759e350bc630a3cb`, the StericX source at
`f04b156f7b8b31700e1c6ab978e337b00cc734e0`, and installed Morfeus 0.8.0.
It does not alter the original failed comparison or approve a replacement gate.

## Convention comparison

The common continuous target is the volume of the evaluation sphere intersected
with the union of selected atomic spheres. A finite lattice additionally needs
an explicit point set, sphere predicate, atom predicate and normalization.
The same continuum target does not imply identical finite estimators.

| Item | StericX | Pinned DBSTEP | Consequence for comparison |
|---|---|---|---|
| Sphere radius | Explicit `sphere_radius`, f32, default 3.5 Å | Explicit Python `r`, f64, default 3.5 Å | Values can be matched; `radius` is not the supported Python keyword |
| Center | Explicit center or virtual center inferred from three donor neighbors | Position of `atom1`, translated to origin; default index 1 | Insert a documented zero-radius Bq at the exact supplied center to match; neither tool independently verifies the chemical center |
| Donor convention | Opposite sum of normalized donor-bond vectors; library distance 2.1 Å, descriptor API 2.28 Å | No analogous default phosphine lone-pair-center construction in ordinary volume mode | Supply the same center; atom-index defaults are not chemically interchangeable |
| Atomic radii | Bondi-style f32 table; H 1.20, Br 1.85 Å; unvalidated unknown fallback 1.80 Å | Named Bondi f64 table; H 1.09, Br 1.83 Å; larger table and fallback paths | Compare explicit common radii separately from different defaults |
| Radius scaling | f32 multiplication, default 1.17 | f64 multiplication, default 1.0; `sambvca=True` sets 1.17 | Match scale and represented radius values explicitly; squaring intermediates can still differ |
| Hydrogens | Occupancy excluded by default, but donor-bound H remains in frame construction | Included by default; `noH=True` excludes occupancy; selected excluded atoms survive as zero-radius Bq | Occupancy and frame identities must be checked separately |
| Metals / ghosts | Ordinary selected atom radii contribute; no real metal required for virtual-center calculation | Metals contribute unless `nometals`; Bq has radius zero and is explicitly skipped | Bq is a reference-center device, not a zero-radius production-atom test for StericX |
| Coordinates / frames | f32 source, translation and projection; negative z points toward donor; three donor-plane orientations | f64; ordinary volume-only call translates but does not rotate; Sterimol/tensor mode can rotate | Public default grids have different orientation/phase; the captured-frame lane explicitly controls this |
| Grid construction | Endpoint-inclusive cube from −R to R; axis count from nominal volume density; f32 arithmetic | Molecular bounding box, extended to cover sphere; requested point spacing; endpoint `linspace`, rounded to eight decimals | Grid extent, phase, axis count and actual spacing differ; nominal point-volume density is not a spacing |
| Sphere membership | `(x*x + y*y) + z*z <= f32(R*R)` when constructing the sphere grid | Tree path: cKDTree radius query; public fast denominator: subtractive residual checks; numerator: sum of squares | The original witness changes the finite sphere set; DBSTEP also has an internal fast-path inconsistency described below |
| Atomic occupancy | Inclusive squared-distance comparison in f32 against rounded squared radius; optimized row reuse preserves operation order | Inclusive cKDTree query or direct f64 sum-of-squares comparison; zero-radius ghosts excluded | Boundaries require explicit point-level diagnosis; identical total counts alone can hide opposing decisions |
| Quadrants / octants | Sign regions, zero belongs to positive side, each region independently normalized; extrema span three planes | No quadrant/octant volume API in this pinned `dbstep/` implementation; residue decomposition is a different quantity | Morfeus/independent sign-region references are needed; do not call residue contributions quadrants |
| Absolute-volume normalization | Occupied point fraction times analytic sphere volume; symmetric mean over three orientations | `bur_vol` is a percentage; `occ_vol` is whole-grid occupied volume, not buried volume; node volumes use nominal spacing³ | An analytic-normalized DBSTEP buried volume is an explicitly derived quantity, not `occ_vol` |
| Boundary / rounding | f32 predicates; zero-plane conventions; symmetric sorted f32 reduction | f64 predicates; eight-decimal grid snapping; public percentage clamped to 100 | Neither percentage identity nor continuum accuracy follows from sharing a named radius table |

Source locations:

* StericX `src/geometry/buried_volume.rs`: defaults at 9–32, center at 482–499
  and 704–765, frames at 849–927, grid at 930–953, occupancy/regions and
  normalization at 955–1062, three-plane summary at 278–354.
* StericX `src/geometry/xyz.rs:351` defines its radius table.
* [DBSTEP `Dbstep.py`][dbstep]: options at 608–672, selected-atom defaults at
  556–586, surface assignment at 247–282, rotation at 284–300, grid at 313–363,
  occupied-point distances at 394–398, and volume calls at 436–447.
* [DBSTEP `sterics.py`][sterics]: fast sphere denominator at 14–26, direct
  occupancy at 29–74, bounds at 83–132, tree occupancy at 135–167, buried volume
  at 383–469, and clipping at 450–455.
* [DBSTEP `constants.py`][constants] and [radii assignment][radii] specify its
  table and fallback behavior. `parse_data.py:251–286` preserves selected
  excluded atoms as zero-radius ghosts.

Morfeus 0.8.0 independently constructs an endpoint-inclusive float64 projection
grid (`geometry.py:358–374`), filters by Euclidean length, and uses one preserved
sphere point set for its occupied fraction (`buried_volume.py:305–335`). Its
regional strict inequalities omit coordinate planes (`buried_volume.py:217–289`),
whereas StericX assigns zero to the positive region. That difference is relevant
on odd grids and does not establish which finite estimator is physically exact.

## A separate DBSTEP public-path inconsistency

The public volume-only fast path computes its denominator using
`R² − x² − y² >= z²`, with successive rounded subtraction. Its occupied numerator
uses `(x² + y²) + z² <= R²`, with rounded addition. These real-arithmetic
equivalences need not produce the same finite point set.

The represented float64 point `(0.6, 0.8, 0)` at R=1 supplies a minimal witness:

* DBSTEP's fast denominator excludes it.
* Its summed numerator predicate and cKDTree sphere query include it.
* Exact rational squared distance exceeds 1 by
  `3602879701896397 / 81129638414606681695789005144064`
  (approximately 4.44e-17). The exact predicate excludes it.

Thus simply declaring all calculations “float64” is insufficient. The two
internal predicates produce different sets on identical represented coordinates.

Twelve actual public calls use a zero-radius Bq at the origin, one carbon of
DBSTEP's unchanged 1.70 Å radius, R=1 Å, and the explicit covering cube
`[-3,3]³`. Four carbon positions are evaluated at spacings 0.2, 0.1 and 0.05 Å.
All settings, axes, captured stdout/stderr and exact-rational counts are in
[the probe results](definition_probe_result.json).

| Carbon position (Å) | Spacing (Å) | Fast sphere denominator | Summed sphere predicate | Occupied numerator | Public % | Same summed-predicate % | Exact-predicate % | Analytical continuum % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| (0,0,0) | 0.2 | 491 | 515 | 515 | 100, clipped | 100 | 100 | 100 |
| (0.8,0,0) | 0.2 | 491 | 515 | 510 | 100, clipped | 99.02912621 | 99.79633401 | 98.42890625 |
| (0,0,1) | 0.2 | 491 | 515 | 461 | 93.89002037 | 89.51456311 | 90.63136456 | 89.048125 |
| (0.8,0,0) | 0.1 | 4,145 | 4,169 | 4,095 | 98.79372738 | 98.22499400 | 98.31121834 | 98.42890625 |
| (0.8,0,0) | 0.05 | 33,377 | 33,401 | 32,882 | 98.51694280 | 98.44615431 | 98.45702130 | 98.42890625 |

For carbon at (0.8,0,0), the coarse public call reports full occupation despite
five unoccupied points in its own numerator's sphere set. Its unbounded ratio is
103.86965377%, then clamped to 100. The exact rational finite-grid occupied count
is 490/491, which also establishes an unoccupied point. The independent analytic
sphere-intersection result is below 100%. This is evidence of inconsistent
numerator/denominator predicates in the reference implementation, not a need to
alter StericX radii or force agreement with the public DBSTEP value.

The inconsistency shrinks in these refinement examples, but it is not erased
from the evidence. Exact rational finite-grid values are not continuum values.
The analytic values use the closed-form two-sphere intersection; π cancels in
the fraction and the represented radii/distances are retained as rationals.
The complete 12 cases include all full/partial and resolution outcomes.

The original ligand-401 strict failure uses DBSTEP's **tree path on supplied
native points**, so this newly identified fast-path inconsistency is not its
cause. The original change is f32 versus f64 sphere membership, already localized
to 32 supplied nodes. Both issues must remain separate in the release decision.

## Implication for the original comparison contract

Calling `occupied` on a supplied finite point set compares atomic occupancy on
that set. Subsequently calling `buried_vol` reselects the sphere in f64; it is no
longer an atomic-occupancy comparison on the original set. A tolerance derived
only for final f32 division/multiplication does not bound jumps in indicator
counts or a changed denominator. The original equality failure remains valid
evidence that the two finite estimators differ, but cannot alone establish a
StericX implementation error.

A defensible revised policy would need separate comparisons of (a) pointwise
occupancy on identical supplied points, (b) each tool's declared sphere/grid
predicate, and (c) independently bounded continuum estimates. Every ambiguous
boundary decision should retain its exact margin; a new tolerance must not be
chosen from the observed 0.024 Å³ residual. The parent investigation must document
and approve any such policy before using it for acceptance. This source review
does not declare that approval or waive any failed gate.

## Reproduction and provenance

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .stericx/paton/reference/venv/bin/python \
  docs/validation/dbstep/resolution_analysis/probe_dbstep_definitions.py \
  --output .stericx/dbstep_resolution/definition_probes/NEW_RUN
```

The script refuses an existing output directory, verifies the official commit
and clean checkout, and checks installed `Dbstep.py`/`sterics.py` bytes against
that checkout. The current scientific result is
`.stericx/dbstep_resolution/definition_probes/run_v2/result.json`. The initial
run completed the same numerical experiment before a documentation-string lint
fix; its output and exact controller remain under `run_v1/`. No production source,
reference source, original evidence or numerical acceptance limit was edited.

[dbstep]: https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/Dbstep.py
[sterics]: https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/sterics.py
[constants]: https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/constants.py
[radii]: https://github.com/patonlab/DBSTEP/blob/4365353e0a2f24d05336e619759e350bc630a3cb/dbstep/radii.py
