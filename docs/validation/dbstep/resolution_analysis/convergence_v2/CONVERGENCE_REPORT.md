# Fixed-resolution buried-volume convergence

This study preserves every scheduled comparison and does not itself score the approved release policy. The Rust production engine and all reference implementations remain unchanged.

The accepted campaign contains **20 geometries, 320 comparisons**, and 0 unsuccessful comparisons. It includes 17 mathematical primitives and three captured native frames of original ligand 401/0 at R=2.5 Å.

## Methods and provenance

Native and Morfeus density controls are 0.08, 0.01, 0.00125 and 0.00015625 Å³. DBSTEP spacings are 0.4, 0.2, 0.1 and 0.05 Å. These controls generate different lattices. All comparisons use identical accepted binary32 atomic coordinates and scaled radii promoted exactly to the reference precision; the native kernel retains its actual rounded squared radii.

Morfeus is called through public `BuriedVolume`, with an excluded zero-radius ghost center at the origin, explicit radii, `radii_scale=1`, and hydrogen inclusion after shared explicit filtering. DBSTEP has no equivalent public custom-radii input: its unmodified `occupied`/`buried_vol` and `occupied_direct`/`buried_vol` routines are therefore identified as support-function lanes. Their actual rounded axes cover an explicit cube [-3,3]³ Å³. Separate unchanged public DBSTEP calls and their reference defect are documented in [DBSTEP_SOURCE_REVIEW.md](../DBSTEP_SOURCE_REVIEW.md).

The [actual-cell derivation](../CONVERGENCE_CELL_BOUNDS.md) includes finite coverage, nonuniform axis widths, direct-reference roundoff, observed mask differences, count/weight conversion, ratio normalization and clamp arithmetic. DBSTEP fast denominator masks are explicitly reconstructed and count-checked, not independently observed point decisions. The bounds are conservative and often broad; they do not establish fine numerical precision.

All raw native requests/responses, reference-returned points, axes, full masks, disputed-node witnesses, source snapshots, dependency hashes and elapsed costs are bound by the campaign manifest and [receipt](convergence_receipt.json). Native timings include observation calls and JSON transfer and are **not production latency benchmarks**.

## Independent reference checks

Adaptive boxes use the prospectively fixed depths 8 (primitives) and 11 (original molecular frames). Analytical box-enclosure failures: 0; non-nested interval failures: 0; analytical cell-bound failures: 0.

Across convergence grids, 801 native disputed-node observations retain exact rational margins and explicit binary32 operation traces; unexplained native disputed nodes: 0. The separate 10-point boundary set contains 3 native/exact disagreements, with 0 unexplained. These witnesses distinguish rounded point inclusion from the continuous integral.

### Original ligand 401 continuous bounds

| Frame | Lower volume (Å³) | Upper volume (Å³) | Interval width (Å³) | Reference cost (s) |
|---|---:|---:|---:|---:|
| original401_plane15 | 16.475986864 | 16.609940067 | 0.133953203 | 67.544 |
| original401_plane13 | 16.475893615 | 16.610025123 | 0.134131507 | 45.774 |
| original401_plane10 | 16.475974859 | 16.609943748 | 0.133968890 | 46.690 |

These intervals concern the continuous volume of the accepted balls, including each frame's stored coordinate rounding. Their widths must be respected when interpreting the historical ~0.024 Å³ difference. Full convergence rows and signed/absolute/relative error intervals are in [convergence_residuals.csv](convergence_residuals.csv); [convergence_results.json](convergence_results.json) retains settings, counts and bound components.

| Frame | Native at density0.01 (Å³) | Proven absolute-error lower bound (Å³) | Native at finest density (Å³) |
|---|---:|---:|---:|
| original401_plane15 | 16.418327332 | 0.057659532 | 16.525009155 |
| original401_plane13 | 16.441848755 | 0.034044860 | 16.525650024 |
| original401_plane10 | 16.430088043 | 0.045886815 | 16.525812149 |

A positive error lower bound proves that the corresponding default-density estimate differs from the continuous integral. Finite-grid correctness does not imply continuum exactness or the historical arithmetic-only 0.00025 Å³ accuracy. The finite-width molecular intervals do not support ranking differences substantially smaller than their widths.

### Maximum analytical error by level

These maxima summarize the fixed synthetic analytical cases; they are not a universal accuracy claim or a molecule benchmark.

| Method | Level | Maximum absolute error upper bound (Å³) | Maximum conservative cell error bound (Å³) |
|---|---:|---:|---:|
| dbstep_fast_common_radii | 0 | 0.744061418 | 249.172501 |
| dbstep_fast_common_radii | 1 | 0.417661942 | 124.194470 |
| dbstep_fast_common_radii | 2 | 0.049944807 | 62.154484 |
| dbstep_fast_common_radii | 3 | 0.019857462 | 31.030638 |
| dbstep_tree_common_radii | 0 | 0.744061418 | 249.172501 |
| dbstep_tree_common_radii | 1 | 0.417661942 | 124.194470 |
| dbstep_tree_common_radii | 2 | 0.052446775 | 62.048262 |
| dbstep_tree_common_radii | 3 | 0.020504004 | 31.018010 |
| morfeus_public_common_radii | 0 | 1.265363708 | 283.496432 |
| morfeus_public_common_radii | 1 | 0.108006922 | 141.267227 |
| morfeus_public_common_radii | 2 | 0.090199789 | 68.946768 |
| morfeus_public_common_radii | 3 | 0.013682266 | 33.716915 |
| stericx_native | 0 | 1.265364287 | 283.496527 |
| stericx_native | 1 | 0.106319160 | 141.925168 |
| stericx_native | 2 | 0.090200064 | 68.947460 |
| stericx_native | 3 | 0.013579729 | 33.719505 |

There are 31 rigorously increasing analytical-error steps among the finer successive grids. Every such step is retained in [convergence_trends.json](convergence_trends.json). Center-point integration is not monotonic; finer grids are not automatically closer for every geometry.

## Boundary/convention findings

For exact external tangency (evaluation sphere R=2.5 Å, atom r=1 Å at x=3.5 Å), the continuous intersection is zero. Some native/DBSTEP grids include one closed-boundary contact node and consequently report a small nonzero finite-grid estimate. Morfeus's public prefilter uses `atom.radius + sphere.radius > norm(atom.coordinates)` and removes the exactly tangent atom, reporting zero. This difference first occurs in Morfeus's prefilter, before point occupancy; it is a finite-grid boundary convention with the same zero continuous-volume limit. It does not justify changing validated Rust defaults.

The original strict 401 comparison failure remains a separate finite-domain reclip issue. This convergence campaign neither erases its failed historical criterion nor fits a tolerance to that residual. No production correction was introduced.

## Remaining limits

The shape tests are mathematical inputs, not characterized experimental molecules. Density and spacing are distinct numerical controls. Continuous box intervals remain finite-width; shell bounds overcount overlapping boundaries and can be much larger than observed analytical errors. Region summaries and molecular rigid motions are evaluated in the separate point-domain/adversarial campaigns. Cross-platform installation and matched performance are separate release requirements.
