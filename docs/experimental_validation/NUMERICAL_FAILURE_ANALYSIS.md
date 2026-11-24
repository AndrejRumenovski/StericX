# Investigation of the eighteen historical failed requests

The locked baseline is `occupancy_precision/c1_validation.json`, not a newly
selected set of difficult cases. Requests and old outputs are preserved, including
every original tolerance and budget. The requests are repeated configurations of
four mathematical families; their denominator remains eighteen.

| IDs / variants | Count | Target Å³ | Original budgets | Mechanism |
|---|---:|---:|---|---|
| three_overlapping_chain, 1023/8191/65535; rigid; permuted | 5 | 1e-5 | 1,023 / 8,191 / 65,535 cells | Smooth chord patches cease at projected circle boundaries and uncertain endpoint order. Largest original bound 0.929395, smallest 0.00772958; all exhaust the budget. The first and third unit spheres touch, so positive-volume triple overlap is absent. |
| three_disjoint, same variants | 5 | 1e-5 | same | Disconnected spheres still trigger general integration. Original bounds 0.186916 to 0.00220648; all exhaust the budget. No overlap correction is necessary. |
| boundary_1e-200, overlap/external_below/equal/above | 4 | minimum positive f64, 2^-1074 | 1,023 | Radius squares underflow. Repeated outward steps manufacture several subnormal ULPs of enclosure slack; analytic dispatch returns floating_point_limit. Contact/noncontact cannot be reliably separated with the old squared-distance interval. |
| boundary_1e-150, same variants | 4 | 2^-1074 | 1,023 | Volume underflows while radius square remains represented. Old bound already [0,2^-1074], but unconditional upward rounding of an exactly represented endpoint difference makes error=2*2^-1074 and prevents success. |

The tiny test spheres are mathematical arithmetic stresses, not molecular
coordinates. Passing their absolute tolerance does not establish relative accuracy:
relative error is undefined when the certified lower volume is zero. No radius or
numerical target is changed.

## G1: explicit geometric proof before general integration

New `sphere_union_volume_certified_geometric` is separately selectable; the old
`sphere_union_volume_certified` remains untouched. Canonical sorting/deduplication
makes permutation behavior deterministic. Check zero-measure separation using exact
axis differences and radius sums, bounded by error-free TwoSum. If these are
inconclusive, use the conservative Euclidean interval. Never infer disjointness
from a rounded squared distance alone.

If every atomic ball is certified inside the integration ball, construct a graph
containing every possible positive-volume overlap. Reject the shortcut if this
graph has a triangle. A triple intersection implies three pair intersections,
so a triangle-free conservative graph proves all triple/higher intersections have
zero volume. Then V(union)=ΣV(ball)-ΣV(pair lens) is exact as a real identity.
Interval-enclose every term, accumulation and subtraction. The graph may reject
valid shortcuts; that loses performance, not correctness. Clipped balls and possible
triple intersections fall back to the original certified method. This is not a
general pairwise approximation.

TwoSum supplies the exact residual of finite nonoverflowing addition in IEEE
round-to-nearest with gradual underflow; its sign identifies the one adjacent
float needed for an enclosure. See [Shewchuk's arithmetic development](https://www.cs.cmu.edu/~glmiller/Publications/ShewPHD.pdf).
6,005 independently generated Python Fraction witnesses test both interval endpoints,
including cancellation, exponent disparity and subnormals. Exact error-distance
subtraction avoids adding an unnecessary ULP when the subtraction is exact.

Isolated G1 resolves **14/18** original requests (ten unions and four exact/external
contacts), with all 167 original independent reference intervals contained. Four
positive subnormal-volume cases remain. The 180-observation extended suite (96
independent analytical unions and 84 molecular convergence requests) agrees with
the initial combined prototype after removing time/allocation diagnostics. No
unresolved contradictory enclosure is suppressed. Budget-limited fallbacks remain
available and independently tested.

## Evidence and interpretation

`.stericx/experimental_validation/geometric_graph_v1/` contains native output,
validation receipt, extended output, API/exact arithmetic test logs, wall/CPU time
and peak RSS. Timings are exploratory under the ordinary shared workstation load;
no controlled speed ratio is claimed. `geometric_v1/` retains the earlier combined
prototype and its tests. The combined prototype is decomposed into independently
reviewable accepted changes rather than retrospectively treated as one experiment.
Further subnormal improvement and regional methods are evaluated separately.

## G2: global bound for positive subnormal volumes

For 0<R<=2^-359 Å, V<=4πR³/3<8R³<=2^-1074 Å³, using π<6. The radius
threshold is constructed as an exact exponent bit pattern. After validated exact
nonoverlap checks, return [0,2^-1074] if no sharper proof exists. The estimate is
zero, but the nonzero upper/error bound explicitly preserves a possible positive
volume. The exact error-distance machinery in G1 reports error=2^-1074 without
unnecessarily widening it. This does not assume that a rounded zero equals the
true volume, and it does not claim a relative certificate.

Isolated G2 resolves the remaining **4/18** requests. The combined explicit method
meets **18/18 original targets within the original budgets**, with all **167/167**
independent analytical enclosures contained. Source-bound results and original
failures are retained in `geometric_subnormal_validation.json`; raw logs and peak
RSS are under `.stericx/experimental_validation/geometric_subnormal_v1/`.
The four API tests include invalid-input handling, triple-overlap fallback,
permutation/duplicate invariance and the exact underflow threshold. Existing
certified API tests continue to exercise the original algorithm and its honest
budget failures.

These finite analytical successes do not establish that arbitrary tangent,
clipped, thin or multiply overlapping molecular geometries meet arbitrary error
targets. Regional stress tests and the complete regression remain separate gates.

## Expanded-budget and profiling replay

The frozen and new binaries were rerun on all eighteen requests at the original
budgets and at four times those budgets, retaining each tolerance. Baseline success
is 0/18 at both budgets; the new method is 18/18 at both. Every independent reference
interval is contained. All new accepted cases use analytical proofs with no adaptive
subdivision. Full per-case geometry, original/expanded work, status, bounds, runtime
and memory diagnostics are retained in `performance_v1/`; see [PERFORMANCE.md](PERFORMANCE.md).
The remaining four deliberately stringent regional requests at 1e-8 Å³ are different
experiments, remain honestly budget-limited, and are not counted as resolutions of
these eighteen historical requests.
