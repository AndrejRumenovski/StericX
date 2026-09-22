# Scientific remediation before performance optimization

This report describes corrections to the audited `b515c4f` implementation. The
original [scientific audit](../scientific_accuracy_audit/SCIENTIFIC_ACCURACY_AUDIT.md),
its negative findings, datasets, models and published predictions remain intact.
The current correction snapshot is **v4**, identified below. It is a baseline for
preserving explicitly tested numerical and behavioral contracts, not a claim
that every descriptor convention or predictive interpretation is verified.

**Admission status: corrected v4 is admitted for exact performance optimization
under the explicit contracts and limits in [ADMISSION.json](ADMISSION.json).**
No performance optimization has yet been accepted. The immutable optimization
freeze separately records the corrected source commit and benchmark evidence.

| Component | Corrected behavior / independent evidence | Important remaining limit |
| --- | --- | --- |
| Sterimol L/B5 | Axial/perpendicular support is evaluated without dropping small real axis tilts. Original worst Kraken witnesses are replayed. | f32 inputs/output; explicit axis and radii convention matter. |
| Sterimol B1 | Validated axis/frame handling; atom permutations preserve all captured values. | The retained 360-direction scan is an approximation, with measurable phase sensitivity. |
| Buried volume | Valid symmetric occupancy is accepted. Scalar total, percentage and near/far summaries average the same three donor-plane grids symmetrically. | Cartesian grid, finite precision, regional boundary normalization and radius-table differences remain explicit. |
| Connectivity / frames | SDF connectivity is authoritative; XYZ inference remains identified. Tied attachment axes require an explicit reference atom. Ambiguous planar coordination requires an explicit center. | Distance thresholds alone do not establish chemical bonding. |
| Pyramidalization | Invalid/degenerate inputs fail explicitly; independent reference comparison is retained. | Ordinary f32 output error and sensitivity close to degeneracy remain. |
| Conformers | Valid finite weights are normalized stably; nonfinite values fail; energy span is max minus min. | Missing weights retain the documented uniform policy; supplied weights are not automatically Boltzmann populations. |
| CREST | Requested temperature is passed and checked; reweighting requires complete individual-state energies. Old unverified population caches are not reused. | Mocked command/cache tests establish interface behavior, not conformer-search completeness or free-energy accuracy. |
| Eyring / selectivity | Difference-only input no longer invents an absolute rate. Stable analytical selectivity, SI constants, checked absolute rates and log rates are tested. | Equal prefactors/irreversible pathways are assumptions. Unsigned ee does not identify R/S. |
| Models / screening | Declared aggregation must match geometry input; malformed applicability metadata and nonfinite evaluation fail explicitly. Independent model checks: 7,624, zero failed comparisons. | Historical poor coverage, conditional feature-selection diagnostics and failed outer folds remain. |
| Database / search | Finite large means and standardized distances avoid f32 intermediate overflow; unavailable/unrepresentable results are explicit errors. | Descriptor-space proximity is not experimental similarity or calibrated reliability. |

## What was corrected

The [pre-edit ledger](remediation_ledger.json) separately inventories 43 negative
audit propositions. A numerical bug, a convention mismatch, an unsupported
statistical inference and an absent experiment require different responses.
Implementation changes are tied to failing inputs and independent equations;
unsupported broader claims are narrowed or withdrawn rather than made to pass
by changing a radius, density, target, dataset or reference result.

Geometry corrections include the demonstrated atom-order/axis failures,
near-parallel Sterimol alignment error, rejection of valid zero anisotropy,
nonfinite/degenerate public geometry results, and failure to honor known SDF
bonds. The scalar buried-volume orientation reduction is an intentional
**scientific convention correction**: it is the symmetric mean of three already
computed orientations. It is not a performance optimization and should not be
compared bitwise with the arbitrary historical first orientation. Quadrant and
octant extrema retain their definitions. Primary and secondary phosphines retain
donor-bound hydrogens in their frame even though hydrogen occupancy is excluded.

Planar sign selection is explicitly limited by representable input coordinates.
Clearance intervals compare the two sides of a fixed computed normal; these are
not uncertainty intervals over all possible normals or physical geometries.
Unresolved sides produce an error instead of an atom-index-dependent answer.
The full 31,721-conformer topology corpus contains no case entering that fallback.

Thermodynamic changes separate energy differences from absolute barriers, reject
invalid weights/descriptors, avoid premature f32 normalization, correct supplied
energy spans, and distinguish exact 100% ee censoring from finite conversion.
Future Ni-hDA response metadata is 353.15 K. Published target values and historical
298.15 K conformer populations are retained with their provenance; populations
were not relabeled or recalculated to improve a model comparison.

Model changes fix nonfinite evaluation, clipped/extremely central Student-t
quantiles, ordinary covariance-distance availability and a malformed-index API
panic. Schema 3 declares descriptor aggregation. Legacy unknown/precomputed
models refuse implicit geometry substitution. Supplied-weight screening uses the
same aggregation and axes as reaction parsing. Cross-coupling labels now use the
primary reference's inclusive activity boundary. Model, interval and applicability
wording describes what is computed without asserting chemical reliability.

New adversarial checks during remediation exposed finite mean overflow in
descriptor/DB aggregation and standardization overflow in search. Those failures
were preserved and corrected before the final snapshot. The preserved large-finite
descriptor/DB/search witnesses now produce finite results or explicit
unrepresentable-result errors.

## Frozen corrected implementation

The source snapshot and native/observer build receipt are
`.stericx/scientific_remediation/validated_build_v4/manifest.json`, SHA-256
`d79fa7ec5c9e6c1254aa3bd1de665fe4338dde2aabe2401b7bf79a36b5a7bdd3`.
Its native release executable SHA-256 is
`102b6883c639bfc9ef211a88abae02391419cbd277c498098fa2bdfc7008206f`.
The observer is an explicitly identified SUT measurement adapter, not an
independent scientific implementation.

Thirteen observation lanes contain **128,021 responses**, including all original
requests and two additional complete topology-aware corpus lanes. Requests,
outputs, errors, IEEE bits, private bins, executables, source, dependency locks,
helpers and reference environment are hashed. The original audit seal remains
`c31b655641a352402b14dc2d4535261a9de9c9004979c5e571d78e45cac82f03`.
Intermediate v1–v3 builds and failed harness attempts are retained separately.

## Independent results and scope

The unchanged audit reference programs have been rerun against both coordinate
inference and authoritative topology. The topology lane has 31,721 successful
conformers; the original inferred lane retains its eight extra-contact errors.
Known connectivity is additional input information, not a tuned bond cutoff.

The independent buried-volume campaign compares all **31,721 conformers,
95,163 orientations and 1,141,956 quadrant/octant bins**. Among 30,779 conformers
with matching radius tables, the largest observed total mean residual is
0.00777188794 Å³ (0.00432913408 percentage points). All 56 matching-table structures
with differing integer octant occupancy were minimized to 58 differing point
decisions. Frozen native grid/aligned-atom data reproduce native total/octant
bits exactly; the unchanged independent kernel reproduces its own counts.
The witnesses isolate arithmetic decisions at occupancy boundaries. These are
measured residuals, **not newly introduced acceptance tolerances**.

Across all elements, the largest total mean discrepancy is 2.84793579 Å³; the
larger differences include documented B, Br, Fe, Se and Sn radius/fallback
conventions. They remain in per-ID results, strata and largest-outlier tables.
Morfeus is an independent implementation, not physical ground truth. No radius
was changed to improve agreement.

All 120 captured atom permutations reproduce reported values/errors exactly.
All 618 rigid variants retain their success/error classification. All 515
nonplanar rigid variants have identical buried-volume outputs in this campaign.
The planar base/rotation/permutation family consistently reports ambiguity.
This finite campaign does not imply universal grid invariance. The largest
captured coordination-B1 rotation change is 0.0199442 Å, retained and disclosed.

The fixed reference union contains **2,517 IDs in each lane**. It includes every
original 590 diagnostic and 2,009 expanded-selection ID, as well as newly
selected discrepancies. Eleven original reference variants retain their error
and missing-result statuses. Reuse is allowed only for identical complete
request/SUT pairs; every missing union result is independently recomputed.
This avoids replacing the original failed-case selection with a newly easier
selection after code changes.

The [model recheck](models/REFERENCE_RECHECK.md) contains 7,175 comparisons plus
a 449-comparison completeness/training-geometry supplement. Of these, 4,400
historical numeric comparisons are exact. Both corrected fits and all three
documented JSON screens are independently checked. Seven original unsuccessful
outer folds remain unsuccessful. No predictive classification is upgraded by
these arithmetic checks.

The final thermodynamic replay contains 744 analytical kinetic scalars,
51 aggregation scalar/rejection comparisons, 14 actual CLI cases, six CREST,
11 MMFF and 24 ee cases. It preserves all published targets over 1,566 source
rows. Maximum relative error for normal representable rates is
5.29218×10⁻⁸; maximum ee absolute error is 3.69286×10⁻⁶ percentage points.
The checked/log APIs separately cover representational extremes.

The behavioral oracle captures 93 CLI cases at each of 1/2/4/6 threads.
All 279 comparisons with the one-thread capture match exactly under the
declared removal of named timing/resource fields and `portable.created_utc`. Required
successful artifacts must exist and be nonempty. Errors, exits, descriptor
values, rankings, exclusions and packed data are not normalized away.

The final engineering checks pass: 306 Rust tests, 113 Python tests and four
additional independent-reference integrity tests, plus Clippy with warnings
denied, rustdoc with warnings denied, Rust formatting and Ruff check/format.
These support software integrity; the separate reference results provide the
scientific evidence.

## Claims retained and withdrawn

The [convention ledger](CONVENTIONS_AND_LIMITS.md) and individual domain reports
distinguish tested implementation behavior from chemical interpretation. Safe
claims describe the chosen radii, attachment axis, finite grid, conformer
aggregation and conditional statistical model explicitly. Native precision
claims must use native/reference residuals, not a Python identity at float64
precision. The implemented subset of Kraken fields must be named individually.

StericX should not claim universal Morfeus/Kraken equivalence, exact continuous
B1, complete conformer populations, fully converged MMFF ensembles, full-pipeline
validation from fixed-feature diagnostics, calibrated predictive coverage from
nominal intervals, or experimental validation from screening rank. Formula
equality does not establish isomer identity. The source conflict for Ni-hDA
ligand 2064 and missing historical thermochemical provenance remain unresolved.

## Four separate conclusions

1. **Optimization correctness:** not yet established here. The subsequent
   performance task must compare the optimized build to this corrected baseline
   using exact complete scientific outputs, errors and records.
2. **Implementation correctness:** independent equations, reference programs,
   full-corpus residuals and retained adversarial witnesses support the scoped
   corrections, subject to the documented numerical and convention limits.
3. **Methodological validity:** radius/axis choices, finite integration,
   approximate populations and conditional statistical assumptions remain
   modeling decisions. Reproducing their arithmetic does not validate them for
   every chemical interpretation.
4. **Predictive/experimental validity:** no new prospective measurements were
   obtained. Unfavorable historical validation remains evidence.

## Reproduction and performance handoff

Follow [REPRODUCE.md](REPRODUCE.md), the domain reproduction programs and their
fixed input hashes. Heavy raw streams/binaries are retained under `.stericx`;
compact reports, witnesses and receipts are published with the repository.
Destinations must be new. Never overwrite a failed attempt or the original audit.

The [corrected benchmark inputs](profile_inputs_v2/README.md) contain ten native
workloads, prepared by the exact v4 executable. The full database includes 1,543
preselected ensembles/31,618 conformers with zero runtime skips. The 23 source
ensembles outside the declared supported input rule are retained in an explicit
preparation selection ledger; no benchmark failure drove selection. SDF bonds
and coordinate tokens are preserved. The largest complete SDF workload contains
86 conformers. These workloads use the coordination-axis contract explicitly;
old timings with different scientific inputs cannot be headline comparisons.

`scripts/scientifically_exact_optimization.py` admits a freeze only after a
reviewed scientific receipt, source/build binding and complete independent gates.
Optimization comparisons require the admitted baseline and exact complete
observations. No scientific tolerance is widened for speed. The native profiling
harness retains every warmup/measured run, wall/CPU/RSS/fault measurements,
thread count, affinity and scientific-output fingerprint. Candidate performance
claims must follow scientific gates and alternating native A/B measurements.
