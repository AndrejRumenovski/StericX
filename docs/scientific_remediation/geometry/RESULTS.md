# Geometry corrections and remaining limits

The original audit remains immutable. `PLAN.md` and `claim_triage.json` record
the pre-edit plan; this document records the implemented contract and focused
evidence. A full corrected-corpus replay is a separate acceptance step. These
results do not turn all historical compatibility or experimental claims into
verified scientific statements.

## Corrected behavior

`SterimolCalculator::compute`, `compute_with_dummy`, and
`PyramidalizationCalculator::compute` now return `Result<Params,
geometry::DescriptorError>`. Invalid indices, undefined axes/planes and invalid
numeric inputs produce errors instead of indistinguishable measured zeroes.
Production callers propagate the errors. Migration is explicit: callers must
handle the `Result`; there is no silent-zero compatibility wrapper.

Sterimol L and B5 use direct axial and perpendicular support. A stable
shortest-arc transverse basis preserves the supplied axis near either pole;
the previous quaternion approximation could snap a nonzero tilt to a parallel
or antiparallel axis. The existing 360-direction, one-degree B1 scan remains
f32. Its angular phase and discretization limit remain relevant, including
after a coordinate rotation. Inputs and outputs remain f32. Wider intermediates
are limited to vector/plane constructions implicated by the finite-overflow,
tiny-axis and nearly-collinear witnesses, plus finite ensemble summation.
The BV grid, occupancy predicate, radius table and regional normalization are
unchanged. BV frame construction uses wider intermediates only on the formerly
rejected small/nonfinite intermediate path.

A finite nonzero buried volume with zero quadrant asymmetry is valid. Frame
validation no longer treats that measured symmetry as a degeneracy. Total,
percentage, near and far volumes now use the arithmetic mean over the **same
three donor-plane grids** already evaluated for regional extrema. Previously
those four values came from the first selected plane. Each scalar's three f32
values are sorted, divided by three and summed in that reproducible order;
percentage is derived from the averaged total and sphere volume. Quadrant and
octant extrema and maximum adjacent-quadrant difference retain their reductions
over all three planes. This is a declared estimator correction for atom-order
dependence. It is not an assertion that any finite grid is continuously rotation
invariant. The observer retains each original plane's individual values.

The geometric center still uses the opposing sum of **unit** donor-bond
directions; it is not relabeled as either an electronic orbital center or the
original Kraken DFT sum of raw bond vectors. Unit directions are ordered before
their sum to avoid index-order accumulation differences. A planar clearance
tie has no index-independent sign and returns an explicit ambiguity error.
The first correction used exact f32 clearance equality; full-corpus review then
showed that 95 of 124 rounded planar variants spuriously selected a sign. The
final planar fallback compares outward-rounded coordinate-rounding clearance
intervals for its fixed represented trial normal. Each input component has half
its adjacent f32 spacing as a rounding halfwidth; overlapping nearest-clearance
intervals require an explicit center. All 124 original planar variants and an
unambiguous off-plane control pass the corrected regression. This propagates
coordinate rounding only; uncertainty in the inferred normal itself is **not**
enclosed. The existing nonplanar branch and planarity threshold are unchanged.
`planar_precision/INTERVAL_CONTRACT.md` gives the equations, and
`outward_interval_evidence.json` retains every conditional interval endpoint.
`compute_with_center` or `compute_with_center_and_neighbors` supplies the
resolution when a center is scientifically specified.

`MolecularFrame`, `parse_coordinate_file_with_topology` and
`parse_sdf_with_topology` retain V2000 connectivity. Bond indices, self bonds,
duplicates and truncated bond records are checked. A supplied zero-bond graph
is authoritative, not a request for distance inference. The legacy coordinate
parsers still return `Molecule` objects, now validating supplied records first.
New `compute_with_neighbors`, `compute_with_center_and_neighbors` and
`coordination_center_with_neighbors` APIs bypass the geometric bond cutoff for
validated supplied topology. XYZ retains the unchanged 1.3-times-covalent-radii
heuristic and its limitations. No bond threshold or radius was tuned.

The descriptor CLI uses SDF connectivity and accepts `--reference-index` for
an explicit bond axis. An unresolved nearest-heavy-bond tie is an error. The
index options apply to a single file and are rejected before opening inputs in
batch mode. Coordination-axis mode needs no unique heavy-bond axis and therefore
can compute PH3 and symmetric donors. Arithmetic conformer means validate finite
fields and accumulate/divide in f64 before a checked f32 result: repeating a
finite conformer can no longer overflow an intermediate sum into JSON `null`.
The same demonstrated overflow in the DB export mean was corrected separately
without changing its grouping convention; `db_mean` preserves the four-file
before/after CLI evidence and four regression tests.

## Retained evidence

| Witness | Historical failure | Corrected focused result |
| --- | --- | --- |
| PH3 and Kraken `1299:54318` | Valid symmetric occupied sphere rejected | PH3 V=31.9372138977 Å³ and asymmetry=0; sealed independent lattice V=31.9372145173 Å³ |
| Controlled `963:48394` | L=6.8599996567 Å from snapped axis | L=6.8582773209 Å; independent 6.8582773365 Å |
| Controlled `1075:49864` | B5=7.1169853210 Å from snapped axis | B5=7.1142158508 Å; independent 7.1142157096 Å |
| Nearly collinear triad | P=alpha=0 sentinel | alpha=29.9999256134°; independent 29.9999249736° |
| Finite 1e38-coordinate atom | B5=Infinity | B5=1.6439899498e37 Å, finite |
| Four identical finite extreme SDF frames | Successful JSON L=`null`; one frame finite | Repeated-frame full and screening means equal the single-frame value |
| Eight graph-checked Kraken 1281/1907 conformers | Inferred fourth P/Si neighbor rejects BV/coordination Sterimol | All eight supplied-topology BV calculations succeed; original inferred requests still expose the documented ambiguity |
| Six tied-axis atom permutations | Automatic physical axis changes with first atom index | All implicit bond-axis cases error; mapped explicit axis gives identical Sterimol and averaged BV values |

`focused_v1` contains 44 original requests, preserved historical outputs,
development observations and independently derived residuals. Its source hashes
identify that development measurement; the later final frozen replay supersedes
it for whole-program acceptance. The original reference equations are loaded
only after checking the pre-remediation audit-manifest SHA and the reference
entry's hash/size. `analyze_witnesses.py` records the complete four-ordinary-case
change distribution, not just selected corrected outliers. Across those four
phosphines, bond L/B changes are at most 7.16e-7 Å, P changes at most 1.20e-7 and
alpha changes at most 4.77e-6°. PH3 alpha changes 5.53e-5° toward its independent
value. The antiparallel controlled B1 changes 0.02185 Å: correcting the supplied
axis also changes transverse scan phase; this is reported separately from the
analytic L/B5 correction and does not create a new B1 tolerance.

`topology_witnesses.jsonl` preserves all nine original operation-failure requests;
its manifest proves their neighbor lists agree with both independent SDF and
original primary connectivity. `topology_v1/independent_v2` independently
recomputes all three planes and all regional values for the nine supplied-graph
cases: 9 compared, 0 reference errors. `recompute_volume_means.py` supports the
full request stream with a required expected row count, ordered ID validation,
bounded workers, original-reference seal verification, package/environment
versions and new output directories. Its grids and `volume_on_grid` function
are the unchanged independent audit implementation, never the SUT kernel.
Captured centers are controlled to isolate volume/frame calculations; this does
not independently validate the geometric center convention.

`recombine_reference_frames.py` separately recombines the preserved 590-case
Kraken reference variants. Missing/error variants stay explicit. It does not
replace their old center, density, radius or regional-boundary conventions.
The final full-corpus report must distinguish captured SUT rows from independently
recomputed reference rows and retain all exceptions.

Focused verification covers 41 geometry-related library tests, 12 descriptor
tests, eight new independent-witness integration tests and three batch-order
integration tests. Clippy with all targets/features and warnings denied passes.
Four helper-contract tests additionally reject changed reference code before
execution, changed audit seals, original-audit output paths, incorrect row counts
and source-ID mismatches beyond a limited computation prefix. Ruff checks pass.
The coordinating replay owns the final full-corpus, engineering and admission
receipts.

## Limits that remain limits

The old odd-grid zero-plane/regional-normalization discrepancy remains. For
`ph3_lens_0.1`, independently recombined total volume agrees to about 1.25e-6 Å³,
while near volume differs by 4.853809 Å³ because of that declared regional
convention. It is not a new allowable global residual. Float32/float64 lattice
boundary decisions, sampled B1, radius-table differences and the distinction
between unit-vector and raw-vector/electronic centers remain separately
classified. Historical ensemble membership, populations and software precision
remain incompletely established. No kernel fit can resolve that missing
provenance or prove experimental effect-size validity.

## Final frozen v4 review

The final source/binary receipt is
`.stericx/scientific_remediation/validated_build_v4/manifest.json` (SHA-256
`d79fa7ec5c9e6c1254aa3bd1de665fe4338dde2aabe2401b7bf79a36b5a7bdd3`).
The observation receipt is
`.stericx/scientific_remediation/observations_v4/manifest.json` (SHA-256
`5120602e43e9d914587f94e0d2ea3caa396db2e189586f269ae36127e180ef30`).
`final_review_v4_2/review.json` compares all ten geometry-related lanes to the
immutable audited baseline, retaining every success/error transition and observed
numeric change. The earlier `final_review_v4/failed_attempt.json` preserves a
request-schema check failure: explicit false dump defaults in the new topology
lane were subsequently recognized without relaxing chemical-input identity.

All 120 original atom permutations reproduce their parent scalar bits and error
status exactly. All 618 rigid-transform requests preserve acceptance status;
the five nonplanar groups (515 variants) reproduce all public BV scalar bits.
All 124 planar base/permutation/rigid requests consistently require an explicit
center. This is the original fixed invariance campaign, not a claim that all
31,721 conformers have been exhaustively permuted. Sampled B1 retains its known
rotation-phase dependence: observed maximum changes are 0.01835823059 Å for the
bond axis and 0.01994419098 Å for the coordination axis. L/B5 changes under
rounded rigid transforms are at most 9.5367432e-6 Å; near-planar pyramidalization
alpha changes reach 0.0002593994 degrees. Those are measured distributions, not
new tolerance limits.

Every one of the 31,721 explicit-topology geometries computes center, both
Sterimol forms, pyramidalization and BV successfully. The historical inference
lane still reports its eight extra-contact topology errors; this is deliberately
not substituted with the supplied graph. The formerly rejected symmetric
`1299:54318` succeeds in both lanes. No successful nonfinite scalar was observed,
and every changed success/error status has a recorded scientific explanation.
The full matched-center analytic support reference has maximum absolute errors
1.3137623913195284e-6 Å for L and 6.983686944295187e-7 Å for B5 across all 31,721
topology cases. It uses observed f32 inputs promoted to f64 and independent
closed-form equations; this isolates support arithmetic from center conventions.

The fresh independent default-grid BV campaign covers **all 31,721 geometries,
95,163 identified planes and 1,141,956 quadrant/octant bins**. Every plane identity
and grid count matches. `all_bins_reference_v4` retains every full-precision
public/frame/bin residual, complete metrics and top-20 lists. The independent
radius table is retained, so convention differences remain separate:

| Radius stratum | Conformers | Maximum absolute mean-volume residual (Å³) |
| --- | ---: | ---: |
| Same table after f32 representation | 30,779 | 0.007771887941594002 |
| B | 241 | 1.8338673756052728 |
| Br | 205 | 0.559481377378475 |
| Fe | 244 | 0.18261729082057343 |
| Se | 19 | 2.785762973034174 |
| Sn | 233 | 2.847935794306096 |

For the same-table stratum, `boundary_witnesses_v4` identifies all 56 rows with
different reconstructed integer octant occupancies. It reproduces the captured
native total and all eight octant IEEE bits for the 56 affected frames, and the
unchanged independent cKDTree occupied counts. There are 58 changed point
decisions. Each minimized point/atom witness retains the actual f32 grid/atom,
independent f64 grid/frame/radius and signed squared-distance margins; maximum
absolute native and independent margins are respectively 3.337860107421875e-6
and 2.5856973193683075e-6 Å². Thus the one-point bin/two-point frame discrepancies
have direct arithmetic-boundary witnesses. No empirical error ceiling is adopted.
This diagnosis covers the changed bin-count cases; it is not an assertion that
all unchanged-count frames have identical internal point-membership masks.

`final_geometry_frame_means_v4.json` independently recombines the three planes
for all 685 available successful original geometry references. Another 144
reference records lack three successful planes and remain explicit. Its largest
residuals retain the original radius and odd-grid convention distinctions. The
unchanged audit comparison still reports its historical first-frame estimator;
those raw outputs are preserved, while the recombination measures the corrected
three-plane mean.

The adaptive Kraken outlier selectors are not treated as a fixed historical
replay. `.stericx/scientific_remediation/fixed_reference_union_v4` explicitly
unites the original 590 and 2,009 case sets with both new selections: **2,517
unique geometries in each of the inferred and topology lanes**. All original
chemical requests and eleven unchanged reference variants are retained, including
BV density 0.01 and 0.001, primary/captured centers and radius variants. A fresh
result is reused only when its sealed complete request/SUT pair is exactly equal;
87 missing inferred and 100 missing topology records were calculated anew.
Original source equations and primary connectivity code are audit-hash verified.
Inferred-center variants are unavailable for the eight documented topology
errors; every topology variant executes successfully.

`fixed_variant_comparison_v4` retains the separately recombined means, all nine
BV field residuals, metrics and top-20 lists. The unchanged historical matched
helper uses distance-inferred neighbor lists, so its four-plane K11 records are
explicitly excluded from a *three-topological-plane* mean, not silently changed:
2,509 matched-plane means are comparable and all 2,517 primary means are retained.
At matched density 0.01, the largest mean-volume residual is
0.0038884576846385244 Å³. At density 0.001, comparison to the native default grid
reaches 0.5791799613123487 Å³; primary-center variants reach 4.666251879688389 Å³
(density 0.01) and 4.852038070301688 Å³ (density 0.001). These remain protocol
comparisons, not implementation errors or newly permitted residual thresholds.

The geometry review supports freezing this corrected implementation with the
explicit contracts above: the demonstrated symmetry, axis-alignment, invalid
sentinel, finite-mean, order-dependence and supplied-topology defects are corrected;
no new unexplained regression was identified in this scope. It does **not** certify
unqualified historical compatibility or global scientific validity. Sampled B1,
odd-grid regional normalization, radius choices, geometric-center interpretation,
XYZ connectivity and historical ensemble provenance remain classified limits.
The coordinator owns the final cross-domain admission decision. DB/search final
v4 checks and their exact f32 serialization comparison are separately documented
in `db_mean/RESULTS.md`. All newly added geometry evidence scripts pass Ruff and
the four reference-integrity contract tests still pass.
