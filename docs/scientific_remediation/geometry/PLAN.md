# Geometry remediation plan and evidence ledger

Status: plan recorded before production changes. The original scientific audit
is immutable. This work corrects demonstrated implementation failures and makes
undefined or ambiguous inputs explicit; it does not declare every historical
compatibility proposition true.

The complete negative G/K claim inventory is retained in `claim_triage.json`.
The source of each historical result remains the sealed audit, including its
qualifications and limits. The geometry report and methods, alignment dependency
source, raw requests/outputs, and independent references have been inspected.

## Evidence and decisions

| Claims | Finding and disposition | Independent witness / required outcome |
| --- | --- | --- |
| G22, K12 | Correct the rejection of a valid nonzero volume solely because quadrant asymmetry is zero. Geometric frame validation remains independent of occupancy anisotropy. | `phosphine_PH3`, the four `ph3_lens_*` cases, full-cover/small-sphere constructions, and `KRAKEN:1299:54318`. Valid frames must return the independently calculated zero asymmetry; PH3 default lattice volume is about 31.9372145 Å³, with a separately documented continuum/grid residual. |
| G39, K10; affected G01/G03/G07/G37 | Correct near-parallel/antiparallel Sterimol axis alignment. The frozen dependency deliberately approximates near-singular rotation arcs; StericX incorrectly treated that approximation as an exact supplied-axis projection. | Four frozen full/minimized `963:48394` / `1075:49864` cases and their exact 90° coordinate permutations; full-corpus extrema `16:31923` and `1407:54385`. Use independent direct dot/perpendicular projections on the identical input floats. |
| G35, G36 | Replace silent zero sentinels with checked `Result` APIs; evaluate finite, nonzero short axes and near-collinear but defined planes without an arbitrary absolute EPSILON rejection. Detect unrepresentable outputs. | Frozen short-axis neighbors of 0.000345267 Å, `nearly_collinear` (independent alpha about 29.999925°), invalid indices/coincident geometry, numerical NaN/Inf and finite 1e38-coordinate/radius cases. Valid cases must yield finite reference-consistent results; undefined/unrepresentable cases must return explicit errors. |
| G06 | A nearest-bond tie is an ambiguous automatic axis, not permission to select a chemically different axis by file order. Add an explicit reference-neighbor option and reject unresolved bond-axis ties deterministically. | All six frozen tied-axis permutations: implicit mode must consistently report ambiguity; mapped explicit axes must produce the same descriptors within the stated arithmetic/discretization limits. |
| K11; related G25/G26/G28 | Retain and use supplied SDF connectivity through additive topology-aware parsing and explicit-neighbor geometry APIs. Keep coordinate-only inference documented and unchanged; never retune the 1.3 cutoff to fit these examples. | All eight independently graph-checked conformers of IDs 1281/1907, plus the 20 historical primary/secondary phosphine cases and ordinary graph controls. Supplied graph has three neighbors even when a nearby nonbonded P/Si passes the distance heuristic. |
| G30; frame/order consequences | The planar fallback is a declared geometric model, not a physically established lone pair. Unresolved sign/frame ambiguity needs an explicit center or an error, not an index-dependent choice. Coordination-axis scalar calculations must not require an unrelated unique heavy-atom bond axis. | Planar permutations, clear and tied clearance cases, PH3 coordination mode, and ordered three-frame observations. Final reduction/ambiguity contract is being coordinated before implementation. |
| G33 | Correct the local API documentation: a tetrahedral triad has P=4/(3√3), while an orthogonal triad has P=1. | Existing analytic tetrahedral and orthogonal witnesses. No equation change is justified by this terminology error. |
| G12/G13, G14/G17/G18/G20, G23/G29/G37 | Retain distinct radius, grid, zero-plane, regional-normalization, center, and numerical-resolution conventions and their quantified limitations. Correct local overclaims; do not alter those conventions to fit Morfeus or published aggregates. | Original per-case table comparisons and convergence results remain evidence. Odd-grid regional nonadditivity and the B1 angular/grid error floor remain explicit numerical limitations. |
| K03/K04/K06/K15 | Historical narrative/precision or protocol-attribution errors, not instructions to change a valid radius table or imitate historical software. Correct current local API overclaims where applicable; global/study documentation belongs to the coordinating agent. | Primary Kraken DFT code/SI and complete native/independent corpus results. |
| K09/K13/K14, G38 | Unresolved historical provenance, missing source populations, or experimental validity cannot be repaired by tuning a descriptor kernel. | Preserve limitations and source-access/identity requirements; do not relabel them as corrected computational defects. |

## Proposed implementation contract

1. `SterimolCalculator::compute`, `compute_with_dummy`, and
   `PyramidalizationCalculator::compute` will return `Result<Params,
   DescriptorError>`. Every production caller must propagate the error. There
   will be no compatibility wrapper silently returning a measured zero.
2. Sterimol evaluates the mathematical axis from the supplied coordinates using
   wider intermediate arithmetic. Axial and radial supports use direct
   projections; the transverse basis for the existing one-degree B1 scan must
   preserve the supplied axis even near ±Z. Output storage stays `f32`, with
   explicit errors if a finite mathematical result cannot be represented.
   The 360-direction B1 discretization is retained and documented.
3. Pyramidalization uses stable intermediate vectors and the independent
   determinant/plane-angle equations on the existing rounded coordinates.
   Exactly undefined axes/planes are errors; a merely small nonzero length or
   normal is not automatically zero. Finite returned values are required.
4. Topology-aware parsing is additive: the existing `Molecule { atoms }` API
   remains available. A new frame wrapper retains V2000 bond records and
   exposes validated adjacency. SDF connectivity is authoritative when supplied;
   XYZ remains explicitly inferred geometry. New buried-volume/center methods
   accept the three supplied donor neighbors and validate them without a
   second, conflicting radius cutoff.
5. The CLI's explicit reference-neighbor index applies to a single ligand file.
   Bond-axis mode requires either that explicit reference or an unambiguous
   nearest bonded heavy atom. Coordination mode does not need a heavy bond axis.
   Public regional-reduction and planar-sign ambiguity changes will be recorded
   precisely after coordination, before implementation.
6. Remove unsupported local wording about an electronic Kraken DFT center,
   physically conservative unknown-element radii, universal distance bonding,
   and unphysical symmetric occupancy. Preserve the actual declared numerical
   conventions rather than silently replacing them.

The thermodynamics agent owns only `BuriedVolumeCalculator::aggregate` and a
separate aggregation-test file for C16/C17. Its finite-input and stable-weight
normalization correction is independently evidenced there; the geometry owner
retains the rest of the shared buried-volume source.

## Before/after evidence and acceptance

Freeze failing outputs before implementation using the original audited/current
source and raw witnesses. Add focused tests that fail on that source for each
justified correction, then retain their passing results after correction.
Compare projection/angle values with independently evaluated equations, not with
newly generated expected values from the corrected implementation. Preserve
every relevant negative case and mark intentional behavior/API changes.

Exercise atom permutations, exact coordinate rotations, ordinary translations,
explicit topology versus coordinate inference, finite/subnormal/extreme values,
and unchanged ordinary graph cases. Region-bin reference equations remain
independent; any changed reduction rule must have its own mathematical and
per-case evidence. Do not promote a newly observed maximum residual to a global
tolerance. Full independent and engineering gates follow focused validation.

The old audit observer cannot compile the changed public return types. A new
remediation adapter will be placed outside the historical audit and will record
its exact adapter changes, current source, inputs, binaries, errors and outputs.
Historical requests, outputs, scripts and classifications will not be edited.

## Sources and access limits

- [Sealed geometry methods](../../scientific_accuracy_audit/geometry/METHODS.md)
  specify the independent support, volume and pyramidalization equations.
- [Sealed geometry report](../../scientific_accuracy_audit/geometry/REPORT.md)
  retains the constructive and numerical witnesses.
- [Sealed Kraken report](../../scientific_accuracy_audit/kraken/REPORT.md)
  records primary DFT code/SI provenance and all eight topology discrepancies.
- [Official Morfeus Sterimol documentation](https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html)
  describes explicit dummy/base axes and the distinct raw versus corrected L.
- [Official Morfeus volume documentation](https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html)
  documents its descriptor conventions; those choices are not universal truth.
- [Official pyramidalization API](https://digital-chemistry-laboratory.github.io/morfeus/api/morfeus.pyramidalization.html)
  and frozen versioned source supplement the independent determinant/angle
  derivation. The original historical papers' access limitations remain those
  recorded in the audit; this remediation does not claim new full-text access.
