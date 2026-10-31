# Paton recommendations: implementation and evidence plan

Baseline source: `c9380d2c24ba5d48cd52ec13c139ee6e57f04991` (0.3.0).
This is an implementation ledger, not a scientific admission or release claim.
The initial working tree has user-owned untracked `docs/media/`; it is excluded
from this work. Existing audit evidence and failed experiments are immutable.

## Architecture and existing capability

`src/geometry/` implements f32-coordinate Sterimol, trivalent-donor buried
volume, pyramidalization, XYZ/V2000 parsing and authoritative SDF connectivity.
`src/descriptors.rs` is a binary-private CLI coordinator: it detects donors and
returns arithmetic ensemble means. It already parallelizes independent files
with ordered Rayon collection. `commands/parse.rs`, `buried_volume.rs`, `db.rs`
and `screen.rs` consume related ensemble reductions. Packed schemas and models
have separate compatibility and scientific contracts; no redesign is needed.
`kinetics/` implements checked thermodynamic/kinetic calculations. Python
scripts prepare RDKit/MMFF or CREST/xTB inputs and carry population provenance.
The scientific engine is Rust; Python will call it rather than duplicate it.

There are **no existing PyO3 bindings**. `pyproject.toml` describes the
`stericx-tools` research environment, not an importable `stericx` distribution.
No density-grid reader or density-volume kernel exists. Existing scan studies
vary integration density, not a public DBSTEP-style radius profile.

Library buried-volume defaults are radius 3.5 Å, nominal point volume 0.01 Å³,
center distance **2.1 Å**, radii scale 1.17, hydrogen occupancy excluded.
Several studies explicitly use 2.28 Å; an older remediation prose statement
calling that the default is superseded by the actual code and forensic report.
Keep defaults and legacy output arithmetic unchanged. Explicit donor-bound H
still participates in topology. Legacy missing-weight uniform aggregation is
documented behavior and must not become a claimed thermodynamic population.

## Feature inventory and intended commits

| Step | Baseline status | Missing work / validation |
|---|---|---|
| 0 Freeze | Extensive existing frozen evidence; new task freeze needed | Source/toolchain/dependency/executable identities, fresh tests and full replay; preserve failures |
| 1 Per-conformer descriptors | Partial: native kernels and volume conformer rows exist | Public ordered records with source text, stable IDs, metadata, independent descriptor errors; expose regional arrays |
| 2 Statistics | Partial: means and some volume extrema exist | Generic mean/min/max/population SD/median, explicit optional weighted mean; no invented weights |
| 3 Populations | Partial: validated preparation helpers, no public Rust energy API | Typed units, provenance and energy kind; stable normalization, degeneracy semantics, independent analytical reference |
| 4 Sensitivity | Partial: forensic ensemble extrema exist | Reusable report with range, SD, weighted mean if justified, all tied extreme identities |
| 5 Ensemble study | Partial: studies 001–011 and full Kraken extrema study exist | Frozen same-input ensemble comparison; representative/arithmetic/weighted/min/max; no refitted or invented historical populations |
| 6 DBSTEP environment | Absent | Official unmodified pinned source, isolated environment, complete package/settings/input/command receipt |
| 7 DBSTEP Sterimol | Absent | Shared atoms/axes/radii/L convention; independently specified criteria; all residuals and seven requested statistics |
| 8 DBSTEP volume | Absent | Compare matched definitions; characterize distinct lattices/normalizations and convergence |
| 9 Discrepancies | Extensive Morfeus/Kraken precedent; DBSTEP absent | Frozen witnesses, atom/frame/radius/integration checks, Morfeus triangulation, earliest divergence classification |
| 10 Density method | Absent | Official source plus primary literature, threshold/units/lattice/boundary semantics; implementation decision |
| 11 Density input | Absent; conditional | Strict compatible volumetric reader preserving origin/vectors/dimensions/units; explicit unsupported cases |
| 12 Density volume | Absent; conditional | Separate mode; explicit threshold, deterministic integration; independent DBSTEP and analytic validation before public exposure |
| 13 Distance profiles | Absent | Document exactly what DBSTEP varies; independent validation at each requested position |
| 14 Python core | Absent bindings; tools already exist | Thin PyO3 bridge, typed Python results and structured errors, optional NumPy/pandas/RDKit adapters and ordered batch API |
| 15 Python advanced | Absent | Expose verified conformer/stats/populations/profile APIs; density only if its gates pass |
| 16 Distribution | Absent | Maturin wheel/sdist, isolated install, package-name check, OS/Python CI matrix; no unrun platform claimed tested |
| 17 Tutorials | Partial CLI/research examples | Executable Python examples and complete notebook with units and interpretation |
| 18 Audit | Extensive existing independent suite | Fresh full replay, Morfeus/Kraken/DBSTEP, invariance/symmetry/precision/convergence/populations/new API tests |
| 19 Performance | Established frozen methodology | Alternating paired native baseline/final runs, single/batch/six-core/ensemble, RSS and allocations, Python overhead |
| 20 README | Broad model/research oriented | Descriptor-engine focus and scoped evidence; distinguish same-input agreement from historical reproduction |
| 21 Release | Version 0.3.0; no Python artifact | Candidate only after applicable gates, including actual cross-platform wheel tests; document blockers otherwise |

Use the user's requested commit messages for each completed step. A conditional
feature that fails scientific admission stays a documented research prototype.
Do not create a release commit implying a gate passed when it did not execute.

## New API contracts, fixed before implementation

* A conformer record has source order, a positional key independent of duplicate
  titles, original title/identifier, original geometry text, atom order and bonds,
  and raw source properties. Malformed records and descriptor failures remain
  records. No successful-subset reweighting is implicit.
* Each descriptor can fail independently. A nontrivalent molecule can still have
  explicitly specified bond-axis Sterimol. Undefined pyramidalization or a
  regional frame is an error, not zero. Index conventions are zero-based.
* Regional arrays carry their orientation and axis definition; regions from
  incompatible conformer frames are not silently averaged as laboratory sectors.
* Arithmetic statistics describe the supplied finite ensemble: mean=sum/n,
  population SD=sqrt(sum((x-mean)^2)/n), median=central value or mean of two
  central values. All ties are retained for extrema. Descriptor coverage and
  failed identities accompany results. Incomplete-ensemble statistics require an
  explicit policy; complete-ensemble statistics remain unavailable otherwise.
* Weights are absent unless supplied with declared interpretation or derived
  from complete comparable energies at positive finite T. Supplied weights are
  not automatically Boltzmann weights. Electronic, force-field and free-energy
  interpretations remain distinguishable. No energy is inferred from a title.
* Explicit degeneracy multiplies a state's partition contribution. Enumerated
  microstates default to one each; equal energy does not itself merge records
  or justify double-counting degeneracy. Use an energy shift/log normalization.
* New wrappers reuse the audited native kernels. Changes to old kernels require
  an independent scientific witness and explicit changed-output ledger.

## Evidence already available and its limits

The audit, remediation, Kraken forensic reports, studies, profiling and later
performance reports are indexed and hashed in the baseline manifest. Important
reading includes the domain geometry/model/thermodynamic reports, convention
ledger, full replay/reproduction instructions, claim inventories, and the
Morfeus benchmark's prospective equivalence policy. Historical documents retain
their at-the-time conclusions; later scoped corrections supersede them.

* Remediation: 128,021 observations in 13 lanes, 31,721 topology conformers,
  1,141,956 bins, 7,624 independent model comparisons. Existing failure and
  unsupported-case classifications are not erased by numerical replay.
* Kraken forensics: 1,541 ligands / 31,611 conformers / 56 historical reductions;
  exact historical-geometry identity established for **zero** ligands. Missing
  historical free energies prevent historical Boltzmann reproduction. Reference
  versions differ for 821/1036; ligand 369 has ensemble-version evidence; 390
  demonstrates discontinuous argmin amplification of a one-cell difference.
* Optional 3600-direction Sterimol already exists. Reuse it. Default B1 remains
  a 360-direction approximation; L/B5 are analytic supports on f32 geometry.
* Morfeus's later admitted benchmark covers nine finite-lattice volume outputs
  on 56 geometries with certified point decisions and independent IEEE models.
  It excludes continuous-volume accuracy, Sterimol and pyramidalization from
  its performance admission. Earlier all-descriptor speed claims were withdrawn.
* Accepted C1/C1b/C2 optimize ordered file work and repeated model quantiles.
  Rejected C3 slowed native workloads and must not be resurrected as a presumed
  optimization. Preserve affinity, thread limits, alternating order, complete
  output fingerprints and every sample; timing must not overlap heavy jobs.
* Historical Ni-hDA 2064 provenance conflict, failed outer model folds, poor
  interval coverage and unfavorable ranking remain unresolved. No new
  experimental predictive validation is available from this task alone.

## Study and acceptance policy

Before new reference outputs are inspected, freeze each comparison's input
selection, equations, conventions, metrics and acceptance rules. Residual
statistics are outcomes, never fitted tolerances. Separate finite-algorithm
agreement from continuum convergence and physical/predictive validity. Preserve
failed reference executions and unavailable comparisons with explicit status.

Use legally retained source geometries and source provenance. Ensemble behavior
can be tested against independent same-input descriptor reductions without
experimental labels. If a predictive ablation is made, use fixed group splits
and fit all preprocessing/selection in training folds; historical shared
scaffolds and unavailable weights remain visible. Do not claim ensembles always
improve accuracy, or compare one reduction to a differently defined target as
though they were the same quantity.

Density adoption requires a primary-source definition and verified reference
grid. Cube orbital fields, ambiguous units, incompatible vectors or incomplete
sampling domains must fail explicitly. Density-based occupation is a distinct
descriptor from hard-sphere occupation. Source thresholds may be offered as
documented conventions, never inferred universal physical constants.

## Execution status

Phase 0 has passing fresh engineering checks (322 Rust, 124 Python) and a fresh
source-bound native/observer build. `baseline/manifest.json` seals the existing
complete 128,021-record capture **after exact equality of all 44 Rust/Cargo
source files and every observation stream is verified**. This is reuse of the
same-source archived outputs, not a claim that a new reference campaign finished.
The simultaneously started fresh replay remains under
`.stericx/paton/baseline/observations/`; its receipt will be added separately on
completion. `baseline/observe.log` is live and is deliberately not part of the
initial immutable seal. Subsequent commits update this ledger with actual
evidence and unresolved gates.
