# Final scientific replay

Status: the final legacy capture, engineering checks, new numerical methods,
advanced references and matched finite-grid competitive gate have passed in
the scopes below. The full independent Kraken/volume reference replay and its
review are complete with zero changed historical scientific values/statuses.
Controlled performance acceptance is separate and pending.

## Frozen source and exact legacy behavior

Original commit: `c03630cb3595eb85b5d05bdbecd3f10305b818bf`.
Final engine implementation: `f9c952e2` (later campaign documentation does not
change production code). Source-bound manifests are retained at
`.stericx/ultimate_engine/{baseline,final}_scientific_build/manifest.json`.

All **128,021 responses across 13 lanes match exactly** between frozen baseline
and final engine, including serialized floating-point bits, ordering and errors.
This includes both inferred and supplied-topology Kraken lanes, all regional
bins, geometry degeneracies/rotations, aggregation, kinetics and model-domain
witnesses. Exact replay preserves existing behavior; independent mathematical
checks are needed to establish its scientific meaning. S1's new small-target
R² witness is separately corrected and is not present in the historical stream.

The 93-case complete CLI capture also matches exactly between 1 and 2/4/6
threads, comparing return codes, ordered stdout/stderr and produced artifacts,
excluding only declared process metrics and model creation timestamps. The
comparison does not silently drop failed inputs or change their ordering.

## Engineering and platform checks

| Check | Final result | Scope or limitation |
|---|---|---|
| Rust all targets/all features | 364 passed | x86-64, Rust 1.97 |
| Rust minimum version | 364 passed | actual Rust 1.85.0; baseline build failure retained |
| AArch64/musl tests | 362 passed | QEMU execution of ARM binaries; two x86-only tests are gated |
| Formatting and strict Clippy | Passed | all targets/features, warnings denied |
| Rust documentation | Passed | warnings denied, isolated frozen scientific source |
| Python binding API | 71 passed | freshly built local wheel, source-bound engine |
| Python research tests | 156 passed, no skips | combined process-local environment, recorded package identities |
| Pinned DBSTEP integration regressions | 28 passed | original pinned reference environment and final observer |
| Upstream DBSTEP tests | 413 passed, 3 skipped | upstream skips retained and not counted as passes |
| Archived Morfeus evidence integrity | Passed | verifies historical evidence, not a fresh scientific comparison |
| Memcheck | Zero errors/lost bytes | final ensemble and separately tested new numerical paths; 544 runtime bytes reachable |

The ARM/MSRV executable code matches the final source. One later documentation
comment clarifies a numerical limitation. The pre-run hashes and post-run
reconstruction are explicitly recorded in [arm_source_provenance.json](arm_source_provenance.json);
the reconstruction is not mislabelled an original archive. Real ARM hardware,
macOS and Windows execution were unavailable in this workspace. No native
performance or universal cross-platform bitwise claim follows from emulation.

The initial root rustdoc attempt was interrupted after Cargo traversed
historical `.stericx` build directories on NTFS for over nine minutes. Its
nonzero receipt/open-file evidence is retained. The unchanged check then passed
against the isolated frozen source. A monolithic run in the pinned reference
environment failed because pandas was absent; neither environment was modified.
The combined process-local replay appends reference packages after project
packages and records their actual identities. Original pinned scientific
comparisons remain separate. These failures are not erased or reported as
scientific failures of the Rust engine.

A separate final formatting check found compact formatting in the campaign's
observer crate (outside the production workspace checks). It was formatted and
the check now passes; [observer_formatting.json](observer_formatting.json)
preserves before/after hashes and the original failure. Timed and validated
executables remain bound to their unchanged archived observer source.

## Independent comparisons

The full independent replay recomputes 31,721 volume references with no reference
errors, 95,163 frames and 1,141,956 regional bins. The fixed 2,517-variant union
is replayed in both inferred- and supplied-topology lanes, preserving the same
membership and coverage. The inferred lane retains eight missing matched-center
cases; the supplied-topology lane has all 2,517. These are retained outcomes,
not silently excluded successes.

The final review compares scientific values and statuses exactly against sealed
historical reference artifacts and finds **zero changes**. It also checks source,
input, binary, command and raw-output hashes. The controller adaptation only
handles current receipt paths, the new exact-replay manifest and the preserved
rustdoc interruption/resume; the reference equations and comparisons are
unchanged. Its diff is [review_final_science.diff](review_final_science.diff).
This proves preservation, not universal agreement with differently configured
Kraken/Morfeus calculations. Historical global-pass fields remain false.

| Family | Evidence |
|---|---|
| Legacy DBSTEP Sterimol | 61 cases succeed with no unresolved cases; backward-axis L convention differences are retained |
| Common-node volume | 69 cases, 207 frames; zero disputed atomic points, unexplained sphere points or arithmetic failures |
| Matched competitive volume | 56 conformers, 168 frames, 2,588,544 certified point decisions; StericX, Morfeus and matched DBSTEP workflow pass unchanged exact-decimal/IEEE rules |
| Competitive ordering | All 56 outputs exact per tool at 1/2/4/6/12 workers |
| Density | All 12 analytic/DBSTEP cases pass; Gaussian convergence residuals retained |
| Radius profiles | Four requests, 20 positions, 16 reference comparisons and four expected per-position errors; no preservation failures |
| Explicit populations | 42 valid and 11 invalid cases; maximum absolute weight error 6.05e-17 against high-precision references |
| Conformer study | 11 ligands, 56 conformers; no fresh reference, conformer or ensemble investigations |
| Models | 7,175 comparisons plus 449 supplement checks; inspect their individual pass fields, not command exit alone |
| Thermodynamics | 744 kinetic scalar and 51 aggregation/rejection comparisons, 14 CLI cases; scoped numerical/domain replay passes |
| Continuous Sterimol | 39 independent 80-digit cases and matched Morfeus/DBSTEP comparisons; exact final-output replay |
| Continuous volume | 244 analytic/adversarial observations, 167 reference bounds, 22 exact permutation checks; 36 final molecular requests replay exactly |

The original hard-sphere DBSTEP gate still fails on the documented ligand-401
R=2.5 case: StericX retains 5,565 sphere nodes, while DBSTEP's extra binary64
clip retains 5,533. Three frame residuals remain approximately -0.024 Å³. The
fresh common-node check accounts for that difference without changing the
production grid, radii or tolerance. The radius-profile historical strict gate
therefore remains false; preservation and scientific explanation pass.

Density Gaussian absolute errors at spacings 0.25/0.125/0.0625 Bohr are
0.03955/0.01206/0.003122 Å³ against the analytic 0.6207147 Å³ target. Their
independently derived boundary-layer bounds pass. These cases do not establish
convergence for arbitrary density fields or validate a physical threshold.

New Reference-volume bounds contain their independent analytical enclosures.
Twenty-three analytic/stress work-budget or underflow targets remain unmet,
and all 12 strict molecular ±0.05 Å³ targets remain unmet at the tested budget.
Every such result reports `target_met=false`. Independent molecular octree
enclosures overlap, but this is not an exact molecular-volume oracle. N1's
high-precision bracket is not a certified floating-point error bound.

## Evidence locations and reproduction

Raw evidence is retained locally under `.stericx/ultimate_engine/`, including
all failed attempts. Compact receipts and figures are tracked in this document
directory. Controllers refuse existing output destinations so later runs do
not overwrite historical evidence.

```bash
CARGO_BUILD_JOBS=2 .venv/bin/python scripts/replay_scientific_remediation.py build \
  --oracle .stericx/profiling/scientifically_validated_optimization/scientific_oracle_v3 \
  --adapter docs/scientific_remediation/observer --output NEW_BUILD
.venv/bin/python scripts/replay_scientific_remediation.py observe \
  --build NEW_BUILD --output NEW_OBSERVATIONS
.venv/bin/python docs/ultimate_engine/compare_replay.py \
  .stericx/ultimate_engine/baseline_observations NEW_OBSERVATIONS NEW_COMPARISON
.venv/bin/python scripts/replay_dbstep_resolution_science.py \
  --build NEW_BUILD --output NEW_SCIENCE --logs NEW_LOGS --phase references
```

The science phase expects observations under `NEW_SCIENCE/observations`.
Use the retained controllers/receipts for exact paths, environment, reference
versions and all remaining phases. No tolerance is enlarged to turn a failed
scientific comparison into a pass. Data-source, model, finite-grid and
floating-point limitations are detailed in [SCIENTIFIC_AUDIT.md](SCIENTIFIC_AUDIT.md).
