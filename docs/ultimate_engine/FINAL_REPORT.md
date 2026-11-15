# StericX scientific accuracy and performance campaign

Status: implementation, profiling, numerical-method validation and full scientific
replay are complete in the scopes below. **Final controlled A/B and competitive
timings remain pending. The campaign is not yet complete.** No final throughput
or universal speed claim is inferred from exploratory measurements.

The post-validation [host check](host_after_science.json) measured 0.351 busy
CPU-core equivalents over 10.012 seconds, above the prospective 0.15 limit;
I/O wait was 0.036 cores. Warp and GNOME Shell were the largest background
users. No final timing was started. Following the user's Phase 10 instruction,
a 20–30 minute quiet-window confirmation has been requested. No background
application was stopped and the admission limit was not relaxed.

After the user reactivated the goal, a [fresh host sample](host_after_resume_1.json)
still failed admission: 0.380 busy core equivalents over 10.013 seconds, with
Warp and GNOME Shell again dominant. The user was asked to minimize Warp and
confirm an idle window. The user then confirmed readiness. The first measurement
controller's [admission check](final_window_1_admission.json) still failed at
0.286 busy core equivalents: GNOME Shell fell to 0.023, but Warp remained at
0.171. The controller stopped before collecting timing samples.

An explicitly scoped temporary suspension of Warp's UI process was proposed,
with automatic resumption and a separate 30-minute watchdog. Approval for this
additional application-control action is pending. This Codex session is not
descended from Warp; Warp's shell/Codex terminal processes would not be signalled.
The [guard checks](quiet_guard_checks.json) cover normal return, command failure
and a killed controller using only newly created synthetic processes. No user
process has been suspended. Final performance remains unmeasured.

## Scientific outcome

The original engine is frozen at `c03630cb3595eb85b5d05bdbecd3f10305b818bf`.
The final production source is `f9c952e2`. All 128,021 historical responses across
13 lanes match exactly, including errors, ordering and floating-point results.
The deliberate R² correction has its own independently derived regression
witness outside that historical response stream. Default steric definitions,
constants, radii, grids and density thresholds remain unchanged.

Two confirmed defects were corrected: R² incorrectly disappeared for small
nonconstant targets because absolute variance thresholds depended on target
units; and existing syntax prevented compilation with the declared Rust 1.85
minimum. Actual minimum-toolchain tests now pass. No chemical constant was
adjusted to fit an external result.

New explicit Rust APIs provide continuous-angle f64 Sterimol, midpoint-slice
continuous sphere-union volume, and conservative continuous-volume reference
enclosures. These are distinct quantities/methods with explicit limitations,
not silent replacements for finite-grid descriptors. The continuous support
method is checked against 39 independent high-precision cases. Volume methods
have 244 analytical/adversarial observations, 167 reference enclosures and
independent molecular interval comparisons. Every unmet reference tolerance
remains reported as `target_met=false`; all 12 strict molecular targets remain
unmet at the tested box budget.

The complete independent replay includes 31,721 volume references, 95,163
frames and 1,141,956 regional bins, plus the fixed 2,517-variant Kraken union
in both inferred- and supplied-topology lanes. Historical negative results
remain negative. The documented DBSTEP extra sphere clip is separately resolved
on 69 common-node cases/207 frames without altering a grid or tolerance.

See [SCIENTIFIC_AUDIT.md](SCIENTIFIC_AUDIT.md),
[NUMERICAL_ACCURACY.md](NUMERICAL_ACCURACY.md) and
[VALIDATION_REPORT.md](VALIDATION_REPORT.md) for definitions, proofs, reference
scope, failures and numerical/physical limitations.

## Profiling and implementation

Profiling preceded production optimization. Occupancy accounted for 81.85% and
grid construction for 9.08% of the single-thread instrumented 10,000-file batch.
The three retained performance changes address measured costs:

| Change | Purpose | Scientific preservation |
|---|---|---|
| U1 borrowed density tokens | Avoid allocating ordinary numeric tokens | Identical parsing, Fortran exponents and errors |
| U2 bounded per-thread grid reuse | Avoid rebuilding identical integration lattices | Exact radius/density bit keys, unchanged node order |
| U3 AVX2 row masks | Test eight atomic occupancies at once | Same multiply/add grouping and ordered predicate; guarded scalar fallback |

In the matched ensemble Callgrind workload, modeled instructions fall from
297,765,939 to 130,723,530 (56.1%). DHAT requested bytes fall from 17,142,207 to
985,652 (94.25%), while peak live bytes **rise** from 460,176 to 605,663 because
the grid remains cached. Density-token allocations fall from 140,636 to 12 per
Ne cube. These are instruction/allocation observations, not hardware-cycle or
elapsed-time speedups. Every unfavorable memory result remains in the report.

[PROFILING_REPORT.md](PROFILING_REPORT.md) includes before/after instruction
flamegraphs, simulated cache/branch events, assembly, dynamic SIMD dispatch,
thread scaling, Python layers and the five ranked measured bottlenecks.
Occupancy remains the largest computational cost. Retained grid capacity,
changing-radius scans and Python object transport remain opportunities.

Hardware performance counters were denied at `perf_event_paranoid=4`.
Valgrind's modeled events are explicitly labelled; cycles, IPC, actual cache
misses and bandwidth are unavailable. LLVM remark attempts produced no useful
loop reports, so assembly and dynamic execution support the SIMD conclusions.

## Performance measurement status

The frozen original exploratory 10,000-file baseline measured 1,977.2 files/s
at one thread and 10,264.5 at six; the older 10,774.5 value is historical context.
Exploratory combined kernel medians improved about 2.93× for volume and 2.53×
for a six-conformer ensemble. These desktop measurements do not establish final
files/s, accepted percentage improvement or achievement of the 15,000/20,000
targets. Final alternating baseline/candidate trials are still required.

The final measurement suite is prepared and functionally preflighted: seven
CLI workloads at 1/2/4/6/12 threads, isolated kernels, native/PyO3/public Python
layers, conversion costs, per-file versus batch calls, accuracy-versus-runtime
curves and matched external workflows. It retains every trial and fingerprint,
warmups, medians, dispersion, CPU time and RSS. Host activity is checked before
timing; a busy host requires the user-requested quiet window.

The competitive volume gate passes 2,588,544 independently certified occupancy
decisions over 56 conformers/168 frames. DBSTEP is measured with an explicitly
matched common-grid kernel workflow, not its different default sphere-clipping
workflow. No blanket combined-descriptor or whole-library speedup is admitted.
See [COMPETITIVE_BENCHMARKS.md](COMPETITIVE_BENCHMARKS.md) for exact scope,
separate process-tree memory observations and the pending timing distributions.

## Engineering validation

Rust 1.97 and 1.85 each pass 364 tests; AArch64/musl under QEMU passes 362
applicable tests. Strict formatting, Clippy and rustdoc pass. The fresh binding
passes 71 tests; the combined research environment passes 156 with no skips;
pinned DBSTEP integration passes 28, and upstream DBSTEP passes 413 with its
three existing skips. Final Memcheck reports no errors or lost bytes in the
declared workloads. Exact CLI output/order equivalence covers 93 cases at
1/2/4/6 threads; final performance preflight also covers 12 threads.

Interrupted rustdoc traversal, the pinned environment's missing pandas, initial
misaligned external conventions, an excluded unpinned adapter build and all
unmet numerical targets are retained as failed attempts with their eventual
resolution or limitation. No failed historical scientific gate is relabelled
a universal pass.

## Accepted commits and files

| Commit | Change |
|---|---|
| `51131ec7` | Freeze baseline and measured profiling before optimization |
| `ab1465be` | U1: `src/density/mod.rs`, borrowed numeric tokens |
| `a286124b` | S1: `src/model/fit.rs`, `src/model/evaluation.rs`, R² correction |
| `e04bc0cf` | U2: `src/geometry/buried_volume.rs`, bounded exact grid reuse |
| `f3658d36` | U3: `src/geometry/buried_volume.rs`, guarded AVX2 row occupancy |
| `81b2bdc9` | Compatible cache syntax and boundary tests |
| `c4dce2d0` | N1: `src/geometry/continuous_sterimol.rs` and public exports |
| `cdc2ddba` | S2: restore actual Rust 1.85 support |
| `f9c952e2` | N2: `src/geometry/continuous_volume.rs` and public exports |
| `06c785da` | Seal final scientific replay, profiling evidence and prepared final measurement controllers |

Use `git show --stat COMMIT` for the full per-commit file set. The prospective
candidate ledger retains hypotheses, acceptance rules, isolated experiments
and the historical rejected row-endpoint optimization. No release, publishing,
GUI or unrelated feature work was performed.

## Reproduction and remaining limits

[BASELINE.md](BASELINE.md) records hardware, compiler flags, source/data hashes
and the original outputs. Raw evidence remains in `.stericx/ultimate_engine/`;
tracked compact receipts bind those files by SHA-256. Each controller uses a
fresh output directory and preserves failed attempts. Reference environments,
source archives, build lockfiles and exact commands remain with their receipts.
The [final scientific evidence seal](final_scientific_evidence.json) binds 65
receipts/logs and verifies unchanged production source. The
[reference review](final_science_review.json),
[exact replay](final_exact_comparison.json) and
[engineering receipt](final_engineering_receipt.json) are also tracked directly.
The [completion audit](COMPLETION_AUDIT.md) maps every requested phase to its
evidence and explicitly identifies the measurements still missing.

After scientific admission and a quiet host, run:

```bash
.venv/bin/python docs/ultimate_engine/run_performance_window.py \
  --output .stericx/ultimate_engine/NEW_FINAL_WINDOW
```

The engine still has binary32 default coordinates/results, phase-sensitive
sampled B1, discrete-grid boundary error, empirical radii and geometric center
assumptions. Density thresholding does not validate an electronic-structure
model. Ensembles depend on conformer coverage, energy provenance and equilibrium
assumptions. Models retain conditioning, feature selection and extrapolation
limits. Continuous-method certification applies to the represented inputs and
documented IEEE arithmetic assumptions, not uncertain molecular geometry.

Only Linux x86-64 native timing is available. ARM emulation supplies functional
checks; real ARM, macOS and Windows hardware tests remain unavailable. The
tested evidence does not establish universal accuracy or fastest-engine status.
