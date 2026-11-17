# Campaign work log

## Current state

Production implementation and full scientific replay are complete. Final source
is `f9c952e2`; all 128,021 baseline responses match exactly. Independent full
reference review finds zero changed historical scientific values/statuses,
including 31,721 volume calculations and both fixed Kraken variant unions.
Rust, Python, pinned references, MSRV, emulated ARM, memory and new numerical
method checks are complete in the scopes recorded in VALIDATION_REPORT.md.
Final controlled timing remains pending; no final files/s claim is made.

The baseline archive contains 220 source/tool files. Dataset hashes cover 557
tracked files; historical-report hashes cover 174 Markdown files. Full tracked
inventory has 21,899 paths. All old failures remain untouched.

Read the requested final reports, the remediation and audit conclusions,
Kraken residual conventions, the rejected C3 row-endpoint experiment, profiling
implementation and core geometry, ensembles, density and profile calculations.
Historical C3 reduced predicate counts but slowed native execution by about 15.6%;
do not repeat that claim as a speedup. Historical DBSTEP mismatch was a second
sphere clip, not an atomic-occupancy defect. The corrected common-node oracle
must be replayed separately from the preserved failed historical lane.

Valgrind 3.26.0 works via `/tmp/stericx-ultimate-valgrind`, a symlink to the retained
local installation. Set `VALGRIND_LIB` to its `usr/libexec/valgrind` and use
`usr/bin/valgrind.bin`; spaces in preload paths otherwise break the loader.
Hardware perf is denied at paranoid=4; no privilege change was attempted.

## Remaining work

1. Sample actual host activity after all independent work has finished. If busy,
   request the user's stipulated quiet window; do not weaken the admission gate.
2. Run the prepared sequential final A/B, profiling, Python, accuracy-cost and
   scientifically matched competitive measurement window.
3. Review every distribution, report unfavorable results, finish final performance
   tables and figures, seal evidence and commit the completed reports.

All seven requested report paths now exist. FINAL_REPORT.md explicitly marks
pending controlled measurements. The campaign remains active; it is not fully
complete until the remaining measurement and reporting work is handled.

## Quiet-window gate

After all independent work and commit `06c785da`, the 10.012-second host sample
recorded 0.351 busy core equivalents (limit 0.15), with Warp 0.149 and GNOME
Shell 0.118. I/O wait was 0.036 (limit 0.05). No final performance measurement
started. A quiet-window confirmation was requested under the user's Phase 10.
The goal remains active; elapsed time is not confirmation. Once the user confirms,
recheck the host and run the prepared controller into a fresh output directory.
Raw and tracked host evidence are preserved as `host_after_science.json`.

## First resumed host check

The user reactivated the blocked goal. A fresh 10.013-second sample still
failed admission: 0.380 busy CPU-core equivalents, including Warp 0.163 and
GNOME Shell 0.117; I/O wait was 0.004. This is a fresh blocked-audit sequence,
not an extension of the previous three-turn count. No benchmark was started.
The user was asked to minimize Warp and leave the computer idle, then confirm
readiness. Original and resumed samples are preserved separately; the quiet
criterion and all scientific/performance protocols remain unchanged.

The user subsequently confirmed that Warp was minimized and the host would be
idle. `final_window_1` stopped in its first admission check (0.286 busy cores),
before any timing. Minimized Warp still used 0.171 cores; GNOME Shell fell to
0.023. The user was asked to authorize temporary suspension of Warp's UI only,
with automatic resume and an independent 30-minute watchdog. This is an
additional application-control permission, not a repeat request for the already
confirmed idle window. No user process has been signalled. Synthetic-process
normal/failure/controller-crash checks all passed; approval remains pending.


## Final controlled windows and campaign closure — 2026-11-17

The user explicitly authorized Warp UI suspension and continuation through the
remaining goal. The earlier pending-approval/status entries above are historical.
Window 6 completes all eight final stages with passing quiet endpoints; all
630 CLI A/B, 80 kernel A/B and 208 Python A/B trials preserve expected scientific
results. Competitive volume completes 432 exact-result launches. Controlled
costs cover 39 support cases/234 configurations and eleven volume geometries/72
configurations. Final source remains `f9c952e2`. Warp was restored.

The final coverage audit added a prospective twelve-molecule Sterimol supplement
(commit `ded41cc6`). Both native methods pass its fixed 0.01 Å target; Morfeus
and DBSTEP each miss one B1 target. All cases and 36 cost trials remain, without
equal-accuracy ratios for failed methods. Its short guarded window passed both
quiet gates and restored Warp; no benchmark process remains.

Final six-thread throughput is 23,783.8 versus 11,242.0 files/s (+111.56%); one
thread is 4,503.3 versus 2,080.3, and twelve threads 28,320.8 versus 13,941.4.
All raw ranges and unfavorable results remain. Final tables, reviewed scientific
plots, reports, phase completion audit and evidence seal are complete. No source
changes after the scientific replay required another production test run.
Hardware-counter and physical cross-platform limitations remain explicit.
