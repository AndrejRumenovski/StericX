# Final campaign completion audit

Completed 2026-11-17 on the available host. This checks the requested phases
against actual evidence and preserves the [pre-timing audit](COMPLETION_AUDIT_BEFORE_TIMING.md).
The campaign has no pending applicable validation, measurement or implementation
work. Unavailable hardware tests/counters and scientifically unadmitted comparisons
are explicit limitations, not passes or work silently substituted by another task.

| Phase | Evidence and outcome | Limits retained |
|---|---|---|
| 0: inspection and freeze | BASELINE.md, source/data/tool manifests, frozen original `c03630cb`; baseline outputs and measurements predate edits | Original working-tree state and historical failures preserved |
| 1: scientific audit | SCIENTIFIC_AUDIT.md covers Sterimol, total/regional volume, ensembles, thermodynamics, density, pyramidalization, radius profiles, models and conventions; R² defect corrected | Empirical chemistry and numerical limitations remain explicit |
| 2: comprehensive profiling | Original/final native profiles, modeled instruction flamegraphs, allocation/DHAT/RSS, assembly/dynamic SIMD, 1/2/4/6/12 threads, Python layers, accuracy/cost | Hardware cycles/IPC/cache misses and bandwidth denied/unavailable; LLVM remarks absent; modeled events not mislabelled |
| 3: numerical accuracy | Explicit N1 continuous support and N2 midpoint/reference Rust APIs, independent analytical/high-precision/interval references, adversarial cases and bounds | Default unchanged; non-certified methods identified; 23 stress and all 12 strict molecular targets remain unmet |
| 4: Rust optimization | Isolated U1 tokens, U2 exact grid reuse, U3 SIMD masks, separate commits, independent reference gates and Memcheck | Cache peaks can rise; no unvalidated floating-point reassociation |
| 5: joint accuracy/cost | Final 234 support and 72 volume configurations, per-input Pareto JSON and standalone reviewed figure, relative/absolute residuals and certified bounds | Observed residual differs from certification; no universal optimum |
| 6: independent validation | Full Morfeus/DBSTEP/Kraken/analytical replay; 128,021 exact legacy responses; regional, thermodynamic, density and numerical evidence | Historical convention failures remain; physical ARM/macOS/Windows unavailable |
| 7: competitive benchmarks | 432 admitted volume launches, 1/56/1,000 input workflows, startup/warm timing, scaling and separate process-tree memory; 36 Sterimol cost/accuracy launches | DBSTEP common-grid scope explicit; two external Sterimol targets fail and no equal-accuracy ratio admitted; no controlled disk-cache-cold claim |
| 8: optimization loop | Measured isolated candidates, three retained performance changes, rejected row-endpoint experiment preserved, final reprofiling | Further occupancy/layout/transport work requires a new measured cycle; no insignificant tuning after full replay |
| 9: full engineering/scientific replay | 364 native and 364 MSRV Rust tests, 362 emulated ARM tests, 71 binding, 156 research, 28 pinned integration, 413 upstream tests; fmt/Clippy/docs/Memcheck | Three upstream skips retained; workload tests do not prove all-platform memory safety |
| 10: final A/B | Window 6 completes all eight stages with every required quiet endpoint true; 630 CLI, 80 kernel A/B, 71 kernel profiles, 175 CLI profiles, 208 Python A/B; exact checks pass | Every unfavorable trial retained; host endpoints cannot rule out all intermediate scheduling noise |
| 11: report and evidence | All requested reports, complete timing distributions/tables, accepted/rejected commits, source/binary/input hashes, reproduction commands and final seal | Claims restricted to measured definitions, accuracy and hardware |

## Final success-criteria disposition

Every implemented scientific family is audited. Confirmed implementation defects
within the campaign were corrected, with independent witnesses. New numerical
methods have independent validation and explicit accuracy contracts. Available
profiling alternatives were exhausted without inventing denied hardware counters.
Measured high-impact bottlenecks received three accepted optimizations, each
validated independently and in the final exact replay. Competitive speed claims
are restricted to the admitted equivalent finite-grid volume task; failed
Sterimol admission is reported rather than hidden. The final frozen-baseline
comparison, reproducible reports and evidence seal are complete. Rust remains
the authoritative scientific implementation.

The final six-thread 10,000-file result is 23,783.8 versus 11,242.0 files/s,
with exact output preservation. The final-report claim is scoped campaign
completion, not universal accuracy, hardware coverage or fastest-engine status.

## Final integrity checks

`seal_final_campaign.py` checks the unchanged production diff from `f9c952e2`,
rehashes all 65 scientific records, verifies complete final measurement receipts
and passing host endpoints, binds every regular raw file in the accepted timing
directories, and binds the final reports, tables and plots. Its `--verify` mode
rehashes the resulting manifest without modifying it. The exact pre-run-hashed
controller archive remains available, including the older plotting controller.
The benchmark service has ended, both authorized UI suspensions have a successful
resume receipt, and no watchdog/benchmark process remains running.

No release, upload, GUI or unrelated feature work was performed. Existing
untracked `docs/media/` content was not included in campaign commits.
