# Occupancy and continuous-volume baseline

Campaign started 2026-11-17. Production baseline is
`ffc698f479cc5b0d75ae4f6856edc81ad5ca9c2a`. The existing untracked
`docs/media/` is unrelated and untouched. Historical scientific reports,
experiments, outputs and failed targets remain preserved in their original paths.

The source archive, per-file SHA-256 values, historical response/input hashes,
hardware, compiler, environment and prospective benchmark settings are bound by
`baseline_manifest.json` and `.stericx/occupancy_precision/baseline/`.
The freeze/build controller is `campaign.py`; every new run refuses to overwrite
an existing destination. Rust remains the scientific implementation.

The prior 23,783.8 files/s at six threads, 2.1156× total throughput improvement,
2.920× volume kernel gain, 2.918× regional gain, 2.560× conformer gain, 128,021
preserved responses and 81.55% occupancy diagnostic share are historical values.
Fresh profiling and controlled final A/B measurements will establish this
campaign's comparisons; historical timings are not its baseline timings.

Scientific gates: preserve exact finite-grid masks and validated output bits,
including operation grouping and inclusive boundaries. Continuous volume is
the Lebesgue volume of the integration ball intersected with the union of
configured atomic balls. It is a separate quantity and explicit Rust API.
Certified targets are prospectively fixed at 0.1, 0.05 and 0.01 Å³. The twelve
previously failing molecular 0.05 Å³ cases and their 262,143-box budget remain
part of the comparison regardless of outcome.

Final timings use one warmup and eight alternating A/B pairs, identical inputs,
portable release settings and 1/2/4/6/12 threads. Native timings, allocation
diagnostics and simulated instruction/cache events are separate measurements.
The inherited quiet-host gate is retained. If unavailable, all independent work
will finish before requesting a quiet window. No application suspension,
publication, release or GUI work is authorized by this campaign.

Fresh checks: 364 Rust tests pass. All 128,021 responses in 13 lanes exactly
match the previous final engine. The new analytical replay reproduces 244
observations, zero enclosure-validation failures and 23 unmet targets. The
twelve molecular 0.05 Å³ targets still fail at 262,143 boxes; their certified
error bounds span 1.176–3.065 Å³. See `baseline_checks.json` for identities.
Two initial test attempts lacked compile-time/runtime fixtures in the isolated
snapshot. Their logs and the supplemental fixture manifests remain preserved;
the original source archive and all production bytes remain unchanged.

The initial native profiling window is exploratory: its host preflight fails
the fixed quiet gate. An exact-output comparison also overlapped the beginning
of kernel sampling; this is recorded in `baseline_timing_scope.json` locally.
No sample is discarded, and none is promoted to final controlled evidence.


Final controlled comparison completed 2026-11-18 on the same CPU after a reboot
onto Linux 7.0.0-43. Both frozen binaries ran together with unchanged compiler
settings and scientific definitions. Fresh baseline six-thread throughput is
23,411.7 files/s; final is 34,296.4 (1.4649×). Initial exploratory timings remain
separate. Final windows 3/4, all failed attempts, and the kernel transition are
recorded in timing_admission.json and FINAL_REPORT.md.
