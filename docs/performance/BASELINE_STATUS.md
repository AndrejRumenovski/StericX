# Corrected optimization baseline

The current baseline is corrected commit
`1ecfbd5be8b354711bdadd8a8d148417e03d8e7a`. Its independently rechecked native release
executable SHA-256 is
`102b6883c639bfc9ef211a88abae02391419cbd277c498098fa2bdfc7008206f`.
The [remediation report](../scientific_remediation/SCIENTIFIC_REMEDIATION.md) and
[scoped admission](../scientific_remediation/ADMISSION.json) record the corrected
scientific contracts and remaining limitations. The original audit remains unchanged.

The immutable freeze is
`.stericx/profiling/scientifically_validated_optimization/accurate_baseline_v1/`.
It binds the source commit, Cargo/Rust/Python source, native and observer binaries,
complete outputs/errors, independent reference evidence, machine/toolchain,
76,714 raw evidence files and corrected benchmark inputs.

The baseline gate includes 128,021 complete observations, all 31,721 topology
conformers and 1,141,956 regional bins, 7,624 model comparisons, the analytical
thermodynamic replay, 93 CLI cases at each of 1/2/4/6 threads, 306 Rust tests,
113 Python tests, four reference-integrity tests, Clippy, rustdoc and formatting.
These counts establish stated coverage; they do not establish universal chemical
or experimental validity.

Fresh profiling and optimization from this baseline retained ordered descriptor
file parallelism (C1/C1b) and exact Student-t multiplier reuse (C2). The final
native SHA-256 is
`b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`.
Fresh scientific/reference/engineering gates and the final 720-launch matrix
passed. The buried-volume row-rejection experiment C3 was exact but slower and
was reverted. The [final report](SCIENTIFICALLY_EXACT_OPTIMIZATION.md) records the
accepted gains, remaining regressions, all evidence and reproduction instructions.

The earlier audited `b515c4f` freeze and its withdrawn file-parallel experiment
remain preserved under `current_audited_b515c4f/` and in the
[historical checkpoint](PRE_REMEDIATION_OPTIMIZATION_CHECKPOINT.md).
Those timings and negative scientific outputs are not the corrected baseline.
