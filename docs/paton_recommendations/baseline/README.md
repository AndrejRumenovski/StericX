# Baseline receipt

Source: `c9380d2c24ba5d48cd52ec13c139ee6e57f04991`, package 0.3.0.
No production source was edited before this freeze.

Fresh commands, run before edits:

```bash
cargo test --locked --all-targets --all-features
.venv/bin/python -m unittest discover -s tests -p 'test_*.py'
.venv/bin/python scripts/replay_scientific_remediation.py build \
  --oracle .stericx/profiling/scientifically_validated_optimization/scientific_oracle_v3 \
  --adapter docs/scientific_remediation/observer \
  --output .stericx/paton/baseline/build
.venv/bin/python scripts/replay_scientific_remediation.py observe \
  --build .stericx/paton/baseline/build \
  --output .stericx/paton/baseline/observations
.venv/bin/python scripts/freeze_paton_baseline.py
```

The test logs retain 322 passing Rust tests and 124 passing Python tests.
The build receipt identifies the exact copied source, dependencies, observer
adapter, executing helpers and build commands. Rust/Cargo are 1.97.0; research
Python is 3.12.13. Installed dependency versions and all source hashes are in
`manifest.json`.

The pre-existing shipping executable is
`b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`.
The freshly built current-source native executable is
`8ed1480f3758e6cf5320549029cf75e10dfd7754e6e9ba6184bb24f8dec95447`;
observer is `f6738ad36017711108cbe291c0b0800eca3c9c6531f5ff46572f5448afcc88ef`.
The old executable predates the optional angular-sampling addition; it must not
be mistaken for the new complete-source binary. Both are retained locally.

The initial output freeze reuses the complete tracked
[angular-sampling observations](../../kraken_residual_forensics/candidates/angular_sampling/full_observations/manifest.json).
The freezer verifies **every Rust source file and both Cargo files** against
that observation provider's source snapshot before reuse, and verifies every
raw output hash. Thus the freeze covers all 13 lanes, all 128,021 rows, their
errors, input IDs, IEEE bits and regional bins without duplicating the large
tracked corpus. The fresh replay is additional evidence; its ongoing log is
excluded from the immutable manifest until completion.

This receipt is not a global scientific PASS. The known limitations, historical
failures, reference-convention differences, missing historical energies and
unfavorable predictive results remain in their original reports. The task plan
indexes what exists and what still requires implementation and validation.

For reproduction, restore the original audit and performance oracle according
to their existing reproduction instructions; use new output directories for
fresh runs. The freezer refuses to overwrite an existing seal. Binaries are
local build artifacts, not committed distribution files.
