# Current-source scientific observer

The initial build in `.stericx/paton/dbstep_build` failed because the frozen
historical observer manifest declared serde only transitively through StericX.
The newly serializable native volume types now name `serde` directly in the
source copied into the observer. The compiler diagnostics are retained in
`native_build_failure_stderr.log`; no numerical comparison was made with that
failed build.

The replay builder now adds a direct serde dependency only when the copied
source uses it. Its version is taken from the current source lock, and only the
observer package's dependency edge is added to the copied lock. Registry crate
versions and checksums remain unchanged. The historical manifest, lock,
observer append, reference code and audit evidence are not edited. Two focused
tests verify the old template is byte-identical for old source and all locked
crate records are preserved for current source.

The fresh build in `.stericx/paton/dbstep_build_v2` succeeds. Its
[manifest](native_build_manifest.json) binds the complete copied source,
byte-identical native volume prefix, unchanged observation append, dependency
adjustment, toolchain, commands, helper sources and executing native/observer
hashes. The build also verifies the live sources have not changed during
compilation. `observer_build_success.log` and `native_build_success.log` preserve
the successful builds.

```bash
.venv/bin/python -m unittest discover -s tests -p test_observer_dependency_binding.py -v
.venv/bin/python scripts/replay_scientific_remediation.py build \
  --oracle .stericx/profiling/scientifically_validated_optimization/scientific_oracle_v3 \
  --adapter docs/scientific_remediation/observer \
  --output .stericx/paton/dbstep_build_v2
```

Use a fresh output directory for replay. This observer captures the system
under investigation; DBSTEP, Morfeus and independent mathematical references
provide validation. A successful build is not a scientific correctness gate.
