# Campaign work log

- Read the six requested historical reports, existing occupancy/SIMD/grid code,
  continuous chord/interval-box implementation, analytical references and
  benchmark controllers. No applicable AGENTS.md found.
- Confirmed production HEAD `ffc698f4`; existing `docs/media/` is untouched.
- Fresh perf attempt denied at paranoid level 4. Hardware cycles, IPC, native
  branch/cache events and bandwidth cannot be reported. Callgrind/Cachegrind
  are simulated events, as documented in the
  [Callgrind manual](https://valgrind.org/docs/manual/cl-manual.html) and
  [Cachegrind manual](https://valgrind.org/docs/manual/cg-manual.html).
- Freeze and fresh builds precede production changes. Timing gate and numerical
  targets fixed in baseline manifest.
- Fresh baseline: 364 Rust tests and all 128,021 exact historical responses pass.
  Two missing-fixture setup attempts are retained; supplemental fixtures are
  hashed, and no original frozen source/archive was rewritten.
- Fresh Callgrind/Cachegrind/DHAT/Memcheck and assembly completed. Occupancy uses
  115,933,006 of 130,723,468 modeled instructions. No FMA occurs in its assembly.
  All 42 algorithmic-counter/mask replay cases pass.
- Analytical replay: 244 observations, zero validation failures, 23 unmet targets.
  Molecular replay: all twelve still miss 0.05 Å³ at 262,143 boxes.
- Initial timing window fails quiet preflight and is exploratory only. An exact
  JSON comparison overlapped early kernel trials; retained and disclosed.
- O1 overlapping-vector-tail candidate is the first isolated production edit.
  It follows profiling and uses the unchanged frozen baseline binaries for A/B.
- O1 passes 23 focused Rust tests and 80 equal-output alternating native trials.
  Exploratory volume/regional speed ratios are 1.115×/1.135×; unchanged Sterimol
  is slightly slower (0.985×), retained. Final quiet acceptance remains pending.
- C1 chord/Taylor certificate prototype implemented outside production. First
  adapter compilation failed on foreign inherent impl/private error constructor;
  v2 uses adapter-local wrapper types while retaining original interval code.
  Both attempts remain. No certificate has yet been admitted.
