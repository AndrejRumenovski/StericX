# Campaign work log

## Current state

Baseline freeze, builds, tests and fresh profiling are complete. The full
128,021-response baseline observation replay completed. U1 density parser
candidate passed exact 23-case replay, eight focused tests, 12 independent
DBSTEP comparisons, Memcheck and an eight-pair exploratory timing test; retain
for full final replay. R² scale-invariance regression witnesses are being tested
against unchanged metric code before correction.

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

## Remaining phases

1. Finish baseline builds/tests, observe current scientific values and run fresh
   CPU/heap/descriptor/SIMD/thread/Python profiling before production optimization.
2. Finish calculation-by-calculation audit, including kinetics/models and Python
   data transfer; investigate defects with independent witnesses.
3. Evaluate explicitly distinct numerical methods, including exact support
   envelopes and continuous sphere-union integration with documented bounds.
4. Implement isolated measured optimizations, validate, benchmark, retain/reject,
   reprofile and commit accepted improvements separately.
5. Run the full existing reference replay (128,021 observations and independent
   Morfeus/Kraken/DBSTEP/thermodynamic/density checks), new tests, memory checks,
   Python tests and available architecture checks.
6. Build accuracy/runtime Pareto evidence and matched competitive comparisons.
7. Complete final controlled A/B measurements (quiet-window request only if
   needed and only after non-timing work), all seven requested reports, evidence
   integrity/reproduction review and requirement-by-requirement completion audit.

Goal remains active until the full requested scope is proved or a real external
blocker satisfies the goal's blocked threshold. No partial-completion claim.
