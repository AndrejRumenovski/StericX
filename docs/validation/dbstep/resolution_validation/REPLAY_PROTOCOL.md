# Fresh scientific replay for DBSTEP discrepancy resolution

This campaign reuses the retained scientific controllers and numerical comparison
rules. It checks preservation of legacy scientific results, including negative
results. Successful execution does not admit a release or establish agreement
between methods with different definitions.

The source-bound build is
`.stericx/dbstep_resolution/baseline/build/manifest.json`. Its native executable
SHA-256 is
`3b9079338386a2532b17ebe844a28d72344a7bdf91d83526d4b3d0fd00438424`;
the observer SHA-256 is
`6b3a5ccf0c8d95275f4c32b3a96207fc4beaecab9eeab054d257161abf24dcef`.
Both match the previous final scientific executables. The observer is a
measurement adapter, not an independent physical reference.

## Commands

Every output destination must be new. The campaign below uses `science_v1`;
replace that suffix consistently for a repeat. The build path identifies the
source and executables being evaluated and must be verified before replay.

```bash
.venv/bin/python scripts/replay_scientific_remediation.py observe \
  --build .stericx/dbstep_resolution/baseline/build \
  --output .stericx/dbstep_resolution/science_v1/observations
.venv/bin/python scripts/compare_paton_replays.py \
  --baseline .stericx/paton/final_observations/manifest.json \
  --candidate .stericx/dbstep_resolution/science_v1/observations/manifest.json \
  --output .stericx/dbstep_resolution/science_v1/legacy_comparison.json
```

Run the new orchestration script once for each phase: `engineering`, `cli`,
`models_thermo`, `references`, and `volume`. The model phase requires complete CLI
and observation receipts; references require observations, and the volume phase
requires the completed references. For example:

```bash
.venv/bin/python scripts/replay_dbstep_resolution_science.py \
  --build .stericx/dbstep_resolution/baseline/build \
  --output .stericx/dbstep_resolution/science_v1 \
  --logs docs/validation/dbstep/resolution_validation \
  --phase references
```

Each phase records exact subprocess arguments, working directory, start/end UTC,
elapsed duration, return code, and stdout/stderr hashes in
`science_v1/execution/<phase>/journal.json`. Executed orchestration bytes are
archived beside that receipt. The first engineering runner originally bound its
live path; `controller_archive.json` binds the byte-identical archived copy before
later orchestration changes. Its successful journal and commands remain unchanged.
The future fixed-union command was corrected to name each reference's `kraken`
subdirectory before that phase ran. No scientific controller or tolerance was
changed.

The archived Kraken helpers hard-code three process workers. Their source is
retained unchanged, and the reference phase is restricted to CPUs 0 and 1. The
full-volume and fixed-union programs accept an explicit worker count and use two.
BLAS, OMP and MKL threads are one. These are scientific replays, not controlled
performance measurements.

After all phases finish, review them independently:

```bash
.venv/bin/python scripts/review_dbstep_resolution_science.py \
  --build .stericx/dbstep_resolution/baseline/build \
  --science .stericx/dbstep_resolution/science_v1 \
  --output .stericx/dbstep_resolution/science_v1/review
```

The review checks raw source/build bindings, complete memberships and coverage,
unaltered reference artifacts, fixed union membership, regional results, model
equations, thermodynamics, actual CLI comparisons, and engineering receipts. It
uses the retained comparison implementation and explicitly leaves global
scientific and release admission false. New DBSTEP mathematical analysis,
convergence, advanced Python methods, wheel portability, and performance require
their separate reports.

## Retained scientific scope

The observations retain all 128,021 records in 13 lanes, including unsuccessful
geometries and the distinction between inferred connectivity and explicit source
topology. The full reference campaign covers 31,721 conformers, 95,163 frames and
1,141,956 quadrant/octant bins. The fixed union retains 2,517 cases in each of two
connectivity lanes. The actual 93-case CLI campaign runs at 1, 2, 4 and 6 threads.

Existing finite-angle B1 orientation differences, radius-table strata,
historical Kraken provenance limitations, failed model folds, unfavorable ranking
and predictive-calibration findings remain visible. Numerical preservation does
not upgrade those findings. The original strict DBSTEP hard-sphere failures are
also retained; any revised validation policy belongs to the separate prospective
DBSTEP resolution analysis.
