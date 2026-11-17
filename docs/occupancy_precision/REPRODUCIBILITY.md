# Reproducing the occupancy/precision campaign

Run commands from the repository root. Local raw evidence is under
`.stericx/occupancy_precision`; tracked receipts and controllers are in this
directory. Existing output directories are never overwritten. Use a new name
for a replay. The frozen baseline is commit
`ffc698f479cc5b0d75ae4f6856edc81ad5ca9c2a`, with its original source archive,
per-file hashes, data identities, compiler and host receipts in
`baseline_manifest.json`. The historical reports remain unchanged.

The command examples below name the actual retained run directories. They are
an execution record, not instructions to overwrite them. A fresh rerun must use
new output names or a separate copy of the workspace. Some reference controllers
expect the retained historical `.stericx` datasets and dependency caches; this
workspace contains those inputs, but a bare Git clone does not.

## Frozen builds and exact finite-grid replay

```sh
.venv/bin/python docs/occupancy_precision/campaign.py freeze
.venv/bin/python docs/occupancy_precision/campaign.py build --name baseline
.venv/bin/python docs/occupancy_precision/campaign.py build --name final_profile
CARGO_BUILD_JOBS=2 .venv/bin/python scripts/replay_scientific_remediation.py build \
  --oracle .stericx/profiling/scientifically_validated_optimization/scientific_oracle_v3 \
  --adapter docs/scientific_remediation/observer \
  --output .stericx/occupancy_precision/final_scientific_build
.venv/bin/python scripts/replay_scientific_remediation.py observe \
  --build .stericx/occupancy_precision/final_scientific_build \
  --output .stericx/occupancy_precision/final_science/observations
.venv/bin/python docs/ultimate_engine/compare_replay.py \
  .stericx/occupancy_precision/baseline_observations \
  .stericx/occupancy_precision/final_science/observations \
  .stericx/occupancy_precision/final_exact_comparison
```

The fresh baseline replay was independently compared with the previous campaign's
final 128,021 responses before any production edit. Baseline fixture supplements
and two failed setup attempts are retained; they did not rewrite production
source or the original archive.

## Scientific and engineering checks

```sh
CARGO_BUILD_JOBS=2 cargo test --locked --offline --all-targets --all-features
cargo clippy --locked --offline --all-targets --all-features -- -D warnings
.venv/bin/python docs/occupancy_precision/engineering_replay.py \
  --build final_profile --output final_engineering
.venv/bin/python docs/occupancy_precision/arm_replay.py \
  --build final_profile --output final_arm
.venv/bin/python scripts/replay_dbstep_resolution_science.py \
  --build .stericx/occupancy_precision/final_scientific_build \
  --output .stericx/occupancy_precision/final_science \
  --logs .stericx/occupancy_precision/final_science/logs --phase references
.venv/bin/python docs/occupancy_precision/volume_replay.py compute
.venv/bin/python docs/occupancy_precision/volume_replay.py compare
.venv/bin/python docs/ultimate_engine/review_final_science.py \
  --build .stericx/occupancy_precision/final_scientific_build \
  --science .stericx/occupancy_precision/final_science \
  --output .stericx/occupancy_precision/final_science_review
```

The CLI and models/thermodynamics phases use the same replay controller with
`--phase cli` followed by `--phase models_thermo`. Volume comparison waits for
both volume computation and both reference lanes. Reviewing unchanged known
reference discrepancies is separate from checking process exit codes.

Rust 1.85 and target checks use installed toolchains. AArch64 functional tests
use the retained QEMU executable and a private CLI-child shim; they are not
physical ARM benchmarks. Python wheels are built and extracted locally for
tests; nothing is published or installed into the user's global environment.

## Certified integration

```sh
.venv/bin/python docs/occupancy_precision/build_certified_observer.py c1_public_v1
.venv/bin/python docs/occupancy_precision/validate_chords.py \
  --engine .stericx/occupancy_precision/c1_public_v1/bin/certified \
  --output c1_public_extended
.venv/bin/python docs/occupancy_precision/validate_predicates.py \
  --build c1_public_v1 --output c1_predicates_v1
.venv/bin/python docs/occupancy_precision/validate_coaxial.py \
  --engine .stericx/occupancy_precision/final_profile/bin/certified \
  --output c1_coaxial
.venv/bin/python docs/occupancy_precision/validate_molecular_certified.py \
  --engine .stericx/occupancy_precision/c1_public_v1/bin/certified \
  --reference-engine .stericx/occupancy_precision/baseline/bin/continuous_volume \
  --output c1_molecular_broad_v3
```

Original analytical requests, original twelve molecular requests and the expanded
budget requests are frozen JSON inputs with hashes in their validation receipts.
The old continuous midpoint/reference API is replayed independently using
`docs/ultimate_engine/validate_continuous_volume.py`; its original targets stay
unchanged. Prototype, public-API and ablation binaries are separately identified.
Integer-square-root/rational references do not reuse production predicates.
Coaxial references support higher-order overlaps; triangle-free references
prove their overlap restriction before using pairwise inclusion–exclusion.

## Profiling and timing admission

`profile_campaign.py --phase valgrind` collects Callgrind, Cachegrind, DHAT and
Memcheck against a named frozen build. The retained local Valgrind distribution
needs `VALGRIND_LIB=/tmp/stericx-ultimate-valgrind/usr/libexec/valgrind` and the
documented symlink to its preserved tools directory. Full call-context profiles
use `--separate-callers=40`; `callgrind_flamegraph.py` renders exclusive modeled
instructions. Static assembly comes from `objdump -Cd -Mintel` on the exact
bound observer binary. None of these yields native cycles, IPC or bandwidth.

The matched competitive gate uses the unchanged exact-decimal/IEEE reference,
Morfeus public buried volume and DBSTEP's occupied-point kernel on common nodes.
The separate radius-profile gate fixes twelve conformers, radii 2.5/3.5/4.5 and
`density=0.01*(R/3.5)^3` before observing results. It is a common-32-side lattice
profile, not a constant-density profile. Both gates bind inputs and binaries.

```sh
.venv/bin/python docs/occupancy_precision/accuracy_campaign.py preflight --output accuracy_preflight
.venv/bin/python docs/occupancy_precision/accuracy_campaign.py memory --output accuracy_memory
.venv/bin/python docs/occupancy_precision/radius_benchmark.py preflight --output competitive_radius_preflight
.venv/bin/python docs/occupancy_precision/performance_window.py --output final_window
```

The final controller requires scientific admission and untimed preflights. It
retains the predeclared 0.15-busy-core/0.05-I/O-wait-core quiet gate and stops on
failed endpoints. It never stops unrelated applications. Every timing attempt,
including failed admission and unfavorable samples, remains recorded. Warm
API/workflow timing and fresh-process timing are separate scopes; OS file caches
are warm, and no cold-disk claim is made. Allocation and memory diagnostics are
separate from native speed measurements.

`EXPERIMENTS.md` identifies rejected implementations, setup failures and methods
that were analyzed but not implemented. Source archives, hashes, compiler
commands and raw observations are the reproducibility basis; historical speed
figures are context rather than a replacement for fresh controlled A/B data.
