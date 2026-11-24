# Reproducing this campaign

The production baseline is commit `e0f1dc79d23bf6638614e85cfc0c9f2282b6b6eb`.
[BASELINE.md](BASELINE.md) and `baseline_manifest.json` bind the original source,
constants, toolchains and historical outputs. Raw evidence is retained under
`.stericx/experimental_validation/`; it includes failed attempts. No historical
file is replaced. Experimental inputs and principal outputs are additionally
versioned in the [deterministic bundles](data/README.md).

## Restore and inspect experimental data

```sh
.venv/bin/python docs/experimental_validation/archive_chemical_data.py restore
```

Restoration checks every SHA256, refuses different existing bytes and permits only
paths within the campaign evidence directory. Original SDF headers/coordinates,
identities, source failure responses, generated candidates and actual computed
energies are preserved. The original and corrected SI PDFs can be re-downloaded from
recorded primary URLs; the primary table was visually and textually compared.
Source DOI, retrieval time, status, file identity and usage terms are in the receipts.

`results/experimental_observations.csv` is the SI transcription; `primary_design.csv`
contains the locked experimental modeling rows. `conformer_descriptors.csv` contains
native per-state descriptors, not invented measurements. `ensemble_statistics.json`
contains all eleven-feature reductions and generated energy distributions.
`statistical_results.json` contains every outer/inner split, prediction, failed-fold
list, model score, coefficient or hyperparameter and uncertainty result. Supplemental
stereo and numerical-method models have separate files.

## Authoritative calculation and independent analysis

The research adapter uses the production Rust library by path. Build with the locked
dependency graph, preserving old binaries and writing new observations separately:

```sh
CARGO_BUILD_JOBS=2 cargo build --locked --offline --release \
  --manifest-path docs/experimental_validation/observer/Cargo.toml \
  --target-dir /tmp/stericx-experimental-reproduction
```

The `chemistry` binary accepts a JSON request array. The archived requests are
`chemistry_preparation/requests.json` and `generated_v2/requests.json`. Each names
original files and, only in the generated lane, explicit per-state energy metadata.
Use a new output directory; do not redirect over archived `native.jsonl` evidence.
Elapsed times vary, while descriptors/certificates and deterministic work counts are
compared independently of timing and allocation diagnostics.

- `prepare_chemistry.py`, `audit_stereochemistry.py`: source identity and mapped P audit.
- `generate_conformers.py`: corrected, separate ETKDGv3/MMFF94 preparation. Seeds,
  budgets, pruning, labels, rejected structures and real energies are serialized.
- `analyze_chemistry.py`: independent reductions, locked fits, nested ridge, conditional
  intervals, resampling and all model sensitivities. Python fits models; Rust supplies
  molecular descriptors, ensemble reductions and thermodynamic weights.
- `validate_chemical_references.py`: independent Morfeus/DBSTEP envelopes and SciPy fits.
- `conformer_sensitivity.py`, `stereo_sensitivity.py`,
  `chemical_numerical_sensitivity.py`: disclosed sensitivity analyses.
- `render_chemical_report.py`: derived tables and standalone SVG/PNG scientific figures.

Controllers normally refuse pre-existing output directories. Their default names
identify the retained run, not an invitation to overwrite it. To repeat in isolation,
use a fresh checkout/evidence root or the separate-output replay described in the
final validation report. The archived first generation controller is deliberately
preserved alongside the corrected controller.

## Numerical methods and regression

`validate_numerical.py` evaluates all 167 frozen analytical requests, including the
eighteen original failures, against the original independent reference intervals.
`validate_regions.py` derives 96 regional cases using exact rational cross-section
integration and a rational Machin-π enclosure. `validate_profiles.py` derives exact
analytical profiles and checks the frozen molecular geometries, standalone identity,
work accounting and monotonicity consistency. The original tolerances are unchanged.

`tests/certified_geometric.rs`, `tests/certified_regions.rs` and
`tests/certified_profiles.rs` exercise public APIs; the exact-addition test reads
6,005 independently generated rational-oracle witnesses. Full engineering checks are
logged by `run_engineering.py`, including actual Rust 1.85, modern Rust, Python,
research, external DBSTEP and available cross-target compilation.

`scripts/replay_scientific_remediation.py` builds a source-bound observer and captures
all thirteen historical lanes. The campaign replay adapter only supplies unused new
certified types/error bridges required to compile the original private-module prefix;
the historical observation path and appended private probes are unchanged.
`docs/ultimate_engine/compare_replay.py` requires exact parsed JSON agreement for all
128,021 outputs. Independent reference controllers use the sealed original equations
and conventions, never the optimized native kernel as their reference.

## Timing and portability limits

All new times are exploratory measurements on the ordinary shared workstation.
`profile_methods.py` records CPU/wall time, RSS, allocation traffic, work counts,
expanded budgets and deterministic parallel outputs. `profile_valgrind.py` uses a
whitespace-free tool alias because `LD_PRELOAD` splits paths containing spaces; the
initial invalid probe remains recorded. Callgrind counts simulated instructions,
not hardware events. Hardware counters are denied by `perf_event_paranoid=4`; no
system permissions were changed. No controlled speedup or new release is claimed.

Package/toolchain versions and binary/source hashes are retained. Reproducing
floating-point certificates assumes IEEE binary64 round-to-nearest with gradual
underflow; geometry generation additionally depends on the specified RDKit version.
Cross-target compilation is not a physical cross-platform accuracy/timing experiment.
