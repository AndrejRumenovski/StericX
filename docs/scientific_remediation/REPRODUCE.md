# Reproducing scientific remediation

The coordinated final results are in [SCIENTIFIC_REMEDIATION.md](SCIENTIFIC_REMEDIATION.md)
and the scoped [ADMISSION.json](ADMISSION.json). A completed capture is not
a scientific PASS. Acceptance requires the domain reports, unchanged independent
references, intentionally changed behavior ledger and full engineering checks.
Every destination below must be new; use a fresh suffix for a repeated run.

First restore and verify the original audit payload using its
[reproduction instructions](../scientific_accuracy_audit/REPRODUCE.md). Keep the
pinned Python environment (`uv sync --frozen --all-extras`) and Rust dependency
lock. The original audit, including failed cases, remains immutable.

## Inputs and corrected-system observations

Prepare a separate copy of the historical requests and reference provenance:

```bash
.venv/bin/python scripts/check_current_scientific_equivalence.py prepare \
  --output .stericx/scientific_remediation/requests_v1
.venv/bin/python scripts/replay_scientific_remediation.py build \
  --oracle .stericx/scientific_remediation/requests_v1 \
  --adapter docs/scientific_remediation/observer \
  --output .stericx/scientific_remediation/recheck_v1/build
.venv/bin/python scripts/replay_scientific_remediation.py observe \
  --build .stericx/scientific_remediation/recheck_v1/build \
  --output .stericx/scientific_remediation/recheck_v1/observations
```

The build snapshots and hashes the Rust/Python source, adapter, dependency files
and helper programs, then preserves native and observer executables. The observer
is a measurement adapter calling StericX kernels, not an independent reference.
Its private-volume module contains the current source followed by the unchanged
historical visibility adapter. The additional explicit-topology adapter is
identified separately.

All eleven original request lanes are retained, including the full 31,721
conformer public and per-bin campaigns. Two additional full-corpus lanes call
the new explicit-neighbor entry points using the already frozen source graph.
They retain the same structures, atom indices, IDs, radii and configurations.
They do not silently replace the original coordinate-inference observations.
The capture checks complete ordered IDs, schemas, float/IEEE-bit correspondence
and private frame/bin coverage. Scientific errors remain raw observations for
review. The manifest explicitly does not infer scientific success from execution.

## Unchanged independent reference replay

```bash
.venv/bin/python scripts/replay_scientific_remediation.py references \
  --observations .stericx/scientific_remediation/recheck_v1/observations \
  --components geometry kraken \
  --output .stericx/scientific_remediation/recheck_v1/references_inferred
.venv/bin/python scripts/replay_scientific_remediation.py references \
  --observations .stericx/scientific_remediation/recheck_v1/observations \
  --components kraken --kraken-lane kraken_topology \
  --output .stericx/scientific_remediation/recheck_v1/references_topology
```

The reference Python programs are copied unchanged from the sealed audit.
Package and Morfeus source hashes are checked, scientific inputs are checked
against the original seal, and BLAS/OMP threads are pinned to one. Old
first-orientation reference summaries must not be mislabeled as tests of the
new orientation-mean convention; that convention has a separate independent
per-plane recombination check. Historical published-value residuals, source
limitations and input failures stay visible.

The geometry, thermodynamics and model domain directories contain focused
reproduction programs and retained before/after witnesses. Their results must
be combined with this full-corpus replay and all Rust/Python, Clippy, rustdoc and
format checks. A new optimization baseline is frozen only after that review.

## Fresh Ni-hDA inputs without overwriting historical evidence

```bash
python3 scripts/prepare_remediation_inputs.py \
  --output .stericx/scientific_remediation/ni_hda_inputs_v1
```

This creates a new CSV with response temperature 353.15 K, absolute geometry
paths, explicit provenance notes and hashes of all 67 geometry files. Published
targets, supplied historical 298.15 K populations, energies, split labels and
all other scientific inputs are preserved exactly. It does not rerun MMFF or
CREST or claim that the old populations correspond to the reaction temperature.
The ligand 2064 source conflict remains explicitly unresolved.

Recompute descriptors with the corrected binary and fit into new destinations;
declare `--descriptor-aggregation supplied_weight_mean` only for records actually
built with those supplied means. Use `--response-temp-k 353.15` for Ni-hDA
response metadata. Preserve the checked-in historical CSV, sigpack, study models,
prediction files and unfavorable validation results.

## Complete corrected independent volume/reference coverage

After the two unchanged-reference commands above, use the full derived topology
lane for the independent three-plane mean and every quadrant/octant. Replace the
example roots consistently; never reuse an existing output destination.

```bash
.venv/bin/python docs/scientific_remediation/geometry/recompute_volume_means.py \
  --requests .stericx/scientific_remediation/recheck_v1/observations/derived_requests/kraken_topology_all_bins.jsonl \
  --observations .stericx/scientific_remediation/recheck_v1/observations/kraken_topology_all_bins/stdout.jsonl \
  --expected-rows 31721 --workers 4 \
  --output .stericx/scientific_remediation/recheck_v1/volume_reference
.venv/bin/python docs/scientific_remediation/geometry/compare_all_volume_bins.py \
  --references .stericx/scientific_remediation/recheck_v1/volume_reference \
  --observations .stericx/scientific_remediation/recheck_v1/observations \
  --output .stericx/scientific_remediation/recheck_v1/all_bins
.venv/bin/python docs/scientific_remediation/geometry/replay_fixed_variant_union.py \
  --observations .stericx/scientific_remediation/recheck_v1/observations \
  --inferred .stericx/scientific_remediation/recheck_v1/references_inferred \
  --topology .stericx/scientific_remediation/recheck_v1/references_topology \
  --output .stericx/scientific_remediation/recheck_v1/fixed_union
```

The union retains the original 590/2,009 case selections as well as fresh outliers.
All eleven unchanged reference variants and unavailable statuses remain visible.
See the geometry report for minimized boundary witnesses and separate corrected
mean recombination. A different convention is labeled, not accepted using a
newly enlarged tolerance.

The final thermodynamic source/build binding and numerical scope are in
[thermodynamics/FINAL_V4.md](thermodynamics/FINAL_V4.md). The independent model
commands, original failed folds, additional input-contract tests and training/
screen completeness checks are in [models/REFERENCE_RECHECK.md](models/REFERENCE_RECHECK.md).
The additive CLI comparator requires complete actual artifacts, unchanged helper
and input-plan identities and exact scientific fingerprints. Use `--mode threads`
for one build at different thread counts, or `--mode candidate` for different
builds at the same thread count.

## Corrected profiling inputs and immutable optimization baseline

```bash
.venv/bin/python scripts/prepare_corrected_profile.py \
  --build .stericx/scientific_remediation/recheck_v1/build \
  --inputs .stericx/scientific_remediation/ni_hda_inputs_v1 \
  --root .stericx/profiling/corrected_recheck_v1 --threads 1
```

Preparation preserves the predeclared complete eligible source ensembles and
aborts on any selected-input failure. It creates corrected packed records and a
fresh model/database before performance measurement. Its manifest binds the full
scientific build proof and all helper/source copies before and after execution.

`scientifically_exact_optimization.py freeze` requires a newly reviewed admission
receipt naming the actual build, observations, complete gate evidence roots and
remaining limitations. An old receipt cannot admit a different build. It binds
the exact git commit, all current source, both scientific executables, full raw
evidence and audit payload, workload inputs and machine/toolchain information.
`compare-observations` requires that frozen admission and compares every public
value, IEEE bit, error and private bin exactly. Independent reference replays
and native performance measurements remain separate required gates.
