# Corrected native paired measurements

Run only after the candidate has completed the scientific and engineering gates:

```sh
.venv/bin/python scripts/benchmark_corrected_pair.py \
  --snapshot .stericx/profiling/scientifically_validated_optimization/accurate_baseline_v1 \
  --candidate-build /absolute/new-candidate-build/manifest.json \
  --gate /absolute/reviewed-candidate-gate/manifest.json \
  --output /absolute/new-paired-measurements
```

The output directory must not exist. There is no force, reduced-check, changed
repetition, or tolerance option. All ten frozen workloads run at 1/2/4/6 threads
with their frozen CPU affinities. Each configuration starts with one baseline
and one candidate warmup, then eight measured pairs in AB/BA/AB/BA order. This
means **40 configurations, 80 warmups, 320 pairs, 640 measured processes, and 720
total processes**. Run this phase without competing CPU-intensive work.

The baseline executable is `SNAPSHOT/bin/native`; the candidate executable is
the native binary in its verified build receipt. Each launch gets a separate
output directory. Existing `profile_stericx` measurement, launcher, result checks,
and distributions are reused. Stderr is additionally compared exactly. Only the
existing specifically named timing/resource/output-path fields are excluded
from stdout/DB-manifest fingerprints; packed records, DB tables, and all emitted
scientific outputs retain exact fingerprints. There is no new numeric tolerance.

Prediction stdout contains previews and aggregates. A mandatory separate gate
compares all million f32 predictions at each thread count; the benchmark does not
claim that every timed CLI process emits that complete vector.

Every attempted process, raw stdout/stderr/artifact, metrics, warmup and failure
is retained. A mismatch stops further launches and fails the result. No outlier
is removed. Summaries give native wall/CPU/RSS/fault distributions, median/MAD/
range, all eight within-pair wall speedups (`baseline/candidate`), and raw samples.
Input, source/build, executable, helper, launcher, scientific evidence, and
successful output hashes are checked before/after the entire phase. Large frozen
inventories are not rehashed between individual launches. A final inventory also
exists for failed attempts. Postflight failure makes the overall result fail.

## Required gate records

A **record** below is `{"path": ABSOLUTE_PATH, "bytes": INTEGER, "sha256": HEX}`.
All gate/comparison/build/capture receipts have the existing `.json.sha256`
sidecar. Referenced plain result JSON, logs, bitstreams, and native baseline
summaries need only their record. Use `replay_scientific_remediation.manifest` to
write new sealed JSON; preserve all earlier evidence directories.

The top-level schema is:

```json
{
  "kind": "corrected_candidate_performance_gate",
  "complete": true,
  "baseline_snapshot": "RECORD of frozen manifest.json",
  "baseline_build": "RECORD of the original build (copy.original for receipts/build.json)",
  "candidate_build": "RECORD of candidate build manifest.json",
  "observation_comparison": "RECORD of exact_corrected_baseline_comparison",
  "cli_comparisons": {"1": "RECORD", "2": "RECORD", "4": "RECORD", "6": "RECORD"},
  "independent_references": "RECORD of reviewed reference gate described below",
  "engineering": "RECORD of engineering gate described below",
  "prediction_comparison": "RECORD of bitstream comparison described below",
  "baseline_measurements": {"1": "RECORD", "2": "RECORD", "4": "RECORD", "6": "RECORD"}
}
```

The strings marked `RECORD` above are placeholders for record objects. The
scientific gate must be reviewed; this helper verifies its evidence and concrete
pass conditions but does not rerun the independent equations itself.

`observation_comparison` is the existing exact observation comparator output:
all thirteen lanes, nonzero complete matching row counts, zero mismatches, exact
stderr and raw observation stdout, zero subprocess exits, zero tolerance, known comparison helper, and the
admitted frozen observation baseline. Linked observations must name the exact
baseline/candidate builds and observer binaries. All referenced raw files are
rehashed. `cli_comparisons` uses `compare_corrected_cli.py --mode candidate` for
each thread count; every receipt must pass all 93 cases exactly, contain no
differences, and link both builds and their native binaries. The underlying
capture fingerprint maps are recomputed from actual raw captures against the
canonical 93-case plan; required successful artifacts must remain nonempty.

`baseline_measurements` points to the completed original `profile_stericx`
native `summary.json` for each thread count. It must contain all ten workloads,
one warmup, at least seven measured repetitions, the frozen input manifest and
native SHA, and the identical thread settings/affinities. Every original raw
sample is checked against its recorded output fingerprint and stdout/stderr
hashes. Every new A/B launch must match these original fingerprints.

The reference review wrapper has this schema:

```json
{
  "kind": "corrected_candidate_independent_reference_gate",
  "complete": true,
  "passed": true,
  "baseline_build": "RECORD",
  "candidate_build": "RECORD",
  "domains": {
    "geometry": {"passed": true, "evidence": ["RECORD"], "checks": [
      {"evidence": "SAME RECORD", "field": ["reviewed_scope_passed"], "equals": true}
    ]},
    "thermodynamics": {"passed": true, "evidence": ["RECORD"], "checks": [
      {"evidence": "SAME RECORD", "field": ["failures"], "equals": []}
    ]},
    "models": {"passed": true, "evidence": ["RECORD"], "checks": [
      {"evidence": "SAME RECORD", "field": ["failed_comparisons"], "equals": []}
    ]}
  }
}
```

Use the **actual field paths in the reviewed candidate result JSON**, not the
illustrative field names above. Every domain needs nonempty evidence and checks;
each checked record must appear in that domain's evidence. Paths are key/index
arrays into JSON. Assertions permit only strict `true`, integer `0`, or `[]` and
are recomputed from the referenced results. All JSON evidence and its directly
inventoried raw source/input/output records are rehashed, including legacy
`inputs`/`artifacts` path-to-SHA maps (relative artifact paths resolve beside the
result JSON). Every **checked result** must itself contain the exact candidate
build record, not merely borrow the outer gate's build declaration. A new sealed
scoped review result may bind the candidate build, actual independent results,
raw inputs/outputs, and reviewed pass conditions where an older result schema
lacks such a binding; never mutate old results to add provenance. The review must retain
the baseline's explicit scientific limitations and distinguish them from a new
unexplained discrepancy; universal scientific correctness is not asserted.

Engineering uses:

```json
{
  "kind": "corrected_candidate_engineering_gate",
  "complete": true,
  "passed": true,
  "baseline_build": "RECORD",
  "candidate_build": "RECORD",
  "sources": ["EXACT candidate build sources + python_sources list"],
  "commands": [
    {"name": "rust-tests", "argv": ["cargo", "test", "--locked", "--all-features"],
     "returncode": 0, "log": "RECORD"}
  ]
}
```

Required distinct command names are `rust-tests`, `python-tests`, `clippy`,
`rust-format`, `rustdoc`, `ruff-check`, and `ruff-format`. Every command must have
an argv, zero return code, and a hashed log. Additional named commands are
allowed but must also succeed. Record the real command and complete source
binding; do not create a pass wrapper around tests of a different source tree.

Finally, the complete prediction comparison uses:

```json
{
  "kind": "corrected_prediction_bitstream_comparison",
  "complete": true,
  "exact_equivalence": true,
  "baseline_build": "RECORD",
  "candidate_build": "RECORD",
  "encoding": "ieee754-f32-le",
  "records": 1000000,
  "data": "RECORD of frozen predict --data input",
  "weights": "RECORD of frozen predict --weights input",
  "threads": {
    "1": {"baseline": "RECORD of raw stream", "candidate": "RECORD of raw stream"},
    "2": {"baseline": "RECORD of raw stream", "candidate": "RECORD of raw stream"},
    "4": {"baseline": "RECORD of raw stream", "candidate": "RECORD of raw stream"},
    "6": {"baseline": "RECORD of raw stream", "candidate": "RECORD of raw stream"}
  },
  "exporter_builds": {"baseline": "RECORD of exporter builds.json", "candidate": "RECORD of exporter builds.json"}
}
```

Each stream must be exactly 4,000,000 bytes, with equal baseline/candidate SHA at
the same thread count. Generate the streams using the existing
`examples/profile_predictions.rs` from each corresponding frozen source/build,
and retain capture/exporter build evidence. `exporter_builds` uses the existing
`corrected_frozen_diagnostic_builds` receipt schema: complete, successful native
build command, native executable identical to its admitted build, all source
bindings, and four native export captures. The verifier checks the compiled
source files against the corresponding admitted build, requires identical
exporter adapter bytes, and checks each capture's executable/input/output argv,
thread affinity, zero exit status, and vector record. Additional diagnostic
builds/exports may be retained but never count as native timing evidence.
That exporter writes native-endian
bytes; this host is little-endian. Never label a big-endian native stream as LE.

Focused harness tests (mock measurements only):

```sh
.venv/bin/python -m unittest discover -s tests -p test_corrected_pair.py -v
.venv/bin/ruff check scripts/benchmark_corrected_pair.py tests/test_corrected_pair.py
```
