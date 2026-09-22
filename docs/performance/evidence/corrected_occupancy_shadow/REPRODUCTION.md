# Corrected occupancy shadow evidence

This is a compact, lossless publication of the diagnostic study under
`.stericx/profiling/scientifically_validated_optimization/corrected_occupancy_shadow_v1`.
It supports the [C3 proposal](../../candidates/C3_ROW_REJECTION_PROPOSAL.md).
It contains no production optimization or measured C3 speedup.

The original report is retained byte for byte as [RAW_REPORT.md](RAW_REPORT.md),
including its cosmetic `and33.600897%` spacing typo. No historical report or receipt
was edited during publication. [ASSEMBLY.json](ASSEMBLY.json) maps each published
file to its original absolute path, byte count and SHA-256, and records the stored
encoding/hash. Controllers use `.py.txt` so their exact executed bytes remain
readable without entering Python lint discovery. Rename/copy them only when running
an explicitly new replay; do not rewrite the retained originals.

## What is retained

- [summary.json](summary.json), [completion_binding.json](completion_binding.json),
  [source_manifest.json](source_manifest.json), [build_manifest.json](build_manifest.json),
  pre-edit identities/hypothesis, build/test logs and environment records.
- Exact executed `run_corrected_probe.py`, `legacy_run_probe.py`,
  `finalize_evidence.py` and profiling helpers under `controllers/`, plus the exact
  Rust append, original/instrumented occupancy functions and patches.
- `instrumented_source.tar.gz`: all 128 private source files with unchanged contents.
  Reversing `probe.patch` recovers every original baseline source hash.
- All 30,168 raw frame records under `captures/*/counts.jsonl.xz`, not samples;
  capture commands/results and complete stdout/stderr. Equal native/probe streams
  are stored once, with both original identities retained in their receipts.
- Complete baseline/workload manifests compressed under `dependencies/`, plus the
  corrected native summary used for fingerprint comparison. These manifests also
  identify the frozen input files; the input dataset itself is not duplicated here.
- [static_operation_model.json](static_operation_model.json), which labels
  source-derived transforms, radius lookups and vector construction sites separately
  from measured loop counters.

[RAW_ARCHIVE_INVENTORY.json](RAW_ARCHIVE_INVENTORY.json) inventories the original
local archive, including retained native/probe binaries. It excludes disposable
Cargo `target/` and Python `__pycache__/` directories. Binaries are not duplicated in
Git; their source/build and pre/post execution hashes remain in the receipts. The
local original `bin/probe` and `bin/baseline_native` are verified copies of the
executed bytes. Preserved absolute paths identify the historical capture environment;
they are not portable instructions to overwrite another checkout.

## Verify without replaying science

From the repository root, with Python 3 and the standard `patch` utility:

```sh
python3 docs/performance/evidence/corrected_occupancy_shadow/verify_package.py.txt
```

This reads the retained streams, checks compressed and decompressed identities,
verifies receipt links, reverses the source patch in a temporary directory, and
checks all 128 recovered baseline files. It recomputes every measured counter total,
row histogram, complete file/orientation coverage and the ordered 452,520-value bit
fingerprint. It also checks full stdout/stderr identities and native workload
fingerprints. It performs no Rust build, scientific replay or native benchmark.
The captured output is [verification_result.json](verification_result.json).
`PACKAGE_MANIFEST.json` inventories the final publication, including this new
verifier and its result; the verifier is separate from the original executed probe.

The original runtime assertions compared every would-skip predicate and all 15
per-frame output bits. Offline verification checks their retained evidence and
coverage; it does not independently rerun those floating-point operations.

## Rebuild the diagnostic source

Run builds and replays only outside native benchmark windows. Use a fresh directory;
never build inside a retained archive or overwrite an original capture. For example:

```sh
shadow_replay=$(mktemp -d)
tar -xzf docs/performance/evidence/corrected_occupancy_shadow/instrumented_source.tar.gz -C "$shadow_replay"
cargo build --offline --locked --release --bin stericx --manifest-path "$shadow_replay/source/Cargo.toml" --target-dir "$shadow_replay/target"
cargo test --offline --locked --release --lib row_range --manifest-path "$shadow_replay/source/Cargo.toml" --target-dir "$shadow_replay/target"
```

The exact original compiler/Cargo versions and commands are in `build_manifest.json`;
allowlisted build overrides are in `build_environment.json`. Offline Cargo requires
the pinned dependency cache. Rebuilding with a different toolchain or absolute source
location can produce a different executable hash; retain a new build receipt and
compare scientific outputs rather than asserting historical binary identity.

The source archive contains the additive shadow hooks. Its ordinary CLI still emits
original scientific results, but its runtime includes extra replay/counter work and
must never be used for native performance timings. Setting
`STERICX_ROW_RANGE_PROBE_PATH` to a new sidecar path records the counters; the probe
refuses to overwrite that sidecar. The test filter executes both focused shadow
boundary/fallback test groups.

## Replay the exact retained workloads

Restore/verify the frozen input files named by `dependencies/workloads.json.xz`
using their original byte hashes. All selected XYZ inputs must be present; the
10,000-file workload contains 56 distinct byte contents repeated under distinct
paths. The published command receipts retain the exact input sequence and flags.
The original run used `RAYON_NUM_THREADS=1`, `LC_ALL=C`, `TZ=UTC`,
`OMP_NUM_THREADS=1`, and no `STERICX_PROFILE_PATH`. Copy each command's `argv`, replacing
only its executable with the newly built diagnostic binary. Set a fresh
`STERICX_ROW_RANGE_PROBE_PATH`; retain argv/environment, binary hash, exit status,
complete stdout/stderr and the sidecar before checking results.

Use argument arrays, not shell expansion of the 10,000 paths. The 56-file command
is readable at `captures/conformers_56/probe/command.json`; the larger command is
`captures/descriptors_10000/probe/command.json.gz`. A short Python loader is:

```python
import gzip, json
from pathlib import Path

package = Path("docs/performance/evidence/corrected_occupancy_shadow")
command = json.loads(
    gzip.decompress(
        (package / "captures/descriptors_10000/probe/command.json.gz").read_bytes()
    )
)
argv = command["argv"]
# Set argv[0] to the newly built probe; execute this argument array in a new
# capture directory after verifying each input against the workload manifest.
```

At the original input paths, successful full stdout fingerprints must be:

| Workload | stdout SHA-256 |
| --- | --- |
| `conformers_56` | `8098e16239bd7740ae4061868d83b45cd0f3989fa950bd869fdef72a16ea606b` |
| `descriptors_10000` | `df2815f9b1411e7d77a62de7bd1074747ff5cface781f62857e278553da37037` |

Relocating inputs changes emitted filename fields and capture identities. Record
that as a new controlled replay with an explicit path mapping; do not silently
claim the historical raw-byte fingerprint. Recompute coverage and every counter,
retain all 15-bit fields and original-vs-shadow assertions, and preserve failures.
The exact original controller is archived for review, but its identity checks bind
the historical frozen roots and its phase outputs refuse overwrites; running it
against this publication directory is not a portable shortcut.

A reproducible shadow result establishes operation counts and differential output
identity. C3 still requires independent scientific/reference/engineering gates and
quiet paired native measurements against its newly accepted C1/C2 baseline.
