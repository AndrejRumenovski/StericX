# Corrected baseline counter evidence

This compact package preserves 84 exact gzip copies (1,223,963 compressed bytes)
from the completed corrected-baseline Cachegrind/DHAT investigation. The relative
archive paths, compressed hashes, and original uncompressed hashes and paths are
in [manifest.json](manifest.json), with its SHA-256 sidecar. Gzip timestamps are
fixed to zero; decompression restores the original bytes.

- `runs/`: all six raw counter streams, per-run commands/status, profiler logs,
  result summaries, and Cachegrind function annotations.
- `source/`: the frozen corrected Rust source and Cargo configuration used for
  the baseline symbol executable.
- `helpers/`: the unchanged profiling harness and native resource launcher, plus
  the publication helper.
- `receipts/`: the counter/report manifests, model assumptions, hardware/tool
  capability receipts, diagnostic/symbol build receipt, report helper and parser
  tests, and compact counter data.

The investigation covers the full unchanged 56-conformer, 1,000-row screen, and
10,000-file descriptor workloads, one thread on CPU2. All six scientific output
fingerprints matched native. The [counter report](../../CORRECTED_COUNTERS.md)
distinguishes modeled events and intercepted heap traffic from native timings.

This package is **not a standalone scientific admission**. Large input, native
oracle, scientific stdout/stderr, and full controller inventories remain in the
hash-bound original `.stericx` evidence paths. The complete controller receipts
repeat the full input inventories; their hashes remain in the counter manifest
rather than duplicating those inventories six times here. Tool binaries and
preload libraries are identified by the capability receipt, not redistributed.
The original complete evidence and sealed profile reports remain unchanged.

For example, `gzip -dc source/Cargo.toml.gz` recovers the frozen Cargo manifest.
Verify compressed files against `archive.sha256`, then verify decompressed bytes
against `original.sha256` in the publication manifest.
