# Frozen baseline

Status: source and dataset freeze, fresh baseline builds, tests, 128,021 scientific observations and pre-optimization profiling complete. See [PROFILING_REPORT.md](PROFILING_REPORT.md) and [BASELINE_TABLES.md](BASELINE_TABLES.md). The baseline archive remains immutable.

Commit: `c03630cb3595eb85b5d05bdbecd3f10305b818bf`. Existing untracked `docs/media/` is untouched. The new campaign directories account for the additional status entries in the manifest.

Source snapshot: `.stericx/ultimate_engine/baseline/source.tar.gz`; SHA-256 `9d8f3c88dd0dda5cfb67b90f4edbfea298131ec8a9b3f6364419e48f7d54ac4d`. Source, dataset and historical report hashes are in [baseline_manifest.json](baseline_manifest.json). All tracked repository paths and sizes are retained in the local repository inventory.

Rust 1.97.0, LLVM 22.1.6; Linux x86_64 AMD Ryzen 5 5600G, six cores/twelve threads. Existing release profile uses thin LTO and one codegen unit. No ambient RUSTFLAGS or target override. The campaign will distinguish portable default builds from native-target experiments.

The 10,774.5 files/s historical result is context only. Fresh measurements are required for all new performance claims. Historical failures and reference conventions remain preserved.

Fresh checks: 348 Rust tests passed; Python unittest discovery ran 156 tests (143 passed, 13 explicit external-reference skips). The skips are pending separately executed pinned reference lanes, not passes. All 13 scientific observation lanes completed and retain raw responses, build hashes and input hashes in `.stericx/ultimate_engine/baseline_observations`.

The freshly measured exploratory six-thread 10,000-file throughput was 10,264.5 files/s; single-thread 1,977.2 files/s. These desktop measurements are candidate-selection evidence. Final quiet-window A/B measurements are still required. Profiling produced 175 CLI observations, 71 kernel/batch observations, Callgrind/Cachegrind/DHAT/Memcheck output, assembly, Python layers and analytical volume convergence evidence.
