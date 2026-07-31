# Reproducing StericX

This document walks a reviewer from a clean `git clone` to reproduced results.
Every study is deterministic on a given toolchain; the checked-in numbers below
are the targets to match. If a number differs, the difference is itself the
finding — no gate is hidden.

## Prerequisites

- Linux x86-64 (the buried-volume and benchmark numbers are Linux-measured).
- `git`, `curl`, and a C-capable environment. Everything else is installed by
  the bootstrap script: a current stable Rust toolchain and
  [`uv`](https://docs.astral.sh/uv/) for the pinned Python environment.

## One-command setup

```bash
git clone https://github.com/AndrejRumenovski/StericX.git
cd StericX
./bootstrap.sh                 # uv + Rust + build + binary smoke test
./bootstrap.sh --with-quantum  # ALSO fetch CREST 2.12 / xTB 6.4.0 (Study 003)
```

`bootstrap.sh` is idempotent. It builds a **portable** binary by default; pass
`--native` for a `target-cpu=native` build (faster, but it will not run on a
different CPU — the cause of an `Illegal instruction` crash if a binary is
copied between machines).

The quantum toolchain and Python virtualenv live under the git-ignored
`.stericx/` and `.venv/` directories and are never committed.

## Chemical fidelity — Sterimol vs. morfeus

```bash
uv run --extra science python scripts/validate_stericx.py
```

| Parameter | Expected R² | RMSE |
|---|---:|---:|
| L  | 1.000000 | 0.000000 Å |
| B1 | 0.999959 | 0.010539 Å |
| B5 | 1.000000 | 0.000001 Å |

## Study 001 — Ni-hDA reproduction

```bash
uv run --extra science python studies/study_001_ni_hda.py --offline
```

Published Kraken-descriptor OLS: training R² = 0.8193, LOO Q² = 0.7521,
LOO RMSE = 0.3430 kcal/mol, historical-blind MAE = 0.3730 kcal/mol. The native
compact-descriptor ablation is intentionally weaker (LOO Q² ≈ 0.002) and is
retained, not hidden.

## Study 002 — Coordination-aware buried volume

```bash
uv run --extra science python studies/study_002_buried_volume.py
```

| Validation target | Expected |
|---|---|
| Native voxel geometry vs Morfeus | R² = 1.000000 (worst mean rel. error 0.000008 %) |
| Native ensemble descriptor vs Kraken | R² = 0.8626, RMSE 2.8740 Å³ (below preregistered target) |
| Native buried-volume Ni-hDA model | train R² = 0.7693, LOO Q² = 0.6549 |
| Historical blind ligand 723 | abs. error 0.1529 kcal/mol |

## Study 003 — Quantum geometry (needs `--with-quantum`)

Phase A isolates the xTB LMO coordination centre on the existing RDKit/MMFF
conformers:

```bash
uv run --extra science python scripts/prepare_quantum_data.py --mode lmo
uv run --extra science python studies/study_003_quantum_geometry.py --no-build
```

Expected phase-A: descriptor R² = 0.8517 vs Kraken (a regression from Study
002, isolating conformational sampling as the next variable).

The full eleven-ligand production CREST ensemble (hours of CREST/GFN2-xTB;
cache-resumable):

```bash
uv run --extra science python scripts/prepare_quantum_data.py \
  --mode crest --threads 6 --lmo-workers 6 \
  --output-csv data/reactions_crest.csv \
  --provenance data/quantum/crest_provenance.json
uv run --extra science python studies/study_003_quantum_geometry.py --no-build \
  --reactions-csv data/reactions_crest.csv \
  --quantum-provenance data/quantum/crest_provenance.json
```

| Validation target | Expected | Gate |
|---|---|---|
| Native descriptor vs Kraken | R² = 0.9254, RMSE 2.8354 Å³ (+0.0627 over Study 002) | improves over Study 002 |
| Native descriptor R² > 0.99 | R² = 0.9254 | fail |
| Fixed-feature Ni-hDA model | LOO Q² = 0.5941 | below published target |
| Historical replay ligand 723 | abs. error 0.1107 kcal/mol | historical only |

CREST conformer counts per ligand are cached under `.stericx/cache/`; a warm
cache replays the full ensemble in seconds.

## Benchmarks (optional)

```bash
./benchmark_linux.sh
```

Rebuilds with `target-cpu=native`, generates a one-million-record mapped
matrix, and records `/usr/bin/time -v` metrics to
`docs/benchmark_results.json`. These numbers are machine-specific by design.
