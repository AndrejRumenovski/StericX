# StericX

**A native-Rust engine for physical-organic molecular featurization and reaction-selectivity inference.**

[![CI](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license)

StericX computes 3D **Sterimol** and coordination-aware **buried-volume**
descriptors directly from Cartesian coordinates, stores reaction observations in
cache-aligned binary matrices, evaluates interpretable linear models with
parallel SIMD kernels, and converts predicted energy differences into product
distributions through the Eyring equation — as a single, portable binary.

It exists because the standard steric-descriptor toolchain is Python-based and
research-specific. StericX is an **independent, from-scratch reimplementation**
whose results are validated two ways: numerically against `morfeus-fsu`
(\(R^2 \geq 0.9999\)) and against Kraken's *own published values* across the
**full 1,541-ligand library**. Every number below is a generated measurement,
and every documented failure is kept in view.

📄 Manuscript-style write-up: [`docs/REPRODUCTION_REPORT.md`](docs/REPRODUCTION_REPORT.md)
· 🖼️ One-page visual overview: [`docs/results.html`](docs/results.html)
· 🔁 Clone-to-results walkthrough: [`REPRODUCE.md`](REPRODUCE.md)

---

## Why StericX?

Most steric-descriptor tooling (Kraken, morfeus) is Python and tuned to a
specific research workflow. StericX is different in four defensible ways:

- **Native Rust, single binary.** No Python runtime, no atom-index bookkeeping —
  point the tool at an `.xyz`/`.sdf`/`.mol` file and it auto-detects the donor
  from the geometry.
- **Validated against two independent references.** Numerical fidelity vs
  `morfeus-fsu` on identical geometries, *and* reproduction of Kraken's published
  DFT descriptor values at full library scale (1,541 ligands / 31,611
  conformers).
- **Failures stay visible.** Where a cheaper conformer pipeline or a compact
  descriptor underperforms, the weaker result is kept and explained rather than
  hidden — the project's entire purpose is an honest, defensible reproduction.
- **Systems-engineered core.** Cache-aligned 64-byte records, zero-copy
  memory-mapped reads, and an explicit AVX2 inference kernel.

StericX is an independent reproduction — not affiliated with, endorsed by, or
produced by the Sigman or Reisman groups. See [References](#references-and-attribution).

---

## Research Highlights

- **Descriptor validation** — Sterimol \(L\), \(B_1\), \(B_5\) reproduce
  `morfeus-fsu` to \(R^2 \geq 0.9999\); the buried-volume voxel kernel matches
  Morfeus to \(R^2 = 1.000000\) on identical geometries (worst mean relative
  error \(8 \times 10^{-6}\%\)).
- **Dataset scale** — Kraken's published `vbur_max_delta_qvbur_min` reproduced
  across **1,541 chemically diverse ligands / 31,611 DFT conformers**.
- **Reproduction accuracy** — \(R^2 = 0.9986\) on the 11-ligand DFT set;
  \(R^2 = 0.9852\) (median absolute error **0.11 Å³**) at full scale.
- **Descriptor families** — whole buried-volume family (8 descriptors) mean
  \(R^2 = 0.9925\); Sterimol on the coordination axis mean \(R^2 = 0.9887\);
  pyramidalization (`pyr_P`, `pyr_alpha`) mean \(R^2 = 0.99998\).
- **Performance** — 62k molecules/s Sterimol extraction, 317.8 M evals/s SIMD
  inference over a memory-mapped million-record matrix.
- **Reproducibility** — one-command bootstrap, checksum-pinned CREST 2.12 / xTB
  6.4.0, content-addressed caches, frozen prediction hashes, green CI, dual
  MIT/Apache licensing.
- **Scientific studies** — eight studies (Ni-hDA model, buried-volume fidelity,
  quantum geometry, Kraken DFT reproduction, pyramidalization, residual
  localization, an independent cross-coupling reaction model, and a head-to-head
  speed benchmark vs morfeus), each with passed *and* failed gates recorded.
- **Transparency at scale** — a genuine kernel bug surfaced by the full library
  was fixed **without discarding a single ligand**, and the remaining residual is
  fully characterized by phosphine class.

---

## Architecture

The pipeline runs from raw geometry to a validated selectivity prediction. The
descriptor engine feeds both the modeling path and the two independent
assurance paths (validation and benchmarking).

```mermaid
flowchart LR
    A["Geometry<br/>.xyz · .sdf · .mol<br/>RDKit · CREST · xTB · DFT"]
    B["Descriptor Engine<br/>Sterimol L/B1/B5<br/>buried-volume family"]
    C["Binary Storage<br/>64-byte .sigpack<br/>memory-mapped, zero-copy"]
    D["Model<br/>BIC forward selection<br/>AVX2 SIMD inference"]
    E["Kinetics<br/>Eyring ΔΔG‡ → R:S, ee"]
    F["Validation<br/>vs morfeus · Kraken"]
    G["Benchmarking<br/>throughput · RSS · faults"]

    A --> B --> C --> D --> E
    B --> F
    C --> G
```

---

## Quick Start

For a fresh checkout, one script installs the toolchains, builds the binary, and
smoke-tests it:

```bash
./bootstrap.sh                 # uv + Rust + build + smoke test
./bootstrap.sh --with-quantum  # also fetch CREST 2.12 / xTB 6.4.0 (Study 003)
```

Then featurize any phosphine geometry — the phosphorus donor and its
substituents are detected from the structure, so there is no CSV to build and no
atom indices to look up:

```bash
./target/release/stericx descriptors data/xyz/SIG-NIHDA-401_9d42bff1.xyz
```

```text
data/xyz/SIG-NIHDA-401_9d42bff1.xyz
  donor          P (atom 14)
  substituents   C, C, C
  conformers     1
  Sterimol      L 7.76   B1 1.73   B5 8.75   Å
  buried volume  Vbur 24.3%   (43.6 Å³)
                 qvbur_min 9.73   qvbur_max 13.53   max_delta_qvbur 3.80   Å³
  pyramidalization pyr_P 0.932   pyr_alpha 17.10°
```

The full command set is in the [CLI Reference](#cli-reference); manual builds and
the reproducible Python environment are covered below.

### Manual build

```bash
cargo build --release
cargo test --all-targets
```

For a machine-native optimized binary:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### Python tooling

The reproduction studies use Python, installed reproducibly with
[uv](https://docs.astral.sh/uv/):

```bash
uv sync --extra science
```

Then prefix Python workflows with `uv run` (for example `uv run scripts/prepare_data.py`).
Study 003 optionally uses the exact public Kraken-era quantum toolchain; the
installer verifies the release-archive checksums and keeps the executables and
calculation cache under the ignored `.stericx/` directory:

```bash
./install_quantum_tools.sh
```

---

## Features

### Scientific capabilities

- **Compact physical-organic features** — the model combines \(L\), \(B_1\),
  \(B_5\), donor NBO charge, IR frequency, and interpretable steric–electronic
  interaction terms.
- **Coordination-aware sterics** — a deterministic metal-centred voxel engine
  calculates total, quadrant, octant, near/far, and conformer-ensemble buried
  volumes using the public Kraken/Morfeus convention.
- **Pyramidalization** — Radhakrishnan `pyr_P` and the mean out-of-plane angle
  `pyr_alpha` from the donor's three bond vectors, matching Kraken/morfeus.
- **Transition-state kinetics** — Eyring calculations convert
  \(\Delta\Delta G^\ddagger\) into rate constants, R:S distributions, and
  enantiomeric excess.

### Systems engineering

- **Zero-copy mapped reads** — reaction observations use an exact 64-byte
  `#[repr(C, align(64))]` layout. `bytemuck` provides POD byte casting and
  `memmap2` exposes `.sigpack` files as borrowed record slices without
  per-record deserialization.
- **Hardware SIMD inference** — the eight-feature RegressX dot product uses an
  explicit eight-lane AVX2 kernel on supported x86/x86-64 processors, with a
  portable unrolled fallback.
- **Parallel work-stealing** — `rayon` distributes independent inference rows
  across CPU workers. The prediction hot loop has no shared mutable model state
  or per-record locks.

---

## Validation

Two independent references anchor every descriptor claim: `morfeus-fsu` for
numerical fidelity on identical geometries, and Kraken's own *published* values
for chemical accuracy at library scale.

### Numerical fidelity vs `morfeus-fsu`

[`scripts/validate_stericx.py`](scripts/validate_stericx.py) evaluates every structure in
`data/xyz/` with both StericX and `morfeus-fsu`, using identical attachment
indices and atomic radii, and renders 400-DPI correlation plots.

```bash
python scripts/validate_stericx.py
```

The current checked-in validation run contains 11 structures.

| Parameter | \(R^2\) | Slope | Intercept (Å) | RMSE (Å) | Status |
|---|---:|---:|---:|---:|---|
| \(L\) | 1.000000 | 1.000000 | 0.000001 | 0.000000 | Validated |
| \(B_1\) | 0.999959 | 1.000249 | 0.008294 | 0.010539 | Validated |
| \(B_5\) | 1.000000 | 1.000000 | -0.000001 | 0.000001 | Validated |

The comparison uses Morfeus's uncorrected \(L\) value because both tools then
report the raw geometric extent \(\max_i(z_i + r_i)\). StericX scans \(B_1\) at
one-degree intervals while Morfeus uses a denser angular search; the remaining
discretization difference is 0.0105 Å RMSE. Separately, the buried-volume voxel
kernel matches Morfeus to \(R^2 = 1.000000\) (worst mean relative error
\(0.000008\%\)) on identical structures — the geometry engine is exact; the open
question is always the input geometry, which the studies below isolate.

| ![Sterimol L correlation](docs/sterimol_l_corr.png) | ![Sterimol B1 correlation](docs/sterimol_b1_corr.png) | ![Sterimol B5 correlation](docs/sterimol_b5_corr.png) |
|:---:|:---:|:---:|
| \(L\) | \(B_1\) | \(B_5\) |

### Reproducing Kraken's published values (full library)

Running the unchanged StericX kernel on Kraken's own DFT geometries at Kraken's
documented 2.28 Å convention reproduces the published descriptors across the
whole public library. Full derivation is in
[Scientific Studies](#scientific-studies) (Study 004).

| Target vs published Kraken (DFT geometries) | Scale | Result |
|---|---|---:|
| `vbur_max_delta_qvbur_min` | 11 ligands | \(R^2 = 0.9986\) |
| `vbur_max_delta_qvbur_min` | 1,541 ligands / 31,611 conformers | \(R^2 = 0.9852\) · medAE 0.11 Å³ |
| Buried-volume family (8 descriptors) | 1,541 ligands | mean \(R^2 = 0.9925\) |
| Sterimol \(L\)/\(B_1\)/\(B_5\) (coordination axis) | 1,541 ligands | mean \(R^2 = 0.9887\) |
| Pyramidalization `pyr_P`, `pyr_alpha` (min/max) | 1,541 ligands | mean \(R^2 = 0.99998\) |

---

## Benchmarks

[`benchmark_linux.sh`](benchmark_linux.sh) builds StericX with
`-C target-cpu=native`, constructs a 10,000-structure workload, generates a
64 MB one-million-record mapped matrix, and collects Linux process metrics with
`/usr/bin/time -v`.

```bash
./benchmark_linux.sh
```

| Workload | Items | Phase latency | Operations/sec | Throughput | Peak RAM |
|---|---:|---:|---:|---:|---:|
| Sterimol extraction | 10,000 molecules | 0.160964 s | 62,125.7 ops/s | 62,125.7 mol/s | 7.75 MB |
| Binary `.sigpack` packing | 10,000 records | 0.000176 s | 56,818,181.8 ops/s | 56.82 M records/s | 7.75 MB |
| RegressX SIMD prediction | 1,000,000 records | 0.003147 s | 317,762,948.8 ops/s | 317.76 M evals/s | 58.94 MB |

The extraction and packing latencies come from separate internal phase timers
inside one `parse` process; consequently, their GNU `time` CPU, RSS, and page
fault measurements are shared. The mapped prediction workload recorded zero major
page faults and 3,176 minor page faults. Machine-readable results (including CPU
utilization and page faults) are in
[`docs/benchmark_results.json`](docs/benchmark_results.json). Run the benchmark
on any target workstation to replace these with results from that system.

---

## Scientific Studies

Four preregistered studies build from a small published reaction family up to the
full Kraken library. Each writes a complete model card, frozen prediction hashes,
raw comparisons, plots, and machine-readable results under `docs/study_00N/`.
**Passed and failed gates are both retained.**

### Dataset preparation

[`scripts/prepare_data.py`](scripts/prepare_data.py) downloads the public Sigman Group
Ni-catalyzed homo-Diels–Alder/Kraken table when available. It preserves the
complete 1,566-ligand source table and exact-download SHA-256 provenance,
identifies the ten published training ligands and historical holdout, and
generates deterministic ETKDGv3/MMFF94 conformer ensembles for measured rows
(optimized, energy-windowed, and Boltzmann-weighted before export).

```bash
python scripts/prepare_data.py            # or --offline for the embedded 100-reaction benchmark
```

```text
data/official/ni_hda_kraken.csv   Exact complete public source
data/official/provenance.json     Source URL, hash, counts, and split IDs
data/conformers/                   56 retained conformer geometries
data/reactions_raw.csv             10 training + 1 historical-blind reaction
```

### Study 001 — Ni-hDA enantioselectivity model

[`studies/study_001_ni_hda.py`](studies/study_001_ni_hda.py) reproduces the published ten-ligand
enantioselectivity relationship using the preregistered Kraken descriptor
`vbur_max_delta_qvbur_min`, writing and hashing the ligand-723 prediction before
revealing the experimental target.

```bash
python studies/study_001_ni_hda.py --offline
```

| Model | Training \(R^2\) | LOO \(Q^2\) | LOO RMSE | Historical-blind MAE |
|---|---:|---:|---:|---:|
| Published Kraken descriptor OLS | 0.8193 | 0.7521 | 0.3430 kcal/mol | 0.3730 kcal/mol |
| Curated descriptors, nested ridge | — | 0.6652 | 0.3986 kcal/mol | — |
| Curated descriptors, nested LASSO | — | 0.6558 | 0.4042 kcal/mol | — |
| Native StericX ensemble descriptors | 0.3625 | 0.0020 | 0.6882 kcal/mol | 0.6978 kcal/mol |

The native result is an **intentional descriptor ablation**: it shows that
compact Sterimol/NBO features do not replace the published coordination-aware
buried-volume descriptor for this small reaction family. Full model card,
limitations, frozen prediction hashes, residuals, correlations, and VIF table:
[`docs/study_001/STUDY_001.md`](docs/study_001/STUDY_001.md).

![Ni-hDA historical holdout parity](docs/study_001/ni_hda_parity.png)

### Study 002 — Coordination-aware buried-volume fidelity

[`studies/study_002_buried_volume.py`](studies/study_002_buried_volume.py) freezes the descriptor
definition from the public Kraken implementation, calculates a trusted Morfeus
reference on every retained conformer, executes the native Rust implementation,
validates `.sigpack` v2 aggregation, and reruns the unchanged Ni-hDA
train/holdout partition.

```bash
uv run --extra science python studies/study_002_buried_volume.py
```

| Validation target | Result | Status |
|---|---:|---|
| Native voxel geometry vs Morfeus | \(R^2 = 1.000000\), worst mean relative error \(0.000008\%\) | Pass |
| Native ensemble descriptor vs official Kraken | \(R^2 = 0.8626\), RMSE \(2.8740\) Å³ | Below preregistered target |
| Native buried-volume Ni-hDA model | Train \(R^2 = 0.7693\), LOO \(Q^2 = 0.6549\) | Below published LOO target |
| Historical blind ligand 723 | Absolute error \(0.1529\) kcal/mol | Pass |

The distinction matters: the Rust geometry kernel matches Morfeus to numerical
precision on identical structures, while inexpensive RDKit/MMFF conformers and an
inferred lone-pair direction do not fully reproduce Kraken's CREST/xTB/DFT and
localized-orbital workflow. Failed gates remain visible. Full specification and
model card: [`docs/study_002/STUDY_002.md`](docs/study_002/STUDY_002.md).

| ![Buried-volume anisotropy parity](docs/study_002/max_delta_qvbur_parity.png) | ![Native buried-volume Ni-hDA model](docs/study_002/ni_hda_native_buried_volume_parity.png) |
|:---:|:---:|

### Study 003 — Quantum geometry & prospective validation

Study 003 adds a checksum-pinned CREST 2.12/GFN2-xTB 6.4.0 backend, immutable
content-addressed caches, executable hashes, full command logs, exact
Kraken-compatible phosphorus LMO selection, and explicit per-conformer
coordination centres.

**Phase A** keeps the 56 existing ETKDGv3/MMFF94 conformers and isolates the
effect of replacing the inferred centre with the xTB LMO centre:

```bash
uv run --extra science python scripts/prepare_quantum_data.py --mode lmo
uv run --extra science python studies/study_003_quantum_geometry.py --no-build
```

| Validation target | Result | Gate |
|---|---:|---|
| All conformers have explicit xTB LMO centres | 56 / 56 | Pass |
| Native descriptor vs official Kraken | \(R^2 = 0.8517\), RMSE \(3.0317\) Å³ | Fail |
| Fixed-feature Ni-hDA model | LOO \(Q^2 = 0.6314\), RMSE \(0.4182\) kcal/mol | Fail |
| Historical replay of ligand 723 | Absolute error \(0.1902\) kcal/mol | Historical only |
| Frozen target-free prospective deck | 10 candidates, targets unaccessed | Pass |

The negative phase-A result is retained: changing the centre alone *reduced*
agreement from Study 002's \(R^2 = 0.8626\), isolating conformational sampling as
the next controlled variable. **Phase B** replaces the ensemble with 322 freshly
sampled CREST conformers using the same CREST profile Kraken recorded:

```bash
uv run --extra science python scripts/prepare_quantum_data.py \
  --mode crest --threads 4 --lmo-workers 4 \
  --output-csv data/reactions_crest.csv \
  --provenance data/quantum/crest_provenance.json
```

CREST ensembles and xTB LMO calculations occupy independent content-addressed
cache stages; CREST is committed durably before LMO work begins, LMO jobs run in
a CPU-bounded worker pool and checkpoint after every conformer, and interrupted
jobs resume from cache while PID/start-time-aware locks reclaim dead owners
without stealing live work.

| Validation target | Result | Gate |
|---|---:|---|
| Native descriptor vs official Kraken | \(R^2 = 0.9254\), RMSE \(2.8354\) Å³ (\(+0.0627\) over Study 002) | Improves over Study 002 |
| Native descriptor \(R^2 > 0.99\) | \(R^2 = 0.9254\) | Fail |
| Fixed-feature Ni-hDA model | LOO \(Q^2 = 0.5941\), RMSE \(0.4389\) kcal/mol | Below published target |
| Historical replay of ligand 723 | Absolute error \(0.1107\) kcal/mol | Historical only |
| Frozen target-free prospective deck | 10 candidates, targets unaccessed | Pass |

Better conformers raised descriptor agreement to the best of the three geometry
pipelines and lowered the ligand-723 error to 0.1107 kcal/mol — but the small
ten-ligand model's LOO \(Q^2\) still fell to 0.5941, so better descriptor
fidelity did **not** translate into better held-out kinetic prediction on this
family. Both outcomes remain visible. Full results:
[`docs/study_003/STUDY_003.md`](docs/study_003/STUDY_003.md).

**Pre-registered, falsifiable prediction.** The frozen ten-ligand deck is
elevated to a full pre-registration —
[`docs/study_003/PREREGISTRATION.md`](docs/study_003/PREREGISTRATION.md) — that
commits (by the deck's SHA-256, dated, before any measurement) to a 95%
prediction interval and applicability-domain judgement for every candidate, an
exact experimental protocol (asymmetric Ni-hDA; ee → ΔΔG‡ via `2RT·atanh(ee)`),
and a pre-registered falsification rule: **corroborated iff** measured ΔΔG‡ lands
within the 95% interval for ≥6 of 8 primary candidates **and** the
predicted-vs-measured rank correlation is positive, **else falsified**. It states
its own honest limit — for 7 of 10 candidates the 95% ee interval spans both
enantiomers, so those near-racemate predictions are not falsifiable on ee, and
the test runs on ΔΔG‡. Unmeasured by design: the rigor of the open prediction is
the point. Regenerate with
[`scripts/preregister_prediction.py`](scripts/preregister_prediction.py).

![Pre-registered Ni-hDA forecast](docs/study_003/prospective_prediction_forecast.png)

### Study 004 — Reproducing Kraken's published descriptors on DFT geometries

Studies 002 and 003 left the native descriptor below the official values
(\(R^2 = 0.8626\) and \(0.9254\)) without isolating whether the shortfall came
from the StericX kernel or the cheaper conformer geometries.
[`studies/study_004_reproduction.py`](studies/study_004_reproduction.py) settles it
by downloading Kraken's own DFT-optimized geometries (PBE/6-31+G(d,p), GD3BJ)
from the public MolSSI descriptor-library REST API and running the unchanged
kernel on them.

```bash
uv run --extra science python studies/study_004_reproduction.py --no-build
```

**Step 1 — isolate the geometry.** Holding StericX's 2.1 Å geometric-centre
convention fixed and changing only the geometry source:

| Geometry source (2.1 Å centre) | Native descriptor \(R^2\) vs official Kraken |
|---|---:|
| Study 002 — RDKit/MMFF | 0.8626 |
| Study 003 — CREST/GFN2-xTB | 0.9254 |
| Kraken's own DFT geometries | \(R^2 = 0.9937\), Pearson \(r = 0.9993\), RMSE 0.5682 Å³ |

Agreement jumps to \(R^2 = 0.9937\) with a near-constant offset, localizing the
earlier shortfall to conformer geometry — not the voxel kernel.

**Step 2 — resolve the residual.** Kraken's descriptor code places the reference
metal 2.28 Å from phosphorus (`PL_dft_library_201027.py`), not 2.1 Å. Adopting
that documented convention closes the offset:

| Reference-metal distance | 11-ligand \(R^2\) | RMSE (Å³) | Slope |
|---|---:|---:|---:|
| 2.1 Å (geometry-isolating baseline) | 0.9937 | 0.5682 | 0.93 |
| 2.28 Å (Kraken's documented value) | **0.9986** | **0.2725** | 0.98 |

At Kraken's convention the kernel reproduces the published descriptor to
\(R^2 = 0.9986\), Pearson \(r = 0.9998\) — confirming the residual was a
coordination-centre convention difference, not the structures or the kernel.
Per-ligand table, parity plot, and provenance:
[`docs/study_004/STUDY_004.md`](docs/study_004/STUDY_004.md).

![Buried volume on Kraken DFT geometry](docs/study_004/kraken_dft_parity.png)

**Scaled to the entire public Kraken set.**
[`studies/study_004_scaled.py`](studies/study_004_scaled.py) runs the identical
kernel (2.28 Å convention) on every Kraken ligand with a published value and DFT
geometry:

```bash
uv run --extra science python studies/study_004_scaled.py --no-build
```

Across **1,541 chemically diverse ligands** (31,611 DFT conformers) the kernel
reproduces `vbur_max_delta_qvbur_min` with \(R^2 = 0.9852\), Pearson
\(r = 0.9927\), and a median absolute error of **0.11 Å³**. The error is
heavy-tailed, so the robust summary (correct for such a distribution, *not*
outlier removal) is a trimmed \(R^2 = 0.9897\) excluding the worst 1 % of
ligands. Details: [`docs/study_004/STUDY_004_SCALED.md`](docs/study_004/STUDY_004_SCALED.md).

![Scaled parity across the full Kraken set](docs/study_004/kraken_dft_scaled_parity.png)

**A genuine kernel bug, found and fixed at scale.** Scaling exposed a bug the
eleven trisubstituted Ni-hDA ligands could never trigger: the quadrant frame took
a donor's three *nearest heavy atoms* as its substituents, which silently
mis-framed primary and secondary phosphines (R–PH₂, R₂P–H) by discarding their
bonded hydrogens. Switching to covalent-radius bonding (hydrogens included)
removed every spurious zero and lifted \(R^2\) from 0.9649 to 0.9852 **without
discarding a single ligand** — the validated count in fact rose. The residual is
then fully characterized: tertiary phosphines (98.4 %) are unbiased, and the
entire remaining bias is confined to 24 primary/secondary phosphines, growing
~0.7 Å³ per P–H bond — the signature of the geometric lone-pair centre standing
in for Kraken's exact xTB centre
([`studies/study_004_frame_residual.py`](studies/study_004_frame_residual.py),
[`STUDY_004_RESIDUAL.md`](docs/study_004/STUDY_004_RESIDUAL.md)).

**The whole descriptor family, not one contrast.**
[`studies/study_004_vbur_family.py`](studies/study_004_vbur_family.py) compares StericX
against Kraken's published values for the *entire* `vbur` family — buried volume,
quadrant and octant extrema, near/far hemispheres, eight descriptors — across all
1,541 ligands (mean \(R^2 = 0.9925\)).
[`studies/study_004_sterimol.py`](studies/study_004_sterimol.py) does the same for Kraken's
published Sterimol once its coordination-axis convention is matched (a virtual
metal at the same 2.28 Å centre, with the +0.40 Å Verloop \(L\) correction; mean
\(R^2 = 0.9887\)). Two independent classical steric-descriptor families, each
reproduced to ~0.99 against Kraken's own numbers over the whole library. (These
are the headline rows summarized under [Validation](#validation).)

| ![Buried-volume family parity](docs/study_004/kraken_vbur_family_parity.png) | ![Sterimol parity](docs/study_004/kraken_sterimol_parity.png) |
|:---:|:---:|

### Study 005 — Pyramidalization descriptors

Kraken publishes two geometric pyramidalization descriptors for the phosphorus
donor: `pyr_P` (Radhakrishnan's dimensionless pyramidalization) and `pyr_alpha`
(the mean out-of-plane angle in degrees), both defined by morfeus's
`Pyramidalization` class. StericX reimplements the identical definitions
natively in Rust, reducing each to a closed form on the donor's three unit bond
vectors — `pyr_P = |det[â, b̂, ĉ]|` (with morfeus's `2 − P` acute correction) and
`pyr_alpha` as the mean signed out-of-plane angle. Both were verified against
morfeus to machine precision (4.4×10⁻¹⁶ and 2.8×10⁻¹⁴) before implementation.

[`studies/study_005_pyramidalization.py`](studies/study_005_pyramidalization.py) runs the
native kernel on Kraken's own DFT conformers and compares the `min` and `max`
conformer reductions (both weight-independent) against Kraken's published values
across the same 1,541 ligands.

```bash
uv run --extra science python studies/study_005_pyramidalization.py --no-build --workers 8
```

| Descriptor | Reduction | \(R^2\) | RMSE | Median abs. err |
|---|---|---:|---:|---:|
| `pyr_P` | min | 0.999983 | 0.000204 | 0.000028 |
| `pyr_P` | max | 0.999977 | 0.000174 | 0.000028 |
| `pyr_alpha` (°) | min | 0.999979 | 0.022463 | 0.001932 |
| `pyr_alpha` (°) | max | 0.999968 | 0.030468 | 0.001748 |

Mean \(R^2 = 0.99998\) — higher than the buried-volume family or Sterimol, and
expected rather than tuned: pyramidalization depends only on the three bond
directions, so it is insensitive to the virtual-metal centre, sphere radius, and
lone-pair conventions that bound the buried-volume agreement. The residual (RMSE
~2×10⁻⁴ for `pyr_P`, ~0.03° for `pyr_alpha`) is consistent with the 4-decimal
coordinate precision of the cached DFT SDFs and `f32` arithmetic — a
geometry-input floor, not a method difference. Full results:
[`docs/study_005/STUDY_005.md`](docs/study_005/STUDY_005.md).

![Pyramidalization parity](docs/study_005/kraken_pyramidalization_parity.png)

### Study 006 — Localizing the residual to the coordination center

Study 004 found the buried-volume residual is confined to 24 primary/secondary
phosphines and grows ~0.7 Å³ per P–H bond, and attributed it to StericX's
geometric lone-pair center standing in for Kraken's xTB localized-orbital center.
[`studies/study_006_residual_localization.py`](studies/study_006_residual_localization.py) turns that
attribution into a controlled test. It splits the six descriptors by their
dependence on the coordination center — buried volume and Sterimol `L`/`B1`/`B5`
are anchored on it; pyramidalization (`pyr_P`, `pyr_alpha`) is computed purely
from the three bond vectors and never references it — and measures the residual's
dependence on P–H count for each, on the *same* 1,541 ligands.

```bash
uv run --extra science python studies/study_006_residual_localization.py --no-build --workers 8
```

| Descriptor | Center | Std. residual slope per P–H |
|---|---|---:|
| buried volume | coupled | +1.51σ |
| Sterimol L | coupled | −1.86σ |
| Sterimol B1 | coupled | −1.12σ |
| Sterimol B5 | coupled | +1.67σ |
| pyramidalization `pyr_P` | free | +0.02σ |
| pyramidalization `pyr_alpha` | free | +0.05σ |

Every center-coupled descriptor shifts systematically with P–H count (mean 1.54 σ
per bond; signs differ because a mis-placed axis lengthens some measures and
shortens others), while both center-free pyramidalization descriptors are flat
(mean 0.04 σ) — an order-of-magnitude separation. Since pyramidalization shares
the donor, geometries, covalent-radius frame, and `f32` kernel with the buried
volume and differs *only* in never placing the coordination center, this rules
out the kernel, the frame, and the geometries: the residual is specifically the
geometric lone-pair center diverging from Kraken's xTB center. A bounded,
understood limit of the free pipeline (1.6% of the library) — not a bug, not
tuned away. Full results:
[`docs/study_006/STUDY_006.md`](docs/study_006/STUDY_006.md).

![Residual localization](docs/study_006/residual_localization.png)

### Study 007 — An independent second reaction model (cross-coupling)

Studies 001–006 validate StericX's descriptors and one reaction model (Ni-hDA).
Study 007 puts a descriptor to work on a **second, independent reaction study**:
Newman-Stonebraker, Sigman & Doyle (*Science* **2021**, *374*, 301,
[10.1126/science.abj4213](https://doi.org/10.1126/science.abj4213)) classify
monodentate phosphines as active/inactive across Ni cross-coupling reactions
using a single-node threshold on **%Vbur(min)** — a descriptor StericX already
reproduces (Study 004). [`studies/study_007_crosscoupling.py`](studies/study_007_crosscoupling.py)
tests whether StericX reproduces both the descriptor and the classifier on the
authors' own datasets (Reactions I–V, RS1).

**1. Descriptor fidelity.** Across **479 ligand–reaction points**, StericX's
independently-computed %Vbur(min) reproduces the published values at
**R² = 0.9992** (MAE 0.144%).

**2. Classifier reproduction (vs the paper's Table S11).** A single-node tree
(the paper's method) fit on StericX's descriptor recovers the same threshold,
direction, accuracy, and Matthews correlation (MCC) per reaction. `baseline` is
the majority-class accuracy (always predicting the larger class):

| Reaction | n | active | baseline | StericX thr / acc / MCC | Paper acc / MCC |
|---|---:|---:|---:|---|---:|
| I | 34 | 18 | 0.53 | 32.5 / 0.76 / 0.59 | 0.79 / 0.62 |
| II | 89 | 35 | 0.61 | 32.0 / 0.72 / 0.56 | 0.70 / 0.53 |
| III | 89 | 27 | 0.70 | 31.8 / 0.64 / 0.47 | 0.67 / 0.50 |
| IV | 89 | 30 | 0.66 | 32.0 / 0.66 / 0.50 | 0.64 / 0.45 |
| V | 89 | 50 | 0.56 | 50.8 / 0.65 / 0.36 | 0.66 / 0.36 |
| RS1 | 89 | 34 | 0.62 | 32.0 / 0.71 / 0.55 | 0.70 / 0.54 |

StericX's mean accuracy is **0.69** (mean MCC **0.50**), matching the paper's
**0.69 / 0.50**. Read these honestly: accuracy is a poor lens for imbalanced
data — for Reactions III/IV it sits at or below the majority-class baseline,
because the paper's 20:1 active weighting trades raw accuracy to avoid missing
active ligands. The honest metric is MCC (positive throughout, **0.36–0.59** — a
real but *moderate* signal; the bootstrap 95% CIs are wide at n = 34–89, and
Reaction V's reaches ~0.00). That's expected, not a shortfall: a deliberately
univariate model (one steric number) can't see electronics, substrate, or
conditions.

**3. Out-of-sample transferability.** Sections 1–2 are in-sample; this tests the
paper's actual claim — that ~32% %Vbur(min) is a *transferable* ligation cliff.
A single universal threshold pooled across all 479 points sits at **32.8%**
(acc 0.68, MCC 0.49). Leave-one-reaction-out (fit on five reactions, predict the
sixth out-of-sample) keeps the trained threshold near ~33% and OOS MCC positive
throughout (**0.41–0.54**) — the cliff genuinely transfers. The honest wrinkle is
Reaction V, whose *own* best-fit threshold jumps to **51%** (matching the paper's
51.5); yet the shared ~32% cliff still predicts V out-of-sample at MCC 0.47. V is
an outlier in threshold-space, exactly as documented, not a transfer failure.

**4. Independent geometry (removing the circularity).** Section 1's fidelity used
Kraken's *own* cached coordinates. Running StericX instead on the **paper's own
DFT free-ligand geometries** (18 ligands, matched to Kraken IDs by molecular
formula) reproduces the published %Vbur(boltz) at **R² = 0.9735** (offset
+0.12%), with **16/18** landing inside StericX's own Kraken-conformer %Vbur range
— a genuinely independent path (their structures → StericX's kernel, no shared
coordinates).

Together this extends the project from "reproduces Kraken's numbers" to
"reproduces, and stress-tests out-of-sample and off Kraken's own geometry, a real
reaction model beyond Ni-hDA." Full results:
[`docs/study_007/STUDY_007.md`](docs/study_007/STUDY_007.md).

The paper's supplementary data is copyrighted (AAAS) and is **not** redistributed
here; the study reads it locally (gitignored) and commits only StericX's own
computed values and the comparison.

![Cross-coupling %Vbur parity](docs/study_007/crosscoupling_vbur_parity.png)
![Independent-geometry %Vbur parity](docs/study_007/crosscoupling_independent_geom.png)

### Study 008 — A head-to-head speed benchmark vs morfeus

Fidelity (Studies 002/005) is only half of a reproduction's practical value; the
other half is cost. Study 008 times StericX against **morfeus**, the reference
Python implementation, computing the flagship **buried-volume** descriptor on the
same **1,546-ligand** library, on the **same single CPU core**, both timed
end-to-end from a warm cache.

| | Wall time | Throughput | Relative |
|---|---:|---:|---:|
| morfeus (Python) | 19.2 s | 81 / s | 1× |
| **StericX** (Rust binary) | **1.4 s** | **1,110 / s** | **≈14×** |

StericX is **~14× faster on a single core** — and the comparison is
*conservative*: StericX computes Sterimol and pyramidalization in the same timed
pass, while morfeus is timed for buried volume alone. The speedup is not from
computing something cheaper: of the 1,534 phosphines morfeus can frame, **1,518
agree to R² = 0.999999** (max |diff| 0.42 %Vbur). The 16 that differ are
frame-topology cases — morfeus's naive "3 nearest heavy atoms" vs StericX's
covalent-radius bonding (the Study 006 frame issue) — shown separately, not
averaged away. StericX also ships as a single **1.9 MB** binary with zero runtime
dependencies. Full results:
[`docs/study_008/STUDY_008.md`](docs/study_008/STUDY_008.md).

![Speed benchmark vs morfeus](docs/study_008/speed_benchmark.png)

---

## Roadmap

```mermaid
graph TD
    subgraph Done["✅ Completed"]
        D1[Sterimol + buried-volume kernels]
        D2[morfeus fidelity · R² ≥ 0.9999]
        D3[Full-set Kraken reproduction · 1,541 ligands]
        D4[Frame bug fixed · residual characterized]
        D5[Ni-hDA selectivity model · Study 001]
        D6[CREST/xTB quantum backend · Study 003]
        D7[SIMD inference · .sigpack v1/v2 storage]
        D8[CLI descriptors tool · benchmarks · green CI]
        D9[Pyramidalization pyr_P / pyr_alpha · mean R² 0.99998 · Study 005]
        D10[Residual localized to the centre · Study 006]
        D11[Second reaction model · Ni cross-coupling · Study 007]
        D12[Speed benchmark vs morfeus · ≈14× single-core · Study 008]
    end
    subgraph Next["🔭 Planned / optional"]
        P1["Prospective 10-candidate deck<br/>frozen; awaits experimental measurement"]
        P2["Zenodo DOI / v0.1.0 release<br/>metadata prepared, currently deferred"]
    end
    Done --> Next
```

Planned items are scoped but not started. The prospective candidate deck stays
target-free until lab data exists; a Zenodo DOI / release remains available if
the metadata is ever published.

---

## CLI Reference

All commands are subcommands of the release binary
(`./target/release/stericx`).

### `descriptors` — featurize a ligand (start here)

Point `descriptors` at any phosphine geometry — an `.xyz`, or an `.sdf`/`.mol`
(one or many conformers) — and get its Sterimol and buried-volume descriptors.
The donor and its substituents are detected from the geometry, so there is no
reaction CSV to build and no atom indices to look up.

Descriptors use Kraken's convention by default (3.5 Å sphere, Bondi radii ×1.17,
2.28 Å reference-metal distance); override with `--sphere-radius`,
`--radii-scale`, or `--center-distance`. A multi-model SDF is treated as a
conformer ensemble, and per-file values are averaged over conformers with
Kraken's `max_delta_qvbur_min` (the minimum over the ensemble) reported as the
headline buried-volume descriptor.

Featurize a whole folder into a spreadsheet-ready table, or emit JSON:

```bash
./target/release/stericx descriptors --format csv  ligands/*.xyz > descriptors.csv
./target/release/stericx descriptors --format json ligands/*.sdf
```

By default Sterimol is measured along the donor→substituent bond. For a
metal-bound ligand, pass `--sterimol-axis coordination` to measure along the
coordination axis instead (a virtual metal 2.28 Å from the donor on the lone
pair, with the +0.40 Å Verloop \(L\) correction) — the convention Kraken and the
ligand-descriptor literature use, which StericX reproduces to a mean
\(R^2 = 0.9887\) against Kraken's published values across 1,541 ligands.

Non-phosphine donors work too — pass `--donor-element N`, or name the atom
directly with `--donor-index` when a structure has more than one candidate. Files
that are not valid trivalent donors are reported on stderr and skipped, so a
batch run over a messy folder still completes.

### `parse` — pack reaction structures

```bash
./target/release/stericx parse \
  --csv data/reactions_raw.csv \
  --xyz-dir data \
  --output data/reactions.sigpack
```

The reaction CSV contains:

```text
Reaction_ID,Ligand_XYZ_Path,Attach_Atom_Idx,Primary_Bond_Vector_Idx,NBO_Charge,IR_Frequency,Temp_K,Exp_ddG_kcal_mol
```

StericX parses each XYZ structure, aligns the attachment vector, calculates
Sterimol \(L\), \(B_1\), \(B_5\), Boltzmann-averages conformer ensembles, stores
their minimum/maximum envelopes in the reserved binary slots, and exports one
64-byte `PackedReactionRecord` per reaction.

### `buried-volume` — coordination-aware buried volume

```bash
./target/release/stericx buried-volume \
  --csv data/reactions_quantum.csv \
  --xyz-dir data \
  --output data/reactions_v2.sigpack \
  --per-conformer-output data/buried_volume_conformers.csv \
  --require-explicit-centers
```

The command places a virtual metal centre 2.1 Å from the donor, scans all three
donor-substituent quadrant orientations, and writes total, quadrant, octant,
near/far, and maximum-adjacent-quadrant descriptors. Per-conformer values are
reduced to Boltzmann, minimum, maximum, range, and
minimum-buried-volume-conformer statistics. When the CSV supplies
`Conformer_Coordination_Centers_Angstrom`, each conformer can use an xTB
localized-molecular-orbital centre; `--require-explicit-centers` prevents an
accidental fallback to a geometrically inferred direction.

### `predict` — parallel SIMD regression

Weights are an eight-number JSON array; the eight model features are
`[1, L, B1, B5, nbo_charge, B1*nbo_charge, B5*nbo_charge, ir_freq]`.

```bash
./target/release/stericx predict \
  --data data/reactions.sigpack \
  --weights weights.json
```

The command reports mapping latency, prediction latency, throughput, Rayon
thread count, RSS memory, and MSE against the packed experimental
\(\Delta\Delta G^\ddagger\) labels.

### `fit` / `evaluate` — freeze and reveal a scientific model

```bash
./target/release/stericx fit \
  --data data/reactions.sigpack \
  --metadata data/reactions_raw.csv \
  --output docs/study_001/stericx_model.json \
  --predictions docs/study_001/stericx_frozen_predictions.csv
```

The fitter learns scaling from training rows only, limits model complexity to
fewer than one term per three observations, rejects descriptor pairs with
\(|r| > 0.95\), and performs BIC-constrained forward selection. The model
artifact includes fixed-feature leave-one-out diagnostics, nested ridge and LASSO
baselines, bootstrap coefficient intervals, response permutation tests, VIF
values, correlations, applicability ranges, and leave-one-scaffold-group-out
validation using `Ligand_Group`. Non-training predictions are written without
experimental targets; reveal them separately:

```bash
./target/release/stericx evaluate \
  --data data/reactions.sigpack \
  --metadata data/reactions_raw.csv \
  --model docs/study_001/stericx_model.json \
  --predictions docs/study_001/stericx_frozen_predictions.csv \
  --output docs/study_001/stericx_evaluation.json
```

### `simulate` — selectivity from kinetics

```bash
./target/release/stericx simulate --ddg 1.82 --temp 298.15
```

Output includes the Eyring rate constant, major/minor percentages, R:S product
distribution, enantiomeric excess, execution time, and memory metrics.

---

## Internal Architecture

### Source layout

```text
src/
├── geometry/     XYZ/SDF parsing, Sterimol, buried-volume, and pyramidalization
├── storage/      Cache-aligned schema and memory-mapped .sigpack I/O
├── model/        Scientific fitting, feature interactions, and SIMD inference
├── kinetics/     Eyring rates and enantiomeric distributions
└── main.rs       clap command-line interface (incl. `descriptors`)

studies/          Reproduction study drivers (Python → docs/study_00N/)
├── study_001_ni_hda.py                 Ni-hDA enantioselectivity model (Study 001)
├── study_002_buried_volume.py          buried-volume fidelity vs morfeus (Study 002)
├── study_003_quantum_geometry.py       CREST/xTB geometry ensemble (Study 003)
├── study_004_reproduction.py           11-ligand Kraken DFT reproduction (Study 004)
├── study_004_scaled.py                 full 1,541-ligand scaled reproduction
├── study_004_vbur_family.py            whole buried-volume family vs Kraken
├── study_004_sterimol.py               Sterimol vs Kraken, coordination axis
├── study_004_frame_residual.py         residual anatomy by phosphine class
├── study_005_pyramidalization.py       pyr_P / pyr_alpha vs Kraken (Study 005)
├── study_006_residual_localization.py  residual localized to the centre (Study 006)
├── study_007_crosscoupling.py          independent cross-coupling model (Study 007)
└── study_008_speed_benchmark.py        speed benchmark vs morfeus (Study 008)

scripts/          Support utilities (data prep, validation, quantum backend)
├── prepare_data.py                     reaction-data preparation
├── prepare_quantum_data.py             CREST/xTB geometry preparation
├── validate_stericx.py                 Sterimol fidelity vs morfeus
├── stericx_quantum.py                  xTB/CREST backend wrapper
├── freeze_prospective_deck.py          frozen prospective candidate deck
└── preregister_prediction.py           pre-registered falsifiable forecast
```

Per-study results, tables, and parity figures live under `docs/study_00N/`; a
manuscript-style write-up is in
[`docs/REPRODUCTION_REPORT.md`](docs/REPRODUCTION_REPORT.md) and the release
history in [`CHANGELOG.md`](CHANGELOG.md).

### Binary record layout

Every legacy `.sigpack` v1 row occupies exactly one cache line:

```rust
#[repr(C, align(64))]
pub struct PackedReactionRecord {
    pub l: f32,
    pub b1: f32,
    pub b5: f32,
    pub nbo_charge: f32,
    pub ir_freq: f32,
    pub temp_k: f32,
    pub exp_ddg: f32,
    pub reserved: [f32; 9],
}
```

The format is a flat native-endian POD matrix intended for high-throughput
processing on the machine that generated it. Version two adds a 64-byte
`SigPackHeaderV2` and stores `PackedReactionRecordV2` as two cache lines (the
unchanged v1 record followed by one 64-byte buried-volume block); dedicated v1
and v2 memory-mapped readers preserve backward compatibility.

```rust
#[repr(C, align(64))]
pub struct PackedReactionRecordV2 {
    pub reaction: PackedReactionRecord,
    pub buried_volume: PackedBuriedVolumeRecord,
}
```

---

## References and Attribution

StericX is an **independent** reproduction and reimplementation. It is not
affiliated with, endorsed by, or produced by the Sigman or Reisman groups, and it
reuses only their publicly released data and descriptor definitions. If you use
the descriptors, datasets, or reaction models reproduced here, please cite the
original scientific work:

- Gensch, T.; dos Passos Gomes, G.; Friederich, P.; Peters, E.; Gaudin, T.;
  Pollice, R.; Jorner, K.; Nigam, A.; Lindner-D'Addario, M.; Sigman, M. S.;
  Aspuru-Guzik, A. "A Comprehensive Discovery Platform for Organophosphorus
  Ligands for Catalysis." *J. Am. Chem. Soc.* **2022**, *144* (3), 1205–1217.
  DOI: [10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718). (The
  "Kraken" descriptor platform.)
- Cadge, J. A.; Lozano, C.; Merriman, M. T.; Oblad, P.; Sigman, M. S.; Reisman,
  S. E. "A Data Science-Guided Approach for the Development of Nickel-Catalyzed
  Homo-Diels–Alder Reactions." *J. Am. Chem. Soc.* **2025**, *147* (34),
  31175–31186. DOI: [10.1021/jacs.5c09948](https://doi.org/10.1021/jacs.5c09948).
  (The Ni-hDA reaction and enantioselectivity dataset reproduced in Study 001.)
- Newman-Stonebraker, S. H.; Smith, S. R.; Borowski, J. E.; Peters, E.; Gensch,
  T.; Johnson, H. C.; Sigman, M. S.; Doyle, A. G. "Univariate Classification of
  Phosphine Ligation State and Reactivity in Cross-Coupling Catalysis."
  *Science* **2021**, *374* (6565), 301–308.
  DOI: [10.1126/science.abj4213](https://doi.org/10.1126/science.abj4213). (The
  cross-coupling reactivity classifier reproduced in Study 007; its supplementary
  data is copyrighted by AAAS and is **not** redistributed in this repository.)

Descriptor reference values are compared against
[`morfeus`](https://github.com/digital-chemistry-laboratory/morfeus) and the
public [Kraken](https://kraken.cs.toronto.edu) tabulation.

### How to cite StericX

A machine-readable citation is in [`CITATION.cff`](CITATION.cff) (GitHub renders a
"Cite this repository" button from it). To cite the software directly:

> Rumenovski, A. *StericX: physical-organic molecular featurization and
> reaction-selectivity inference*, version 0.1.0, 2026.
> https://github.com/AndrejRumenovski/StericX

<!-- ARCHIVE-DOI: replaced with the Zenodo `Cite all versions` DOI when the
     v0.1.0 GitHub release is archived; see RELEASING.md. -->

Once archived, an immutable, version-specific DOI (Zenodo) is preferred for
reproducibility. StericX is an independent reproduction; if you use the
descriptors, datasets, or reaction models it reproduces, please also cite the
original scientific works above.

---

## License

Licensed under either of

- Apache License, Version 2.0 ([`LICENSE-APACHE`](LICENSE-APACHE) or
  <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([`LICENSE-MIT`](LICENSE-MIT) or
  <http://opensource.org/licenses/MIT>)

at your option (SPDX: `MIT OR Apache-2.0`). Unless you explicitly state
otherwise, any contribution intentionally submitted for inclusion in the work by
you, as defined in the Apache-2.0 license, shall be dual licensed as above,
without any additional terms or conditions.
