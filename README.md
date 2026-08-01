# StericX

**A native-Rust engine for physical-organic molecular featurization and reaction-selectivity inference.**

[![CI](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21726666.svg)](https://doi.org/10.5281/zenodo.21726666)

StericX computes 3D **Sterimol** and coordination-aware **buried-volume** descriptors
straight from Cartesian coordinates and turns them into interpretable reaction models —
as a single, dependency-free binary. It is an **independent, from-scratch reproduction**
of the Kraken / morfeus steric-descriptor toolchain, validated against both `morfeus`
(numerically) and Kraken's *own published values* across the full 1,541-ligand library.

📄 [Manuscript write-up](docs/REPRODUCTION_REPORT.md) · 🖼️ [One-page visual overview](docs/results.html) · 🔁 [Clone-to-results walkthrough](REPRODUCE.md)

---

## What is StericX? (the 1-minute version)

**What it is.** Point it at an `.xyz`/`.sdf`/`.mol` file and it auto-detects the donor
atom and prints Sterimol, buried-volume, and pyramidalization descriptors — no Python
runtime, no reaction CSV, no atom indices to look up. The same binary fits interpretable
linear models and converts predicted ΔΔG‡ into product ratios via the Eyring equation.

**Why it matters.** The standard steric-descriptor stack (Kraken, morfeus) is Python and
tuned to one research workflow. StericX reimplements it from scratch in Rust — faster,
portable, and, above all, *checked*. The project's entire purpose is an **honest,
defensible reproduction**: every number below is a generated measurement, and every
failed gate is kept in view rather than hidden.

**Results.**

| What | Result |
|---|---|
| Descriptor fidelity vs `morfeus` (identical geometries) | Sterimol R² ≥ 0.9999 · buried volume R² = 1.000000 |
| Reproducing Kraken's **published** values | R² = 0.9852 across **1,541 ligands / 31,611 DFT conformers** |
| Two published reaction models reproduced | Ni **and** Pd cross-coupling cliff — *both* directions |
| Speed vs `morfeus` (single core) | **~14× faster** — and it holds at 31,721 conformers |
| Honest by construction | compact native descriptors *fail* Ni-hDA (LOO Q² ≈ 0.002) — reported, led-with |

---

## Quick Start

One script installs the toolchains, builds the binary, and smoke-tests it:

```bash
./bootstrap.sh                 # uv + Rust + build + smoke test
./bootstrap.sh --with-quantum  # also fetch CREST 2.12 / xTB 6.4.0 (Study 003)
```

Then featurize any phosphine geometry — the donor and its substituents are detected from
the structure:

```bash
./target/release/stericx descriptors data/xyz/SIG-NIHDA-401_9d42bff1.xyz
```

```text
data/xyz/SIG-NIHDA-401_9d42bff1.xyz
  donor          P (atom 14)
  substituents   C, C, C
  Sterimol      L 7.76   B1 1.73   B5 8.75   Å
  buried volume  Vbur 24.3%   (43.6 Å³)
                 qvbur_min 9.73   qvbur_max 13.53   max_delta_qvbur 3.80   Å³
  pyramidalization pyr_P 0.932   pyr_alpha 17.10°
```

Batch a folder into CSV/JSON with `--format`, target non-phosphine donors with
`--donor-element`, or switch to the coordination axis with `--sterimol-axis coordination`.
Full command set and the reproducible Python environment: see [Documentation](#documentation).

---

## Highlights

- **Exact geometry kernel.** Sterimol L/B₁/B₅ reproduce `morfeus` to R² ≥ 0.9999; the
  buried-volume voxel engine matches to R² = 1.000000 (worst mean relative error 8×10⁻⁶ %);
  pyramidalization matches to machine precision (~10⁻¹⁴).
- **Validated at library scale.** Kraken's published `vbur_max_delta_qvbur_min`
  reproduced across 1,541 ligands / 31,611 DFT conformers (median abs. error 0.11 Å³);
  the whole buried-volume family (8 descriptors) at mean R² = 0.9925.
- **Two reaction models, both cliff directions.** The Newman-Stonebraker cross-coupling
  classifier is reproduced for Ni (small-active) *and* Pd (bulky-active) — descriptor
  R² = 0.999, including datasets from other groups.
- **Failures stay visible.** Compact native descriptors *don't* model Ni-hDA (Q² ≈ 0.002);
  a real kernel frame-bug was fixed at scale **without dropping a single ligand**; the
  remaining residual is fully characterized and localized to one convention difference.
- **Fast, portable core.** ~14× faster than `morfeus` single-core; 318 M SIMD evals/s over
  a memory-mapped million-record matrix; ships as one 1.9 MB binary with zero runtime deps.
- **Reproducible & falsifiable.** One-command bootstrap, checksum-pinned CREST/xTB,
  content-addressed caches, a SHA-256-frozen **pre-registered** prediction, green CI,
  Zenodo DOI, dual MIT/Apache licensing.

---

## Validation

Two independent references anchor every descriptor claim: `morfeus` for numerical fidelity
on identical geometries, and Kraken's own *published* values for chemical accuracy at
library scale.

**Numerical fidelity vs `morfeus`** ([`scripts/validate_stericx.py`](scripts/validate_stericx.py), 11 structures):

| Parameter | R² | RMSE (Å) | Status |
|---|---:|---:|---|
| Sterimol L | 1.000000 | 0.000000 | Validated |
| Sterimol B₁ | 0.999959 | 0.010539 | Validated (1° angular scan) |
| Sterimol B₅ | 1.000000 | 0.000001 | Validated |
| Buried-volume voxel kernel | 1.000000 | — (rel. err 8×10⁻⁶ %) | Exact |

**Reproducing Kraken's published values** (unchanged kernel on Kraken's own DFT geometries, 2.28 Å convention):

| Target vs published Kraken | Scale | Result |
|---|---|---:|
| `vbur_max_delta_qvbur_min` | 11 ligands | R² = 0.9986 |
| `vbur_max_delta_qvbur_min` | 1,541 ligands / 31,611 conformers | R² = 0.9852 · medAE 0.11 Å³ |
| Buried-volume family (8 descriptors) | 1,541 ligands | mean R² = 0.9925 |
| Sterimol L/B₁/B₅ (coordination axis) | 1,541 ligands | mean R² = 0.9887 |
| Pyramidalization `pyr_P`, `pyr_alpha` | 1,541 ligands | mean R² = 0.99998 |

The geometry engine is exact; the open question is always the *input geometry*, which the
studies below isolate. Derivations, per-ligand tables, and parity plots are in each
[study's docs](#scientific-studies).

---

## Scientific studies

Nine preregistered studies build from a small published reaction family up to the full
Kraken library and two independent reaction models. Each writes a complete model card,
frozen prediction hashes, raw comparisons, plots, and machine-readable results under
`docs/study_00N/`. **Passed and failed gates are both retained.**

<details>
<summary><b>Expand the nine studies</b></summary>

| # | Study | What it shows | Full results |
|---|---|---|---|
| **001** | Ni-hDA enantioselectivity model | The published Kraken descriptor reproduces the selectivity (LOO Q² = 0.75); StericX's *compact native* descriptors deliberately do **not** (Q² ≈ 0.002) — the project leads with this ablation. | [STUDY_001](docs/study_001/STUDY_001.md) |
| **002** | Coordination-aware buried-volume fidelity | The Rust voxel kernel equals `morfeus` to R² = 1.000000 on identical geometries; cheap RDKit/MMFF conformers fall short (R² = 0.86) — a failed gate kept visible. | [STUDY_002](docs/study_002/STUDY_002.md) |
| **003** | Quantum geometry & prospective validation | A checksum-pinned CREST/xTB backend, plus a **frozen, falsifiable 10-ligand prediction** committed by SHA-256 before any measurement. | [STUDY_003](docs/study_003/STUDY_003.md) · [PREREGISTRATION](docs/study_003/PREREGISTRATION.md) |
| **004** | Reproducing Kraken's published descriptors on DFT geometries | On Kraken's own DFT geometries at the 2.28 Å convention: R² = 0.9986 (11 ligands) → **0.9852 across 1,541 ligands**. A genuine kernel frame-bug was found and fixed at scale *without dropping a ligand*; the residual is fully characterized. | [STUDY_004](docs/study_004/STUDY_004.md) · [scaled](docs/study_004/STUDY_004_SCALED.md) · [residual](docs/study_004/STUDY_004_RESIDUAL.md) |
| **005** | Pyramidalization descriptors | `pyr_P` and `pyr_alpha` reduced to closed forms, verified vs `morfeus` to ~10⁻¹⁴; mean R² = 0.99998 across 1,541 ligands. | [STUDY_005](docs/study_005/STUDY_005.md) |
| **006** | Localizing the residual to the coordination centre | A controlled test proving the buried-volume residual is a *coordination-centre artefact*: centre-coupled descriptors carry it (1.5σ/P–H), centre-free pyramidalization does not (0.04σ). | [STUDY_006](docs/study_006/STUDY_006.md) |
| **007** | Independent second reaction model — Ni cross-coupling | Reproduces the Newman-Stonebraker classifier (descriptor R² = 0.9992); the ~32% %Vbur(min) ligation cliff transfers *out-of-sample* (LORO MCC 0.41–0.54) and off Kraken's own geometry (R² = 0.9735). | [STUDY_007](docs/study_007/STUDY_007.md) |
| **008** | Head-to-head speed benchmark vs `morfeus` | ~14× faster single-core, same numbers (R² = 0.999999). Holds at **31,721 conformers** (~20× the set): 13.8× faster, R² = 0.99999983. | [STUDY_008](docs/study_008/STUDY_008.md) · [scale check](docs/study_008_all_conformers/STUDY_008.md) |
| **009** | The other direction of the cliff — Pd cross-coupling | The *opposite* (bulky-active) regime reproduced for **6/6** reactions (mean MCC 0.64 vs 0.67), including **two datasets from other groups**. Honest about the failures the paper itself flags (VII, XII). | [STUDY_009](docs/study_009/STUDY_009.md) |

A manuscript-style narrative of all nine is in [`docs/REPRODUCTION_REPORT.md`](docs/REPRODUCTION_REPORT.md).

</details>

---

## Documentation

- 📄 **[Manuscript write-up](docs/REPRODUCTION_REPORT.md)** — the full narrative, honesty led-with.
- 🖼️ **[One-page visual overview](docs/results.html)** — the whole scope at a glance.
- 🔁 **[REPRODUCE.md](REPRODUCE.md)** — clone to results in one pass.
- 📦 **[RELEASING.md](RELEASING.md)** · 📝 **[CHANGELOG.md](CHANGELOG.md)** — release runbook and history.

<details>
<summary><b>CLI reference</b> — all subcommands of <code>./target/release/stericx</code></summary>

### `descriptors` — featurize a ligand (start here)

Point it at any `.xyz`/`.sdf`/`.mol` (one or many conformers); the donor and substituents
are detected from the geometry. Kraken's convention by default (3.5 Å sphere, Bondi radii
×1.17, 2.28 Å reference metal); override with `--sphere-radius`, `--radii-scale`,
`--center-distance`. Multi-model SDFs are treated as conformer ensembles.

```bash
./target/release/stericx descriptors --format csv  ligands/*.xyz > descriptors.csv
./target/release/stericx descriptors --format json ligands/*.sdf
```

Sterimol is measured along the donor→substituent bond by default; `--sterimol-axis
coordination` uses the metal-bound coordination axis (a virtual metal 2.28 Å from the donor
with the +0.40 Å Verloop L correction). Non-phosphine donors: `--donor-element N`, or
`--donor-index` to name the atom directly.

### `parse` — pack reaction structures

```bash
./target/release/stericx parse --csv data/reactions_raw.csv --xyz-dir data --output data/reactions.sigpack
```

Reads a `Reaction_ID,Ligand_XYZ_Path,Attach_Atom_Idx,Primary_Bond_Vector_Idx,NBO_Charge,IR_Frequency,Temp_K,Exp_ddG_kcal_mol`
CSV, computes Sterimol, Boltzmann-averages conformer ensembles, and exports one 64-byte
`PackedReactionRecord` per reaction.

### `buried-volume` — coordination-aware buried volume

```bash
./target/release/stericx buried-volume --csv data/reactions_quantum.csv --xyz-dir data \
  --output data/reactions_v2.sigpack --per-conformer-output data/bv_conformers.csv --require-explicit-centers
```

Places a virtual metal centre, scans the three donor-substituent quadrant orientations, and
writes total/quadrant/octant/near-far/max-adjacent descriptors reduced to Boltzmann, min,
max, and range. `--require-explicit-centers` forces xTB LMO centres from the CSV instead of
a geometric fallback.

### `predict` — parallel SIMD regression

Weights are an eight-number JSON array over features `[1, L, B1, B5, nbo_charge,
B1*nbo_charge, B5*nbo_charge, ir_freq]`.

```bash
./target/release/stericx predict --data data/reactions.sigpack --weights weights.json
```

Reports mapping/prediction latency, throughput, Rayon thread count, RSS, and MSE against the
packed experimental ΔΔG‡ labels.

### `fit` / `evaluate` — freeze and reveal a scientific model

```bash
./target/release/stericx fit --data data/reactions.sigpack --metadata data/reactions_raw.csv \
  --output docs/study_001/stericx_model.json --predictions docs/study_001/stericx_frozen_predictions.csv
```

Learns scaling from training rows only, caps complexity below one term per three
observations, rejects `|r| > 0.95` descriptor pairs, and does BIC-constrained forward
selection. The model card includes LOO diagnostics, ridge/LASSO baselines, bootstrap
coefficient intervals, permutation tests, VIF, and leave-one-scaffold-group-out validation.
Non-training predictions are written target-free; `evaluate` reveals them separately.

### `simulate` — selectivity from kinetics

```bash
./target/release/stericx simulate --ddg 1.82 --temp 298.15
```

Outputs the Eyring rate constant, major/minor percentages, R:S distribution, and ee.

</details>

<details>
<summary><b>Internal architecture</b> — source layout & binary record format</summary>

```text
src/
├── geometry/     XYZ/SDF parsing, Sterimol, buried-volume, and pyramidalization
├── storage/      Cache-aligned schema and memory-mapped .sigpack I/O
├── model/        Scientific fitting, feature interactions, and SIMD inference
├── kinetics/     Eyring rates and enantiomeric distributions
└── main.rs       clap command-line interface (incl. `descriptors`)

studies/          Reproduction study drivers (Python → docs/study_00N/), study_001 … study_009
scripts/          Support utilities: prepare_data, prepare_quantum_data, validate_stericx,
                  stericx_quantum, freeze_prospective_deck, preregister_prediction
```

Reaction observations use an exact 64-byte `#[repr(C, align(64))]` layout; `bytemuck`
provides POD byte casting and `memmap2` exposes `.sigpack` files as borrowed record slices
with no per-record deserialization. An eight-lane AVX2 kernel drives the RegressX dot
product (portable fallback otherwise), and `rayon` distributes inference rows across workers
with no shared mutable state in the hot loop.

```rust
#[repr(C, align(64))]
pub struct PackedReactionRecord {
    pub l: f32, pub b1: f32, pub b5: f32,
    pub nbo_charge: f32, pub ir_freq: f32, pub temp_k: f32, pub exp_ddg: f32,
    pub reserved: [f32; 9],
}
```

Version two adds a 64-byte `SigPackHeaderV2` and stores a second cache line (a buried-volume
block) per record; dedicated v1/v2 mmap readers preserve backward compatibility.

</details>

<details>
<summary><b>Performance benchmarks</b> (<a href="benchmark_linux.sh">benchmark_linux.sh</a>)</summary>

| Workload | Items | Throughput | Peak RAM |
|---|---:|---:|---:|
| Sterimol extraction | 10,000 molecules | 62,126 mol/s | 7.75 MB |
| Binary `.sigpack` packing | 10,000 records | 56.82 M records/s | 7.75 MB |
| RegressX SIMD prediction | 1,000,000 records | 317.76 M evals/s | 58.94 MB |

Machine-readable results (CPU utilization, page faults) in
[`docs/benchmark_results.json`](docs/benchmark_results.json). Numbers are specific to the
machine that generated them; the ratios are the portable quantity.

</details>

<details>
<summary><b>Roadmap</b></summary>

Everything through v0.1.0 is shipped: the descriptor kernels, full-set Kraken reproduction,
the frame-bug fix and residual characterization, two reaction models (Ni + Pd cross-coupling),
pyramidalization, the speed benchmark, and the Zenodo release. The one open item is the
**frozen prospective 10-candidate deck**, which stays target-free until experimental
measurement exists; a measured outcome would seed a v0.2.0 under the existing concept DOI
(never a refit — pre-registration integrity).

</details>

---

## References & citation

StericX is an **independent** reproduction — not affiliated with, endorsed by, or produced
by the Sigman or Reisman groups; it reuses only their publicly released data and descriptor
definitions. If you use the descriptors, datasets, or reaction models reproduced here,
please cite the original work:

- Gensch et al. "A Comprehensive Discovery Platform for Organophosphorus Ligands for
  Catalysis." *J. Am. Chem. Soc.* **2022**, *144* (3), 1205–1217.
  DOI: [10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718) — the "Kraken" platform.
- Cadge et al. "A Data Science-Guided Approach for the Development of Nickel-Catalyzed
  Homo-Diels–Alder Reactions." *J. Am. Chem. Soc.* **2025**, *147* (34), 31175–31186.
  DOI: [10.1021/jacs.5c09948](https://doi.org/10.1021/jacs.5c09948) — Ni-hDA (Study 001).
- Newman-Stonebraker et al. "Univariate Classification of Phosphine Ligation State and
  Reactivity in Cross-Coupling Catalysis." *Science* **2021**, *374* (6565), 301–308.
  DOI: [10.1126/science.abj4213](https://doi.org/10.1126/science.abj4213) — cross-coupling
  (Studies 007/009); its SI is copyrighted by AAAS and is **not** redistributed here.

**To cite StericX** (machine-readable form in [`CITATION.cff`](CITATION.cff)):

> Rumenovski, A. *StericX: physical-organic molecular featurization and reaction-selectivity
> inference*, version 0.1.0, 2026. Zenodo.
> DOI: [10.5281/zenodo.21726666](https://doi.org/10.5281/zenodo.21726666).

The concept DOI [10.5281/zenodo.21726666](https://doi.org/10.5281/zenodo.21726666) resolves
to the latest version; the v0.1.0 version DOI is
[10.5281/zenodo.21726667](https://doi.org/10.5281/zenodo.21726667).

---

## License

Licensed under either of [Apache-2.0](LICENSE-APACHE) or [MIT](LICENSE-MIT) at your option
(SPDX: `MIT OR Apache-2.0`). Unless you state otherwise, any contribution intentionally
submitted for inclusion shall be dual licensed as above, without additional terms.
