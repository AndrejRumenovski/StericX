# StericX

**StericX is a native-Rust engine that reproduces the published molecular steric
descriptors used in modern catalysis research — with an order-of-magnitude speed
improvement over the existing Python implementations.**

[![CI](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/StericX/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21726666.svg)](https://doi.org/10.5281/zenodo.21726666)

Chemists use *steric descriptors* — numbers that capture how big and what shape a ligand
is — to predict how catalysts behave. The reference tools (Kraken, morfeus) are Python and
tuned to specific research workflows. StericX is an **independent, from-scratch
reproduction** of that toolchain as a single, dependency-free binary: it computes the same
**Sterimol** and **buried-volume** descriptors straight from atomic coordinates, and is
validated against both `morfeus` (numerically) and Kraken's *own published values* across
the full 1,541-ligand library — with every failure kept in view.

📄 [Manuscript write-up](docs/REPRODUCTION_REPORT.md) · 🖼️ [One-page visual overview](docs/results.html) · 🔁 [Clone-to-results walkthrough](REPRODUCE.md)

---

## What is StericX?

**What it is.** Point it at an `.xyz`/`.sdf`/`.mol` file and it auto-detects the donor
atom and prints Sterimol, buried-volume, and pyramidalization descriptors — no Python
runtime, no reaction CSV, no atom indices to look up. It also **searches**: point it at a
ligand and a library and it ranks the most sterically similar candidates, under constraints
like `--vbur 30:35 --b5-max 8` or `--less-bulky`; `compare` puts ligands side by side and
`db build` precomputes a reusable descriptor database. The same binary fits
interpretable linear models, **screens** a whole library through a fitted model — predicted
performance, an uncertainty band, and an applicability-domain warning per ligand — and
converts predicted ΔΔG‡ into product ratios via the Eyring equation.

**Why it matters.** The standard steric-descriptor stack (Kraken, morfeus) is Python and
tuned to one research workflow. StericX reimplements it from scratch in Rust — faster,
portable, and, above all, *checked*. The project's entire purpose is an **honest,
defensible reproduction**: every number below is a generated measurement, and every
failed gate is kept in view rather than hidden.

---

## Key Results

Six numbers carry the project. Everything else — RMSE, Pearson, slope, intercept, median
AE, trimmed summaries, per-parameter tables — lives in the [study docs](#scientific-studies).
These are the ones that matter:

Correctness first — speed is the last row on purpose:

| What it means | Number |
|---|---|
| The geometry kernel **equals** the reference tool (`morfeus`) on identical structures | **R² = 1.000000** |
| It reproduces Kraken's **published** descriptors across the full public library | **R² = 0.9852** · 1,541 ligands |
| It reproduces two independent **published reaction models** (Ni + Pd cross-coupling) | **descriptor R² ≈ 0.999** |
| **Honest limit:** compact native descriptors *don't* model Ni-hDA — reported, not hidden | **LOO Q² ≈ 0.002** |
| A **falsifiable prediction**, frozen before any measurement | **10 ligands · SHA-256** |
| Only then — on one CPU core it is far faster than `morfeus`, computing the same numbers | **~14×** |

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

## How it's validated

Two independent references anchor every claim: **`morfeus`** for numerical fidelity on
identical geometries, and **Kraken's own published values** for chemical accuracy across the
full library. The geometry engine itself is exact — so the open question is always the
*input geometry*, which the studies isolate one variable at a time, keeping failed gates in
view. The full per-parameter tables (R², RMSE, slope, intercept, median AE, trimmed
summaries) live in each study's card under `docs/`, not here.

The kernel is element-generic, not phosphorus-only: `--donor-element N` reproduces
morfeus's pyramidalization, buried-volume, and Sterimol descriptors on nitrogen donors
to the same fidelity as phosphorus
([non-phosphorus validation](docs/validation/NONP_DONOR.md)).

---

## Technical challenges

The hard parts weren't writing the code — they were the small, exacting decisions that
decided whether an independent reproduction actually *matched the published truth*.

- **Detecting the donor from raw geometry.** The tool takes no atom indices: it has to find
  the phosphorus donor and its substituents from coordinates alone, using Cordero
  covalent-radius bonding, and reject anything that isn't a valid trivalent donor (on
  stderr, so a messy batch folder still completes). It isn't phosphorus-only: the path
  is [validated on nitrogen donors](docs/validation/NONP_DONOR.md) against morfeus.
- **Telling the kernel apart from the input geometry.** When the native descriptor sat below
  Kraken's published values (R² ≈ 0.86), the real question was whether *my kernel* was wrong
  or *my geometries* were. Resolved by changing one variable at a time — RDKit/MMFF →
  CREST/xTB → Kraken's own DFT structures — which localized the shortfall to geometry and
  proved the kernel exact ([Studies 002–004](docs/study_004/STUDY_004.md)).
- **Matching Kraken's coordinate conventions.** The reproduction lived or died on constants
  that are easy to miss: a virtual metal placed 2.28 Å from phosphorus (not 2.1), and the
  Sterimol coordination axis with a +0.40 Å Verloop *L* correction. These came from reading
  Kraken's *own* source, not from guessing distances.
- **The phosphine frame bug — found only at scale.** The buried-volume frame took a donor's
  three *nearest heavy atoms* instead of its *covalently bonded* neighbors, silently
  discarding bonded hydrogens on primary/secondary phosphines. Eleven test ligands could
  never trigger it; the full 1,541-ligand set did. Switching to covalent-radius bonding
  lifted R² from 0.9649 → 0.9852 **without discarding a single ligand**
  ([Study 004](docs/study_004/STUDY_004_SCALED.md)).
- **Hand-written SIMD that must match its fallback.** The inference hot loop uses an explicit
  unsafe AVX2 kernel guarded by *runtime* CPU-feature detection, with a portable scalar
  fallback for everything else. The engineering burden isn't speed — it's keeping two numeric
  paths and a block of `unsafe` code consistent, so a test asserts the AVX2 result agrees with
  the scalar fallback across thousands of random inputs (to f32 tolerance, since the two sum
  in different orders).
- **Reproducibility, honestly scoped.** Toolchains are checksum-pinned (CREST 2.12 / xTB
  6.4.0), caches are content-addressed, predictions are frozen by SHA-256, and quantum jobs
  resume from checkpoints under PID/start-time locks. The `.sigpack` binary format is
  deliberately *machine-native* (fast, not portable) — so it carries an **endian marker**
  that makes a mismatched file fail loudly instead of silently corrupting.

---

## Scientific studies

Ten studies build from a small published reaction family up to the full Kraken library,
two independent reaction models, and the numerical convergence of the descriptor itself.
Each writes a complete model card, frozen prediction hashes, raw comparisons, plots, and
machine-readable results under `docs/study_00N/`. **Passed and failed gates are both
retained.**

<details>
<summary><b>Expand the ten studies</b></summary>

| # | Study | What it shows | Full results |
|---|---|---|---|
| **001** | Ni-hDA enantioselectivity model | The published descriptor reproduces the selectivity; StericX's *compact native* descriptors deliberately do **not** — the ablation the project leads with. | [STUDY_001](docs/study_001/STUDY_001.md) |
| **002** | Coordination-aware buried-volume fidelity | The Rust voxel kernel equals `morfeus` exactly on identical geometries; cheap RDKit/MMFF conformers fall short — a failed gate kept visible. | [STUDY_002](docs/study_002/STUDY_002.md) |
| **003** | Quantum geometry & prospective validation | A checksum-pinned CREST/xTB backend, and a **frozen, falsifiable prediction** committed by SHA-256 before any measurement. | [STUDY_003](docs/study_003/STUDY_003.md) · [PREREGISTRATION](docs/study_003/PREREGISTRATION.md) |
| **004** | Reproducing Kraken's published descriptors on DFT geometries | Reproduces the published values on Kraken's own DFT geometries at library scale; a real kernel frame-bug found and fixed *without dropping a ligand*, and the residual fully characterized. | [STUDY_004](docs/study_004/STUDY_004.md) · [scaled](docs/study_004/STUDY_004_SCALED.md) · [residual](docs/study_004/STUDY_004_RESIDUAL.md) |
| **005** | Pyramidalization descriptors | Two more descriptors (`pyr_P`, `pyr_alpha`) reduced to closed forms and reproduced to machine precision. | [STUDY_005](docs/study_005/STUDY_005.md) |
| **006** | Localizing the residual to the coordination centre | A controlled test proving the small buried-volume residual is a coordination-centre convention artefact, not a kernel error. | [STUDY_006](docs/study_006/STUDY_006.md) |
| **007** | Independent second reaction model — Ni cross-coupling | Reproduces the Newman-Stonebraker classifier; the ligation cliff transfers *out-of-sample* and off Kraken's own geometry. | [STUDY_007](docs/study_007/STUDY_007.md) |
| **008** | Head-to-head speed benchmark vs `morfeus` | The same numbers, ~14× faster single-core — and the speedup holds at ~20× the scale. | [STUDY_008](docs/study_008/STUDY_008.md) · [scale check](docs/study_008_all_conformers/STUDY_008.md) |
| **009** | The other direction of the cliff — Pd cross-coupling | The *opposite* (bulky-active) regime reproduced, including datasets from other groups; honest about the reactions the paper itself flags as resistant. | [STUDY_009](docs/study_009/STUDY_009.md) |
| **010** | Grid convergence of the buried-volume integrator | Sweeps the integration grid coarse→fine to show the descriptor is *converged* at the default resolution — the earlier agreement isn't a grid-lucky artifact — and quantifies the voxel discretization floor. | [STUDY_010](docs/study_010/STUDY_010.md) |

A manuscript-style narrative of the reproduction studies is in [`docs/REPRODUCTION_REPORT.md`](docs/REPRODUCTION_REPORT.md).

</details>

---

## Documentation

- 📄 **[Manuscript write-up](docs/REPRODUCTION_REPORT.md)** — the full narrative.
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

### `db build` — precompute a ligand database

Featurizing is the expensive step; reading a table back is not. `db build` does the work once
and writes a CSV plus a `.manifest.json` recording the exact settings, counts, and a SHA-256
of the table, so a database is a reproducible artifact rather than an ad-hoc dump.

```bash
./target/release/stericx db build --source ligands/ --output my_ligands.csv
# Kraken cache layout: each directory is one ligand, its files are conformers
./target/release/stericx db build --source .stericx/kraken_dft_cache \
    --output db.csv --group-by-parent --label-from parent --extension sdf
```

`--group-by-parent` aggregates a ligand's conformers into one row; `max_delta_qvbur_min`
takes the **minimum** across them because that is Kraken's own `*_min` convention, so
aggregating any other way would silently redefine the descriptor. `--label-from parent` makes
the directory name the ligand label — for the Kraken cache that is the molecule id, so hits
come back as `723` rather than an opaque path. `--extension` restricts which formats are read,
which matters when a tree mirrors the same conformers in more than one format.

**A database ships with the repo**: [`data/ligand_db/kraken_phosphines.csv`](data/ligand_db/)
— all **1,541** Kraken phosphines aggregated from **31,611** DFT conformers, built in 5.3 s,
215 KB. These are StericX's *own* computed descriptors keyed by Kraken molecule id, not
Kraken's published values; the manifest records the settings that produced them.

### `search` — find sterically similar ligands, or filter by constraints

```bash
# nearest neighbours to a query ligand
./target/release/stericx search --similar-to query.sdf \
    --database data/ligand_db/kraken_phosphines.csv

# a pure constraint query — no query ligand needed
./target/release/stericx search --database data/ligand_db/kraken_phosphines.csv \
    --vbur 30:35 --b5-max 8

# "find me a less bulky ligand of similar shape"
./target/release/stericx search --similar-to query.sdf --database db.csv --less-bulky
```

Searching all 1,541 ligands takes **~5 ms** against a prebuilt database, versus ~5 s to
re-featurize them; `--database` also accepts a plain directory and will featurize on the fly.

Similarity is Euclidean distance in **standardized** descriptor space — every descriptor is
z-scored against the database first, so an Å-scale `L` and a percent-scale `%Vbur` contribute
comparably instead of whichever carries the biggest units dominating. The default space is the
shape envelope (`L`, `B1`, `B5`), `%Vbur`, the quadrant asymmetry `max_delta_qvbur`, and the
donor pyramidalization `pyr_P`; choose your own with `--features l,b5,vbur`. `buried_volume`
is deliberately *not* in the default set — it is `percent_buried_volume` rescaled by a
constant, so including both would double-weight the same quantity.

Constraints narrow the field **before** ranking. Ranges take `LOW:HIGH`
(`--vbur 30:35`, `--l`, `--b1`, `--b5`) and each has `--…-min` / `--…-max` inclusive bounds
(`--b5-max 8`). `--filter` remains available for anything else, accepting `name=LOW..HIGH`,
`name<V`, `name<=V`, `name>V`, `name>=V` over any descriptor. `--less-bulky` /
`--more-bulky` are shorthand for a `%Vbur` bound relative to the query. With no
`--similar-to` the result is ordered by the first constrained descriptor (or `--sort-by`),
and the distance column reads `—` rather than inventing a similarity to nothing.

Descriptors that are constant across the database are reported and excluded, so a ranking is
never silently computed on fewer axes than requested.

> **Two honest caveats.** This ranks *steric similarity, not reactivity* — two ligands close
> in this space occupy space alike; whether they behave alike in a given reaction is an
> experimental question, and `screen` is the command that brings a reaction model to bear.
> And compare like with like: the shipped database holds per-ligand **conformer-ensemble**
> values, so querying with a single conformer compares an ensemble average against one
> geometry. The query's own descriptors are printed so the difference is visible.

### `compare` — put ligands side by side

```bash
./target/release/stericx compare ligand_a.sdf ligand_b.sdf ligand_c.sdf \
    --database data/ligand_db/kraken_phosphines.csv
```

```text
descriptor                   31676     46112     33790    spread        σ
sterimol_l                   6.991     6.883     6.959     0.107     0.07
sterimol_b1                  4.956     3.568     2.841     2.115     2.61
percent_buried_volume       58.612    41.018    30.562    28.050     3.03
```

Every descriptor for every ligand, the spread across them, and — with `--database` — that
spread in **library standard deviations** plus a standardized pairwise distance. The σ column
is what makes a raw number mean something: above, these three ligands are the same *length*
(0.07 σ) while differing enormously in minimum width (2.61 σ) and bulk (3.03 σ). Without a
database the raw differences still print; they just cannot be placed on a scale, and the
output says so.

### `screen` — rank a library with a fitted reaction model

Where `search` ranks by *shape*, `screen` ranks by *predicted performance* under a model
fitted by `stericx fit`, and reports what the prediction is worth.

```bash
./target/release/stericx screen model.json ligand_library/
./target/release/stericx screen model.json reactions.csv --top 20 --inside-domain-only
```

Each ligand gets a predicted ΔΔG‡, the corresponding ee at `--temperature` (signed by the
ΔΔG‡ convention, so a negative value means the same excess of the opposite enantiomer), a
conservative uncertainty band, and an applicability-domain verdict. Ligands outside the
range the model was trained on are listed separately with **how far** outside they fall, as
a fraction of the training range width — `--inside-domain-only` drops them entirely.

**The model decides what the library must supply.** StericX's regression space mixes
geometry (`L`, `B1`, `B5`) with donor electronics (`nbo_charge`, `ir_frequency`) and their
interactions. `screen` reads the fitted weights, works out which inputs actually carry a
nonzero coefficient, and refuses to run when the library cannot provide one — it will not
invent a donor charge to make a number appear:

```text
error: model `mechanistically_constrained_ols` uses B5_x_nbo_charge but the library does
not provide nbo_charge. StericX will not guess a missing input. …
```

So a model whose selected features are geometry-only screens a bare directory of
structures, while a model carrying an electronic term needs those values as CSV columns. A
reaction CSV supplies both at once: `screen` reads `NBO_Charge` / `IR_Frequency` from the
row and featurizes the geometry named by `Ligand_XYZ_Path` to fill the Sterimol terms.

**Uncertainty and out-of-domain detection.** `stericx fit` records the training-set geometry
— `(X'X)⁻¹` in the standardized design frame, `n`, `p`, and the residual standard error `s` —
so `screen` can answer *how much a prediction is worth*, not just what it is. Two
independent signals are reported:

| Signal | What it measures |
|---|---|
| **Leverage** `h = x'(X'X)⁻¹x` vs `h* = 3p/n` | how far the ligand sits from the centre of the training design, in the metric the fit itself defines |
| **Per-feature range check** | whether each selected descriptor lies inside its training min/max (a 1-D box) |

They can disagree, and that is the point: a ligand can sit inside every 1-D range yet still
be far from the training cloud. Both are reported and the worse one governs a graded verdict
— `reliable`, `caution:high_leverage`, `caution:outside_range`, or
`do_not_trust:extrapolation`. Ligands in the last grade get called out explicitly:

```text
1 ligand(s) are outside the training range AND above the warning leverage.
For these the model is extrapolating and its prediction should not be trusted:
  wild_extrap — leverage 3.396 = 5.7x h*
```

The 95 % band is a real Student-t **prediction interval**, `ŷ ± t(0.975, n−p)·s·√(1+h)`, so it
widens automatically with leverage — a distant ligand is reported with an honestly wider
error bar rather than a falsely precise one. On the Ni-hDA fit that is ±1.49 kcal/mol at the
training centroid and ±2.97 kcal/mol for a ligand at 5.7× the warning leverage. `h* = 3p/n`
is the same criterion the project's own [Study 003
pre-registration](docs/study_003/PREREGISTRATION.md) applies, and the Rust implementation
reproduces its `h* = 0.60` exactly.

Models fitted before the geometry was recorded still load: leverage and the prediction
interval report as unavailable, the verdict degrades to `range_only:inside` /
`range_only:outside`, the weaker bootstrap band is printed marked `~[…]` so it can never be
mistaken for a prediction interval, and the output says to refit to enable the leverage
check. Scoping the claim down beats faking it.

The `coefficient_band` remains available as the weaker signal: it propagates the bootstrap
coefficient intervals by interval arithmetic, ignoring coefficient correlation. A band that
spans zero tells you the model cannot commit to a direction for that ligand at all.

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

Add `--portable-model <path>` to also write a schema-2 portable model: the same document plus
the response definition, feature transformations, training-data digests, and creation
metadata, so it can be scored elsewhere without the training data or the fitting code.

```bash
./target/release/stericx fit --data data/reactions.sigpack --metadata data/reactions_raw.csv \
  --output docs/study_001/stericx_model.json --predictions docs/study_001/stericx_frozen_predictions.csv \
  --portable-model docs/study_001/stericx_portable_model.json \
  --reaction-family "Ni-catalyzed homo-Diels-Alder" --catalyst-metal Ni --response-temp-k 298.15
```

Schema 2 is a strict superset of schema 1, so `--output` stays byte-identical and existing
readers keep working. Chemistry context that is not supplied is recorded as `null` and
reported as `portable_model_missing_provenance` rather than guessed. Spec:
[`docs/MODEL_FORMAT.md`](docs/MODEL_FORMAT.md).

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

studies/          Study drivers (Python → docs/study_00N/), study_001 … study_010
scripts/          Support utilities: prepare_data, prepare_quantum_data, validate_stericx,
                  stericx_quantum, freeze_prospective_deck, preregister_prediction,
                  validate_nonp_donor
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

**v0.1.0** shipped the science: the descriptor kernels, full-set Kraken reproduction, the
frame-bug fix and residual characterization, two reaction models (Ni + Pd cross-coupling),
pyramidalization, the speed benchmark, and the Zenodo release.

**v0.2.0** makes it useful for *choosing* between ligands rather than only measuring them:
a precomputed 1,541-ligand descriptor database, similarity search, constraint queries,
side-by-side `compare`, and `screen` — reaction-model ranking with leverage-based
applicability-domain warnings and Student-t prediction intervals.

The one open scientific item remains the **frozen prospective 10-candidate deck**, which
stays target-free until experimental measurement exists. A measured outcome would seed a
future release under the existing concept DOI — and never a refit, for pre-registration
integrity.

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
> inference*, version 0.2.0, 2026. Zenodo.
> DOI: [10.5281/zenodo.21726666](https://doi.org/10.5281/zenodo.21726666).

The concept DOI [10.5281/zenodo.21726666](https://doi.org/10.5281/zenodo.21726666) resolves
to the latest version. Version DOIs: v0.2.0 =
[10.5281/zenodo.21985632](https://doi.org/10.5281/zenodo.21985632), v0.1.0 =
[10.5281/zenodo.21726667](https://doi.org/10.5281/zenodo.21726667).

---

## License

Licensed under either of [Apache-2.0](LICENSE-APACHE) or [MIT](LICENSE-MIT) at your option
(SPDX: `MIT OR Apache-2.0`). Unless you state otherwise, any contribution intentionally
submitted for inclusion shall be dual licensed as above, without additional terms.
