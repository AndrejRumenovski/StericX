# Earlier pre-policy benchmark: retained measurements, withdrawn admission claim

This is the complete **earlier** campaign requested for publication: 240 launches
(volume and pyramidalization), eight separate memory runs, and its reproduction
smoke test. The measurements and raw outputs are retained below.

**The scientific-admission claim in the original report is withdrawn.** Its policy
was not frozen before classification. Statements below that a gate “passed” describe
the original regional-count/rounding checks; they do not mean the newer prospective
scientific policy admitted this campaign. In particular, observed f32-rounded
pyramidalization agreement is not an independently established error bound.

The [new descriptor-specific policy](EQUIVALENCE_POLICY.md) and its
[freeze receipt](policy_freeze.json) define the later work. That campaign is outside
this published snapshot. No new timing is substituted for an earlier measurement.
The original report follows as a historical record, not a current validated
speedup claim. Its labels and numbers must be read with this status correction.

---

# Current StericX versus morfeus: frozen end-to-end descriptor benchmarks

## 1. Executive result

On **10,000 files cyclically repeating 56 real conformers from 11 ligands**, the
current unmodified StericX buried-volume API was **25.22× faster
on one CPU core** and **20.25× faster at equal six-core parallelism**
than morfeus 0.8.0, using end-to-end command-line API drivers on this host.
These are medians of four fresh paired `morfeus wall / StericX wall` ratios.
No historical speed ratio enters the calculation.

The headline covers **nine buried-volume outputs, including all three donor-plane
calculations**, with exactly equal regional point populations and occupied counts.
Every floating-point residual is reconstructed exactly from those independent
counts using the native documented f32 rounding operations. No error tolerance
was invented or inferred from a historical maximum.

Separate pyramidalization timing gives **7.34× at one core**
and **4.51× at six cores** on the same 10,000-file repetition.
Sterimol and combined-descriptor performance are **excluded**: finite angular grids,
transverse phase and B5 algorithms differ. These results do not establish the
speed of the shipping `stericx descriptors` command, which also calculates Sterimol.

## 2. Hardware and execution environment

- AMD Ryzen 5 5600G; 6 physical cores, 12 logical CPUs; 32,179,068 KiB usable RAM
  (30.69 GiB), Ubuntu 26.04.1 LTS, Linux 7.0.0-31-generic x86-64.
- Benchmark date: 2026-09-22 UTC. Exact launch timestamps are retained.
- CPU affinity: 1 worker → CPU 2; 2 → CPUs 2,3; 4 → CPUs 0–3; 6 → CPUs 0–5.
  Each six-core mask uses one logical CPU from each physical core. The same mask
  applies to both tools, including Python parent/workers and native threads.
- Ordinary desktop applications remained active. No competing benchmark builds,
  scientific reference campaigns, or memory sampling ran during timed pairs.
  This is a desktop measurement, not an isolated dedicated host.
- [Full machine captures](environment/), including topology, CPU flags, OS/kernel,
  RAM, Python/Rust versions and NumPy build/runtime configuration.

## 3. Frozen software

StericX commit: `0143f153f974cb66c723e3ad21e42cc51e211d58`; Rust 1.97.0
(`2d8144b78`, LLVM 22.1.6). Production source hashes, exact command and environment
are in [freeze.json](freeze.json). The pre-existing unrelated `docs/media/` files
were left untouched. No production source was changed or optimized for this study.

Build: `cargo build --locked --offline --release --bin stericx --target-dir …`.
Release uses opt-level 3, thin LTO, one codegen unit, no optional Cargo features,
no RUSTFLAGS, and rustc's default x86-64 target CPU. Enabled target features:
`fxsr`, `sse`, `sse2`; host capabilities are recorded separately.

- Shipping executable SHA-256: `b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6`.
- Timed volume driver SHA-256: `09de3b179a9d9ee816c73f405ab31947733a381b8ab6e6c6e55ba35b8189eafb`.
- The timed driver links the unchanged current public API, uses Rayon `par_iter`
  for independent files, and emits the complete nine-field result for every file.
  Its dependencies are locked in [adapter/Cargo.lock](adapter/Cargo.lock).
  [Source snapshot](build/source.tar.gz), production/driver binaries and build logs
  are retained. The untimed grid observer is a byte-identical source copy plus
  a visibility probe; it supplies **no timing samples**.

morfeus: **0.8.0**, the current published package at acquisition, installed into a
fresh Python **3.12.13** environment with `uv pip install --index-url
https://pypi.org/simple morfeus-ml==0.8.0`. Dependencies: NumPy 2.5.3, SciPy 1.18.1,
fire 0.7.1, packaging 26.3, termcolor 3.3.0. Exact requirements, source hashes,
installation output and [PyPI response](environment/pypi-morfeus.json) are retained.
The package identity is `morfeus-ml`, not the unrelated `morfeus` distribution.
[Reference package](https://pypi.org/project/morfeus-ml/0.8.0/).

`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `BLIS_NUM_THREADS`,
`VECLIB_MAXIMUM_THREADS` and `NUMEXPR_NUM_THREADS` are all 1. StericX uses
`RAYON_NUM_THREADS=1/2/4/6`. Ordinary morfeus usage is a one-process Python loop;
the parallel baseline is straightforward `multiprocessing` **spawn** + `Pool.map`,
with its default chunk sizing. Worker startup, IPC, parent and resource tracker
costs are included. No geometry/descriptor memoization or grid cache was added.

## 4. Identical input corpus

[input_manifest.json](input_manifest.json) contains all original 10,000 filenames,
SHA-256s, atom counts and the mapping to the 56 unique conformers, as well as their
ligand IDs, donor identities and configuration. [inputs/](inputs/) preserves the
original 56 XYZ files byte-for-byte. [geometry_configuration.json](geometry_configuration.json)
adds donor-neighbor identities, represented native centers and radii.

The source is the admitted corrected StericX profiling workload, whose manifest
hash is `a1041cc07fc5583c33d2a7cf584dbc638866c64bdd0b64e3720ef83cc85e1e6a`.
All original repeated-file hashes were verified before copying. Both tools read
the **same paths and bytes** on every run. The 1,000-file corpus is the first 1,000
files of that repeated workload; neither it nor the 10,000-file corpus represents
independent new chemistries. Repeated files are regenerated from the preserved
56 files by `benchmark.py prepare`; no hidden original workspace is needed.

## 5. Scientific conventions and equivalent work

The timed buried-volume output is: occupied volume, %Vbur, quadrant minimum and
maximum, maximum adjacent-quadrant difference, octant minimum and maximum, and
near/far occupied volumes. Both tools calculate **three complete orientation
grids per molecule**. Totals and near/far volumes are averaged over those three
planes; extrema and adjacent differences are reduced over all three.

Donor: sole phosphorus atom. Donor neighbors: the unchanged 1.3-times-summed
Cordero covalent-radius rule, including bound hydrogens. All 56 inferred donor
identities/neighbors are checked against the native observations. Virtual metal:
2.28 Å along the negative sum of three unit donor–neighbor directions. Sphere:
3.5 Å; density: 0.01 Å³; explicit native Bondi-style radii × 1.17; hydrogen spheres
excluded from volume, with donor-bound H retained for topology/frame construction.
All corpus elements (C,H,N,O,P,S) have matching reference radii.

The center→donor vector maps to negative Z; each donor substituent defines positive
X in turn. The Python driver independently constructs this frame and presents a
zero dummy center to the normal morfeus `BuriedVolume` API, following the original
audit's workaround for the reference constructor's translated-center aliasing.
No morfeus source is patched. Octant mapping is `[0,1,2,3,7,6,5,4]`.
The default grid has 32 positions per cube edge, 15,408 retained sphere points,
1,926 points per octant, and no zero coordinate planes.

Both drivers parse XYZ, infer topology/center, compute descriptors, retain the
full ordered result array, and write ordinary pretty JSON. The shipping StericX
descriptor CLI cannot select volume alone; using these small public-API CLI
drivers isolates an equivalent scientific workload while retaining end-to-end I/O.
This is **not an in-memory kernel benchmark**, and is not presented as a timing of
the unmodified shipping CLI. The shipping CLI's volume results were checked
bit-for-bit against the linked native API before timing.

Native volume uses f32 represented input coordinates and its documented f32
lattice/volume arithmetic. The primary morfeus driver uses the same XYZ source
coordinates parsed as float64 and independently infers the center. An additional
controlled comparison supplies identical f32-represented coordinates, center and
radii; it is diagnostic, not the timed volume implementation.

Pyramidalization uses the same donor and three neighbors, Radhakrishnan P and
mean signed alpha. Its separate comparison rounds input coordinates to the same
f32 representation before morfeus's float64 calculation. Neither tool's formula
is modified. Native P/alpha exactly match morfeus rounded to f32 in all 56 cases.
The native/reference public APIs' internal intermediate work is left intact.

## 6. Scientific agreement before timing

The audit explicitly says that its observed error maxima are not universal
acceptance tolerances. Accordingly the scientific gate uses **zero-tolerance
discrete equality**, not an invented absolute/relative threshold:

1. All 168 native/reference orientation grids have exactly matching octant
   populations and occupied counts (1,344 occupied-count comparisons).
2. Independent morfeus counts reproduce every native per-plane volume and
   all nine aggregate outputs **bit-for-bit** through the documented native f32
   normalization/reduction operations. All residuals below are thus fully
   accounted for by rounding on this fixed corpus.
3. The actual timing drivers pass all 56 geometries at 1/2/4/6 workers. Every
   output in every timed repeated batch is also checked exactly against its
   tool's complete preflight record. No failed or disagreeing molecule is dropped.

See [scientific_gate.json](scientific_gate.json), [exact integer/rounding evidence](raw/exact_volume_evidence.json),
[preflight.json](preflight.json), and the separate [pyramidalization gate](pyramidalization/scientific_gate.json).
This applies the independently documented [volume equations and grid regions](../../scientific_accuracy_audit/geometry/METHODS.md)
and [corrected three-plane reductions](../../scientific_remediation/geometry/RESULTS.md)
without widening any audit limit. It does not claim continuum accuracy or universal
agreement on other grids, molecules, elements or orientations.

All fields below were compared before timing on the complete 56-conformer corpus,
including the excluded Sterimol fields. Volumes are Å³, percent volume is percentage
points, Sterimol is Å, P is dimensionless and alpha is degrees. R² is against the
identity line; regression is native versus reference. A constant reference has
undefined R²/slope/intercept.

| Descriptor | N | MAE | RMSE | Max absolute | Median absolute | R² (identity) | Slope | Intercept |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sterimol_l | 56 | 5.56959e-07 | 6.60336e-07 | 1.69531e-06 | 5.43283e-07 | 1 | 1 | 1.2002e-06 |
| sterimol_b1 | 56 | 0.00690073 | 0.0089705 | 0.0180639 | 0.00553108 | 0.999762 | 1.00287 | -0.00413621 |
| sterimol_b5 | 56 | 9.51556e-07 | 1.24355e-06 | 2.77774e-06 | 6.79668e-07 | 1 | 1 | 1.94109e-06 |
| pyr_p | 56 | 2.14312e-08 | 2.75563e-08 | 7.5722e-08 | 1.85669e-08 | 1 | 1 | 2.38837e-07 |
| pyr_alpha | 56 | 2.36727e-06 | 2.97751e-06 | 6.65029e-06 | 2.05702e-06 | 1 | 1 | 3.06952e-07 |
| buried_volume | 56 | 1.95246e-06 | 2.58494e-06 | 6.13951e-06 | 1.57181e-06 | 1 | 1 | 2.64391e-06 |
| percent_buried_volume | 56 | 1.38736e-06 | 1.77461e-06 | 4.26496e-06 | 9.76452e-07 | 1 | 1 | -2.36786e-07 |
| qvbur_min | 56 | 3.2366e-07 | 4.41328e-07 | 1.28376e-06 | 2.17816e-07 | 1 | 1 | 9.96205e-07 |
| qvbur_max | 56 | 7.46135e-07 | 9.45649e-07 | 2.4031e-06 | 6.3236e-07 | 1 | 1 | 2.34266e-07 |
| max_delta_qvbur | 56 | 7.21165e-07 | 9.52181e-07 | 2.25736e-06 | 5.39841e-07 | 1 | 1 | 1.46283e-07 |
| ovbur_min | 56 | 0 | 0 | 0 | 0 | undefined | undefined | undefined |
| ovbur_max | 56 | 7.1185e-07 | 8.20349e-07 | 1.79641e-06 | 7.77577e-07 | 1 | 1 | -3.03428e-07 |
| near_vbur | 56 | 2.0633e-06 | 2.73372e-06 | 8.29924e-06 | 1.39934e-06 | 1 | 1 | 2.73637e-06 |
| far_vbur | 56 | 2.10756e-07 | 4.32654e-07 | 1.57598e-06 | 3.11396e-08 | 1 | 1 | -2.10665e-08 |

[Every native/reference value, absolute and relative error](scientific_comparisons.csv)
and [all three comparison-lane statistics](scientific_summary.json) are retained.
Relative error at zero reference is undefined. The independent analytical
[Sterimol diagnosis](sterimol_diagnosis.json) retains every conformer: native B1
uses a 360-direction shortest-arc frame; morfeus uses a 3,600-direction endpoint-inclusive
Kabsch frame. A 361-direction diagnostic matches one-degree spacing but still has
a different phase. Native B5 is a direct radial maximum; morfeus samples angular
support. Those are numerical-method differences, so **no Sterimol or combined
speedup is admitted**. No phase, radius, geometry or tolerance was tuned to pass.

## 7. Benchmark methodology

Warm filesystem cache; one complete warmup per tool/configuration, followed by
four measured pairs in **AB/BA/AB/BA** order. A=StericX, B=morfeus. Matrix:
56, 1,000 and 10,000 files × 1/2/4/6 workers, separately for volume and pyramidalization.
Each experiment retains **120 launches: 24 warmups and 96 measured runs**.
No outlier is removed. Both tools use the same physical-core affinity mask.

The frozen C launcher starts a monotonic wall clock before fork/exec and stops
after `wait4`, including process/interpreter/worker startup, imports, argument
handling, input parsing, descriptor calculation, output serialization/writes and
shutdown. The affinity launcher and opening stdout/stderr files are outside this
interval for both tools. Outputs are real files, not `/dev/null`; compression and
validation happen after timing. No fsync or cold-cache claim is made.

Record wall time, user/system CPU, wait4 peak RSS, faults and context switches.
CPU time includes reaped worker descendants; verify scaling in the raw samples.
Throughput is file count divided by median wall seconds. Paired speedup is
`morfeus wall / StericX wall` for each pair, summarized by median/MAD/full range.
Ratio of separate wall medians is also supplied. MAD is not a confidence interval.
The measured separation exceeds noise: even the slowest paired advantage remains
above 1. All startup-sensitive small batches remain reported.

## 8. Raw timing summary: buried volume

| Files | Workers | Tool | Wall median ± MAD, s | Wall range, s | Files/s | CPU median, s | wait4 RSS median, MiB* |
|---:|---:|---|---:|---:|---:|---:|---:|
| 56 | 1 | stericx | 0.0311 ± 0.0001 | 0.0310–0.0316 | 1,798.2 | 0.0307 | 2.90 |
| 56 | 1 | morfeus | 1.0484 ± 0.0161 | 1.0059–1.0685 | 53.4 | 1.0479 | 71.89 |
| 56 | 2 | stericx | 0.0182 ± 0.0001 | 0.0180–0.0215 | 3,082.1 | 0.0312 | 3.19 |
| 56 | 2 | morfeus | 0.9045 ± 0.0081 | 0.8889–0.9211 | 61.9 | 1.5321 | 71.81 |
| 56 | 4 | stericx | 0.0122 ± 0.0005 | 0.0115–0.0161 | 4,573.1 | 0.0346 | 3.99 |
| 56 | 4 | morfeus | 0.7600 ± 0.0080 | 0.7419–0.7701 | 73.7 | 2.0708 | 71.88 |
| 56 | 6 | stericx | 0.0112 ± 0.0002 | 0.0097–0.0114 | 5,000.7 | 0.0363 | 4.49 |
| 56 | 6 | morfeus | 0.6801 ± 0.0053 | 0.6747–0.7070 | 82.3 | 2.7041 | 71.80 |
| 1,000 | 1 | stericx | 0.5407 ± 0.0044 | 0.5362–0.5730 | 1,849.5 | 0.5305 | 4.25 |
| 1,000 | 1 | morfeus | 14.8484 ± 0.2340 | 14.1567–15.3150 | 67.3 | 14.8351 | 74.22 |
| 1,000 | 2 | stericx | 0.3091 ± 0.0018 | 0.3035–0.3113 | 3,234.8 | 0.5434 | 4.62 |
| 1,000 | 2 | morfeus | 7.9894 ± 0.1099 | 7.7467–8.2088 | 125.2 | 15.5510 | 72.61 |
| 1,000 | 4 | stericx | 0.1869 ± 0.0011 | 0.1845–0.1884 | 5,351.9 | 0.5571 | 4.61 |
| 1,000 | 4 | morfeus | 4.6308 ± 0.0189 | 4.5727–4.6584 | 215.9 | 17.1610 | 72.44 |
| 1,000 | 6 | stericx | 0.1451 ± 0.0012 | 0.1432–0.1495 | 6,893.2 | 0.5726 | 5.06 |
| 1,000 | 6 | morfeus | 3.4359 ± 0.0389 | 3.3812–3.5240 | 291.0 | 18.2969 | 72.39 |
| 10,000 | 1 | stericx | 5.5475 ± 0.0341 | 5.4127–5.6150 | 1,802.6 | 5.4534 | 16.50 |
| 10,000 | 1 | morfeus | 139.9611 ± 0.7571 | 136.5634–140.8717 | 71.4 | 139.8861 | 97.15 |
| 10,000 | 2 | stericx | 3.1500 ± 0.0608 | 3.0741–3.5851 | 3,174.6 | 5.4982 | 16.87 |
| 10,000 | 2 | morfeus | 78.6823 ± 0.6834 | 75.8768–79.5443 | 127.1 | 152.4265 | 97.50 |
| 10,000 | 4 | stericx | 1.8362 ± 0.0237 | 1.8110–1.8900 | 5,445.9 | 5.6372 | 17.75 |
| 10,000 | 4 | morfeus | 41.7912 ± 0.1956 | 41.4109–42.6943 | 239.3 | 160.0830 | 96.82 |
| 10,000 | 6 | stericx | 1.4643 ± 0.0226 | 1.4006–1.5016 | 6,829.2 | 5.9802 | 18.40 |
| 10,000 | 6 | morfeus | 29.7052 ± 0.0728 | 29.5293–29.7914 | 336.6 | 169.6964 | 96.49 |

*For multiprocessing, wait4 RSS is the largest individual process high-water mark, **not summed worker memory**. Full CPU/user/system/RSS MAD, range and samples are in JSON.

## 9. Single-core implementation/runtime comparison

For 10,000 files, StericX median wall time is **5.5475 s**;
morfeus is **139.9611 s**. Paired speedup is
**25.2155×**, MAD 0.3564,
range 24.6136–25.7466×. Both are pinned to CPU 2, with one
native thread or one Python process and numerical-library threads limited to one.
This includes driver and runtime costs, so it is not a pure-language or algorithm-only attribution.

## 10. Equal-core parallel throughput

For 10,000 files and six physical cores, native six-thread median wall is
**1.4643 s**; morfeus six-worker median is
**29.7052 s**. Paired speedup is
**20.2467×**, MAD 0.2664,
range 19.7434–21.2703×. This is an equal-hardware batch comparison,
including Python multiprocessing costs; it is not attributed wholly to Rust.

| Files | Equal workers/cores | Paired speedup median ± MAD | Full paired range | Ratio of wall medians |
|---:|---:|---:|---:|---:|
| 56 | 1 | 33.7502 ± 0.4683× | 31.8612–34.2597× | 33.6645× |
| 56 | 2 | 49.3377 ± 0.7451× | 42.8475–50.1591× | 49.7819× |
| 56 | 4 | 62.7502 ± 1.6402× | 46.8784–64.5503× | 62.0634× |
| 56 | 6 | 62.1664 ± 1.7037× | 58.9677–69.2607× | 60.7354× |
| 1,000 | 1 | 27.4643 ± 0.6600× | 24.7061–28.5616× | 27.4625× |
| 1,000 | 2 | 25.6923 ± 0.0994× | 25.5207–26.6818× | 25.8443× |
| 1,000 | 4 | 24.7817 ± 0.2092× | 24.3728–25.1512× | 24.7833× |
| 1,000 | 6 | 23.5770 ± 0.5547× | 22.8294–24.6093× | 23.6843× |
| 10,000 | 1 | 25.2155 ± 0.3564× | 24.6136–25.7466× | 25.2295× |
| 10,000 | 2 | 24.4532 ± 0.6579× | 22.1872–25.7595× | 24.9788× |
| 10,000 | 4 | 22.7820 ± 0.5206× | 21.9110–23.5349× | 22.7591× |
| 10,000 | 6 | 20.2467 ± 0.2664× | 19.7434–21.2703× | 20.2864× |

## 11. Throughput

The following table is **10,000 repeated input files**, not 10,000 independent
chemistries. The last column uses ordinary single-worker morfeus as a practical
baseline and must not be interpreted as per-core efficiency for multi-core rows.

| Tool | Workers | Wall median, s | Files/s | Relative to ordinary morfeus 1 worker |
|---|---:|---:|---:|---:|
| stericx | 1 | 5.5475 | 1,802.6 | 25.229× |
| morfeus | 1 | 139.9611 | 71.4 | 1.000× |
| stericx | 2 | 3.1500 | 3,174.6 | 44.433× |
| morfeus | 2 | 78.6823 | 127.1 | 1.779× |
| stericx | 4 | 1.8362 | 5,445.9 | 76.221× |
| morfeus | 4 | 41.7912 | 239.3 | 3.349× |
| stericx | 6 | 1.4643 | 6,829.2 | 95.583× |
| morfeus | 6 | 29.7052 | 336.6 | 4.712× |

StericX's six-thread versus ordinary morfeus one-worker wall-median ratio is
95.583×. This is explicitly an unequal-core practical
throughput comparison and is excluded from the equal-core README headline.

## 12. Memory

Timed wait4 peaks are in the raw summary above. For a single-process tool they
provide the kernel high-water RSS; for multiprocessing they **do not sum the tree**.
A separate 10,000-file execution per configuration samples concurrent parent,
worker and resource-tracker RSS, targeted every 10 ms. Sampling overhead does
not affect the headline timings. Peaks below are the largest **simultaneous**
sum; individual process peaks are never added together.

| Tool | Threads/workers | Sampled peak simultaneous sum RSS, MiB | Largest sample spacing, ms |
|---|---:|---:|---:|
| stericx | 1 | 16.71 | 11.64 |
| morfeus | 1 | 97.36 | 11.22 |
| stericx | 2 | 17.13 | 11.65 |
| morfeus | 2 | 260.31 | 16.00 |
| stericx | 4 | 17.98 | 10.49 |
| morfeus | 4 | 407.71 | 11.88 |
| stericx | 6 | 18.92 | 11.35 |
| morfeus | 6 | 551.99 | 12.84 |

Summed RSS counts shared pages in each process, so it is not unique physical
memory or PSS. Sampling can miss brief peaks; the actual maximum spacing is
shown. These are separate sampled runs, not repeated RSS distributions. Full
process/PID/RSS time series, peak snapshots, output checks and commands are in
[memory/](memory/) and [memory_results.json](memory_results.json).

## 13. Separate pyramidalization results and limitations

| Files | Workers | Tool | Wall median ± MAD, s | Wall range, s | Files/s | CPU median, s | wait4 RSS median, MiB* |
|---:|---:|---|---:|---:|---:|---:|---:|
| 56 | 1 | stericx | 0.0034 ± 0.0001 | 0.0033–0.0036 | 16,466.3 | 0.0031 | 3.23 |
| 56 | 1 | morfeus | 0.2623 ± 0.0011 | 0.2605–0.2643 | 213.5 | 0.2619 | 68.06 |
| 56 | 2 | stericx | 0.0033 ± 0.0002 | 0.0029–0.0035 | 17,026.3 | 0.0037 | 3.12 |
| 56 | 2 | morfeus | 0.5086 ± 0.0067 | 0.4974–0.5339 | 110.1 | 0.7416 | 68.07 |
| 56 | 4 | stericx | 0.0031 ± 0.0001 | 0.0027–0.0032 | 18,169.7 | 0.0043 | 3.18 |
| 56 | 4 | morfeus | 0.4997 ± 0.0070 | 0.4893–0.5124 | 112.1 | 1.2181 | 68.00 |
| 56 | 6 | stericx | 0.0032 ± 0.0001 | 0.0030–0.0034 | 17,694.6 | 0.0049 | 3.05 |
| 56 | 6 | morfeus | 0.5089 ± 0.0014 | 0.5068–0.5340 | 110.0 | 1.7264 | 67.83 |
| 1,000 | 1 | stericx | 0.0449 ± 0.0005 | 0.0443–0.0594 | 22,261.9 | 0.0399 | 4.34 |
| 1,000 | 1 | morfeus | 0.5530 ± 0.0014 | 0.5515–0.5668 | 1,808.3 | 0.5506 | 70.17 |
| 1,000 | 2 | stericx | 0.0360 ± 0.0010 | 0.0346–0.0405 | 27,779.2 | 0.0421 | 4.29 |
| 1,000 | 2 | morfeus | 0.6698 ± 0.0035 | 0.6640–0.6817 | 1,493.0 | 1.0478 | 71.01 |
| 1,000 | 4 | stericx | 0.0326 ± 0.0013 | 0.0314–0.0350 | 30,647.4 | 0.0455 | 4.31 |
| 1,000 | 4 | morfeus | 0.6278 ± 0.0026 | 0.6241–0.6577 | 1,592.8 | 1.5509 | 70.75 |
| 1,000 | 6 | stericx | 0.0300 ± 0.0007 | 0.0288–0.0309 | 33,371.2 | 0.0454 | 4.19 |
| 1,000 | 6 | morfeus | 0.6391 ± 0.0037 | 0.6259–0.6441 | 1,564.8 | 2.1288 | 70.20 |
| 10,000 | 1 | stericx | 0.4498 ± 0.0018 | 0.4363–0.4525 | 22,230.2 | 0.3957 | 14.68 |
| 10,000 | 1 | morfeus | 3.3011 ± 0.0126 | 3.2681–3.3259 | 3,029.3 | 3.2639 | 90.40 |
| 10,000 | 2 | stericx | 0.3420 ± 0.0019 | 0.3388–0.3462 | 29,242.1 | 0.3945 | 14.63 |
| 10,000 | 2 | morfeus | 2.1858 ± 0.0504 | 2.1330–2.2656 | 4,574.9 | 3.9004 | 93.25 |
| 10,000 | 4 | stericx | 0.3018 ± 0.0169 | 0.2648–0.3320 | 33,136.9 | 0.4221 | 14.29 |
| 10,000 | 4 | morfeus | 1.4920 ± 0.0201 | 1.4545–1.5142 | 6,702.3 | 4.5601 | 92.79 |
| 10,000 | 6 | stericx | 0.2918 ± 0.0108 | 0.2798–0.3044 | 34,268.3 | 0.4301 | 14.32 |
| 10,000 | 6 | morfeus | 1.3102 ± 0.0017 | 1.3075–1.3242 | 7,632.3 | 5.2368 | 92.34 |

*For multiprocessing, wait4 RSS is the largest individual process high-water mark, **not summed worker memory**. Full CPU/user/system/RSS MAD, range and samples are in JSON.

| Files | Equal workers/cores | Paired speedup median ± MAD | Full paired range | Ratio of wall medians |
|---:|---:|---:|---:|---:|
| 56 | 1 | 76.9806 ± 1.8449× | 73.5972–79.0042× | 77.1128× |
| 56 | 2 | 156.1467 ± 7.4278× | 146.3115–175.4936× | 154.6308× |
| 56 | 4 | 163.6735 ± 4.7062× | 156.9299–180.6279× | 162.1177× |
| 56 | 6 | 164.6520 ± 2.8807× | 148.9065–170.3159× | 160.8039× |
| 1,000 | 1 | 12.4497 ± 0.0690× | 9.2829–12.5201× | 12.3109× |
| 1,000 | 2 | 18.5120 ± 0.5487× | 16.8522–19.3738× | 18.6061× |
| 1,000 | 4 | 19.6542 ± 0.2816× | 17.9931–19.9678× | 19.2408× |
| 1,000 | 6 | 21.4516 ± 0.2528× | 20.6036–21.7117× | 21.3261× |
| 10,000 | 1 | 7.3383 ± 0.0660× | 7.2224–7.6227× | 7.3383× |
| 10,000 | 2 | 6.4079 ± 0.1240× | 6.2637–6.5445× | 6.3918× |
| 10,000 | 4 | 4.8534 ± 0.1977× | 4.5480–5.7183× | 4.9441× |
| 10,000 | 6 | 4.5148 ± 0.1738× | 4.2959–4.6923× | 4.4899× |

Pyramidalization scientific agreement is exact at f32 output rounding on the
matched represented inputs. Its small native runtimes make process/serialization
costs particularly visible. Do not generalize these ratios to volume, Sterimol,
different molecule sizes, hardware, densities, radii or convergence settings.

Only 11 ligand chemistries/56 conformers are represented; the larger workloads
repeat them with fresh parsing/calculation per file. No quantum geometry generation,
ensemble thermodynamic reduction, reaction model or database workload is timed.
The volume grid and region normalization retain the audit's finite-lattice limits.
The exact gate is scoped to these geometries/settings and does not make morfeus
physical ground truth. The count check is regional equality, not a claim that
every occupied voxel's identity is identical at different floating-point precision.

The old approximately 14× morfeus study remains [historical](../../study_008/STUDY_008.md).
The separate 5.46× StericX optimization result concerns different old/new native
batch behavior. Neither historical ratio is combined with, or used to derive,
any result here. The current report supersedes the README's historical morfeus headline.

## 14. Reproduction and evidence inventory

Run from the repository root, with the production sources matching the recorded
commit/hashes, Rust/Cargo, a C compiler, uv, Python 3.12 and Linux `taskset` available:

```sh
python3 docs/benchmarks/morfeus_current/verify.py
python3 docs/benchmarks/morfeus_current/reproduce.py \
  --output docs/benchmarks/morfeus_reproduction --python python3.12
```

Use a **new output directory**. The reproduction command uses the preserved
56 XYZ files, source hashes, dependency lock and original repeated-file manifest;
it does not require the original hidden `.stericx` workspace. It rebuilds the
unmodified code and thin drivers, rechecks all scientific gates, and only then
performs alternating measurements. It refuses changed production sources or
failed gates. On another CPU topology, review/update the declared affinity masks
before creating a new experiment and state that difference. New timings and
binary hashes can differ with hardware/toolchains; do not overwrite this evidence.

A [fresh-directory reproduction smoke test](reproduction_smoke.json) rebuilt the
current source and drivers, installed the pinned environment, and passed the
complete 56-case scientific comparison, exact grid gate and all 1/2/4/6-worker
volume preflights. `--prepare-only` stops there without repeating the timing
campaign. Its verified archive is retained in `raw/reproduction-smoke.tar.gz`.

Individual stages are `compare.py`, `exact_volume_gate.py`, `benchmark.py prepare`,
`benchmark.py run`, `benchmark.py summarize`, `diagnose_sterimol.py`,
`pyramid_benchmark.py`, and `memory.py`. The original acquisition/build routine is
`freeze.py`; rerunning it queries the then-current PyPI version, whereas
`reproduce.py` pins this experiment's environment.

Machine-readable [benchmark_results.json](benchmark_results.json),
[raw_timings.json](raw_timings.json), per-run stdout/stderr/resource receipts,
and the independent [pyramidalization samples](pyramidalization/raw_timings.json)
retain every launch. Compressed stdout decompresses to its recorded SHA-256.
The [evidence inventory](evidence_manifest.json) hashes all published files;
disposable compiler targets, virtual environment and regenerated repeated copies
are omitted, with their reproduction inputs retained. Failed development adapter
attempts are kept in `raw/` and are not timing samples.
