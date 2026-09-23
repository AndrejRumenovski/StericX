# Current StericX versus morfeus: prospective finite-grid volume benchmark

## 1. Executive result

On **10,000 XYZ files cyclically repeating 56 conformers from 11 ligands**, StericX
was **22.77× faster on one CPU core** and **19.91× faster at equal six-core
parallelism** than morfeus 0.8.0 for the specified **nine-output, three-frame,
finite-lattice buried-volume workload**. These are medians of four alternating
paired wall-time ratios, measured after the new scientific gate passed.
End-to-end API drivers include startup, parsing, calculation and complete JSON
output. This is not a timing of the full shipping descriptor CLI.

The policy was frozen at `2026-09-22T21:17:02.244163+00:00` with SHA-256
`b0243f8c566b973f01eb09b7cac7b26186b990c92afb895031c77bc8e34d1a61` before its first application. All **2,588,544** point occupancy
decisions were certified against an independent exact-decimal interval reference;
both tools matched every decision (the native direct predicate probe is connected
to its optimized row cache by the [source-level argument](revision_policy/NATIVE_PREDICATE_PROOF.md)). Every plane and aggregate volume output matched
its own independently evaluated rational IEEE rounding model exactly. Cross-tool
float bit equality is not expected and is not claimed.

**Scope exclusions:** Sterimol L/B1/B5, pyramidalization P/alpha, combined descriptors,
and accuracy of the continuous union-of-spheres integral remain unresolved or
outside this policy. No speedup for these quantities is admitted.

**Correction of prior publication:** the first campaign was run without a dedicated
policy frozen before classification. Its claims were withdrawn; its complete
[publication and evidence](revision_policy/prior_publication.tar.gz) are preserved
as exploratory history. Prior residuals had already been seen, so this is not a
blinded study. The new policy follows exact predicates and independently derived
rounding, not those residuals. **Every timing reported below is a new launch after
scientific admission**; no old timing or historical ratio was reused.

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

## 6. Scientific equivalence and excluded descriptors

The descriptor-specific [frozen policy](EQUIVALENCE_POLICY.md) is authoritative.
Its [freeze receipt](revision_policy/policy_freeze.json), [reference implementation](revision_policy/reference.py),
[execution freeze](revision_policy/evaluation_2/execution_freeze.json), and
[passed scientific gate](revision_policy/evaluation_2/scientific_gate.json) establish
scope and ordering. The audit's observed maxima were never acceptance thresholds.

The reference reads XYZ decimals and constants as exact rationals, constructs
outward binary64 intervals for geometry, and uses integer arithmetic for ideal
lattice membership. Every occupied/free union predicate is certified without an
uncertainty allowance. Both actual grids map bijectively to the same ideal
indices, region signs and occupancy mask on all 168 frames. There are no ambiguous
predicates and no differing points. Exact integer counts then feed an independent
rational model with round-to-nearest, ties-to-even at each declared operation.
Pi comes from rational Machin-series bounds, not from either candidate's output.
The current native optimized region outputs and public API aggregates also match
this model exactly. This proves equivalence of the specified **finite estimator**,
not convergence to the continuous volume or correctness on untested geometries.

The first reference evaluation failed on 51 morfeus near/far rounding checks
(45 plane fields and 6 aggregate fields): the reference mistakenly modeled
Python built-in `sum` as sequential additions. Python 3.12.13 uses compensated
summation. The [source-backed correction](revision_policy/reference_revision.json)
models that operation sequence exactly; it changes no acceptance rule or threshold.
The [failed evaluation](revision_policy/evaluation_1/scientific_gate.json) and its
original reference code remain preserved. The second complete evaluation passes.
See the frozen [CPython implementation](revision_policy/cpython-3.12.13-bltinmodule.c)
and [upstream source](https://raw.githubusercontent.com/python/cpython/v3.12.13/Python/bltinmodule.c).
Reference-only arithmetic tests use 2,000 rational cases, 100-digit square-root
checks, exact halfway rounding cases and 600 synthetic compensated-sum cases.

All 14 descriptor residuals were recomputed untimed on all 56 identical geometries.
[Full molecule-level residuals](revision_policy/evaluation_2/scientific_comparisons.csv)
include absolute and relative differences for source-f64, controlled-f32 and
matched-angular-count diagnostic lanes. Relative difference is undefined when
the reference is zero. [Summary statistics](revision_policy/evaluation_2/scientific_summary.json)
retain every lane. The table uses the actual admitted nine-output expressions
from [admitted_scalar_comparisons.json](revision_policy/admitted_scalar_comparisons.json)
for volume, and the source-f64 diagnostic lane for excluded descriptors. The untimed
diagnostic percentage averages per-plane percentages; the timed expression takes
the percentage after volume averaging, with that distinct rounding sequence checked
explicitly. R² is agreement against
the identity line, not merely squared correlation; constant-series regression
and R² are undefined. These statistics describe agreement and do not set policy.

| Descriptor | N | MAE | RMSE | Median AE | Maximum AE | R² | Slope | Intercept | Admission |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| sterimol_l | 56 | 5.5695878e-07 | 6.6033569e-07 | 5.4328274e-07 | 1.6953073e-06 | 1 | 0.99999984 | 1.2001993e-06 | unresolved/excluded |
| sterimol_b1 | 56 | 0.0069007277 | 0.0089704956 | 0.0055310826 | 0.018063878 | 0.99976225 | 1.0028726 | -0.0041362133 | unresolved/excluded |
| sterimol_b5 | 56 | 9.5155559e-07 | 1.243546e-06 | 6.7966843e-07 | 2.7777423e-06 | 1 | 0.99999988 | 1.9410863e-06 | unresolved/excluded |
| pyr_p | 56 | 2.1431247e-08 | 2.7556289e-08 | 1.8566911e-08 | 7.5721954e-08 | 1 | 0.99999974 | 2.3883739e-07 | unresolved/excluded |
| pyr_alpha | 56 | 2.3672661e-06 | 2.9775105e-06 | 2.0570151e-06 | 6.6502892e-06 | 1 | 0.99999997 | 3.0695183e-07 | unresolved/excluded |
| buried_volume | 56 | 1.9524608e-06 | 2.5849422e-06 | 1.571806e-06 | 6.1395112e-06 | 1 | 0.99999994 | 2.6439083e-06 | finite-volume gate passed |
| percent_buried_volume | 56 | 1.387362e-06 | 1.7746078e-06 | 9.7645158e-07 | 4.264961e-06 | 1 | 1 | -2.3678626e-07 | finite-volume gate passed |
| qvbur_min | 56 | 3.2365998e-07 | 4.4132762e-07 | 2.1781569e-07 | 1.2837633e-06 | 1 | 0.99999988 | 9.962053e-07 | finite-volume gate passed |
| qvbur_max | 56 | 7.4613483e-07 | 9.4564858e-07 | 6.323603e-07 | 2.4031027e-06 | 1 | 0.99999996 | 2.3426627e-07 | finite-volume gate passed |
| max_delta_qvbur | 56 | 7.211646e-07 | 9.5218066e-07 | 5.3984144e-07 | 2.257357e-06 | 1 | 0.99999996 | 1.4628271e-07 | finite-volume gate passed |
| ovbur_min | 56 | 0 | 0 | 0 | 0 | undefined | undefined | undefined | finite-volume gate passed |
| ovbur_max | 56 | 7.1184968e-07 | 8.2034861e-07 | 7.7757653e-07 | 1.7964076e-06 | 1 | 0.99999998 | -3.0342824e-07 | finite-volume gate passed |
| near_vbur | 56 | 2.0632996e-06 | 2.7337191e-06 | 1.3993392e-06 | 8.2992415e-06 | 1 | 0.99999992 | 2.7363709e-06 | finite-volume gate passed |
| far_vbur | 56 | 2.1075604e-07 | 4.3265449e-07 | 3.1139593e-08 | 1.5759838e-06 | 1 | 0.99999998 | -2.1066461e-08 | finite-volume gate passed |

Units: L/B1/B5 Å; alpha degrees; P dimensionless; %Vbur percentage points;
all other volume outputs Å³. Each molecule's diagnostic interpretation is retained
in [descriptor_investigations.json](revision_policy/descriptor_investigations.json).
B1's distinct angular phase/resolution and B5's analytic-versus-sampled convention
are evaluated against the preserved independent analytic support diagnostic.
Neither that diagnostic nor f32-rounded P/alpha agreement supplies a certified
input/axis/trigonometric error budget. Those descriptors remain excluded.

## 7. Benchmark methodology

Only the admitted volume scope is timed. One warmup per tool/configuration is
followed by four alternating AB/BA paired runs. All 56, 1,000 and 10,000-file
workloads use 1, 2, 4 and 6 threads/workers, yielding 24 warmup launches and 96
measured launches; all samples are retained. No trimming or outlier removal.
Affinities and thread-library limits are recorded above and with each sample.
One worker is an ordinary morfeus Python loop; parallel morfeus uses spawn and
Pool.map. Both perform all three grids and return all nine scalar outputs.

The C launcher starts the monotonic wall timer before fork/exec and stops after
wait4. Startup/imports, worker creation/shutdown, parsing, calculation and ordinary
complete JSON output are included. Stdout/stderr files are opened by the harness
outside the timer; output compression and verification occur after timing. User
and system CPU time and peak RSS come from wait4. The corpus is warm-cache after
warmup; every repetition is parsed and calculated again, without descriptor cache.
All timed results must exactly reproduce each tool's admitted 56-row output at
every worker count. Scientific reference computations, builds and RSS sampling
were completed before or performed after the timed campaign. Routine desktop
activity and brief progress/file checks were not isolated from the machine.
At the user’s request, the scheduler was paused between six-core launches to
publish the earlier campaign. Its active measurement child finished normally;
git compression/upload occurred only after that child exited. No timed child was
paused and no samples were dropped. The gap is retained in timestamps and
[publication_pause.json](revision_policy/admitted_campaign/publication_pause.json).

Measurement gate: every paired ratio must exceed one, and wall-time ranges must
be disjoint for each same-size/same-worker comparison. These establish a measured
difference on this host, not a general statistical confidence claim about all
molecules or machines. Median/MAD/range and all individual samples are retained.

## 8. Timing, CPU and RSS summary

| Files | Workers | Tool | Wall median ± MAD, s | Wall range, s | Files/s | CPU median, s | wait4 RSS median, MiB* |
|---:|---:|---|---:|---:|---:|---:|---:|
| 56 | 1 | stericx | 0.0314 ± 0.0000 | 0.0313–0.0314 | 1,785.2 | 0.0307 | 2.95 |
| 56 | 1 | morfeus | 1.0615 ± 0.0035 | 1.0349–1.0660 | 52.8 | 1.0610 | 71.58 |
| 56 | 2 | stericx | 0.0182 ± 0.0000 | 0.0182–0.0185 | 3,079.1 | 0.0314 | 3.20 |
| 56 | 2 | morfeus | 0.9376 ± 0.0063 | 0.9151–0.9458 | 59.7 | 1.5808 | 71.64 |
| 56 | 4 | stericx | 0.0122 ± 0.0003 | 0.0117–0.0127 | 4,606.7 | 0.0339 | 3.98 |
| 56 | 4 | morfeus | 0.7638 ± 0.0064 | 0.7544–0.7752 | 73.3 | 2.1226 | 71.92 |
| 56 | 6 | stericx | 0.0109 ± 0.0004 | 0.0102–0.0116 | 5,119.3 | 0.0362 | 4.67 |
| 56 | 6 | morfeus | 0.7033 ± 0.0070 | 0.6952–0.7168 | 79.6 | 2.7559 | 71.83 |
| 1,000 | 1 | stericx | 0.5480 ± 0.0102 | 0.5319–0.5674 | 1,824.8 | 0.5368 | 4.43 |
| 1,000 | 1 | morfeus | 14.9944 ± 0.0905 | 14.8584–15.1605 | 66.7 | 14.9870 | 74.90 |
| 1,000 | 2 | stericx | 0.3148 ± 0.0033 | 0.3094–0.3218 | 3,176.1 | 0.5503 | 4.64 |
| 1,000 | 2 | morfeus | 8.2128 ± 0.0571 | 8.1017–8.3683 | 121.8 | 15.8760 | 72.50 |
| 1,000 | 4 | stericx | 0.1899 ± 0.0034 | 0.1841–0.1992 | 5,266.3 | 0.5635 | 4.76 |
| 1,000 | 4 | morfeus | 4.6810 ± 0.0408 | 4.6251–4.7690 | 213.6 | 17.3126 | 72.28 |
| 1,000 | 6 | stericx | 0.1554 ± 0.0043 | 0.1507–0.1717 | 6,434.1 | 0.5842 | 5.04 |
| 1,000 | 6 | morfeus | 3.4920 ± 0.0177 | 3.4681–3.5350 | 286.4 | 18.5412 | 72.42 |
| 10,000 | 1 | stericx | 5.6073 ± 0.1306 | 5.3962–5.7655 | 1,783.4 | 5.5042 | 17.38 |
| 10,000 | 1 | morfeus | 128.5149 ± 0.7860 | 127.0535–129.7907 | 77.8 | 128.4484 | 102.72 |
| 10,000 | 2 | stericx | 3.1288 ± 0.0583 | 3.0650–3.2515 | 3,196.1 | 5.5955 | 17.77 |
| 10,000 | 2 | morfeus | 75.7595 ± 1.7821 | 73.4339–78.3486 | 132.0 | 149.0517 | 104.78 |
| 10,000 | 4 | stericx | 1.8758 ± 0.0091 | 1.8318–1.8900 | 5,331.1 | 5.7171 | 18.56 |
| 10,000 | 4 | morfeus | 42.2995 ± 0.3988 | 41.8272–42.8797 | 236.4 | 161.1513 | 104.43 |
| 10,000 | 6 | stericx | 1.4764 ± 0.0245 | 1.4191–1.5165 | 6,773.2 | 5.9622 | 19.43 |
| 10,000 | 6 | morfeus | 29.1137 ± 0.0828 | 28.9655–30.2887 | 343.5 | 166.6694 | 104.37 |

*For multiprocessing, wait4 reports a maximum individual-process high-water mark,
not total process-tree memory. CPU/user/system/RSS median, MAD, full range and
individual samples are in [benchmark_results.json](benchmark_results.json).

## 9. Single-core implementation/runtime efficiency

At 1 core, StericX's median wall time is **5.6073 s**, versus **128.5149 s** for morfeus. Paired speedup is **22.7683×**, MAD 0.3703, full range 22.3418–24.0525×.

This is the primary per-core comparison; it includes runtime and driver costs and does not isolate programming language alone.

## 10. Equal-core parallel throughput

At 6 physical cores, StericX's median wall time is **1.4764 s**, versus **29.1137 s** for morfeus. Paired speedup is **19.9110×**, MAD 0.1910, full range 19.5902–20.4113×.

This compares six native threads with six Python workers on the same six physical cores, including multiprocessing overhead. It measures parallel throughput, not a per-core or Rust-only improvement.

| Files | Equal workers/cores | Paired speedup median ± MAD | Full paired range | Ratio of wall medians |
|---:|---:|---:|---:|---:|
| 56 | 1 | 33.8385 ± 0.0698× | 33.0769–33.9166× | 33.8385× |
| 56 | 2 | 51.1899 ± 0.3939× | 50.2968–51.8184× | 51.5530× |
| 56 | 4 | 62.9305 ± 1.9323× | 59.7668–65.3143× | 62.8315× |
| 56 | 6 | 65.2088 ± 2.4030× | 60.2780–68.4277× | 64.2911× |
| 1,000 | 1 | 27.1983 ± 0.3878× | 26.7196–28.2773× | 27.3610× |
| 1,000 | 2 | 26.1796 ± 0.1528× | 25.5274–26.4792× | 26.0849× |
| 1,000 | 4 | 24.4375 ± 0.5129× | 23.6322–25.9023× | 24.6512× |
| 1,000 | 6 | 22.5374 ± 0.5234× | 20.2731–23.2485× | 22.4679× |
| 10,000 | 1 | 22.7683 ± 0.3703× | 22.3418–24.0525× | 22.9193× |
| 10,000 | 2 | 23.8197 ± 0.2683× | 23.4221–25.4717× | 24.2138× |
| 10,000 | 4 | 22.7550 ± 0.1155× | 22.2083–22.9072× | 22.5501× |
| 10,000 | 6 | 19.9110 ± 0.1910× | 19.5902–20.4113× | 19.7192× |

## 11. Throughput

Throughput is file count divided by median wall time. These are files representing
repeated conformers, not independent chemistries.
The final column uses single-worker morfeus as an ordinary-use baseline. For
multi-core rows it is explicitly an unequal-core practical throughput ratio.

| Tool | Workers | Median wall, s | Molecules/s | Relative to ordinary morfeus 1 worker |
|---|---:|---:|---:|---:|
| stericx | 1 | 5.6073 | 1,783.4 | 22.919× |
| morfeus | 1 | 128.5149 | 77.8 | 1.000× |
| stericx | 2 | 3.1288 | 3,196.1 | 41.075× |
| morfeus | 2 | 75.7595 | 132.0 | 1.696× |
| stericx | 4 | 1.8758 | 5,331.1 | 68.512× |
| morfeus | 4 | 42.2995 | 236.4 | 3.038× |
| stericx | 6 | 1.4764 | 6,773.2 | 87.045× |
| morfeus | 6 | 29.1137 | 343.5 | 4.414× |

## 12. Memory

Separate 10,000-file runs sample the simultaneous sum of parent/worker/resource-
tracker RSS, targeting 10 ms intervals. Individual process peak values are never
added. These runs occur after speed measurements and do not supply headline times.

| Tool | Workers | Sampled peak simultaneous sum RSS, MiB | Largest sample spacing, ms |
|---|---:|---:|---:|
| stericx | 1 | 17.57 | 10.94 |
| morfeus | 1 | 102.68 | 12.76 |
| stericx | 2 | 18.02 | 10.44 |
| morfeus | 2 | 267.23 | 11.23 |
| stericx | 4 | 18.88 | 10.54 |
| morfeus | 4 | 415.01 | 15.06 |
| stericx | 6 | 19.80 | 10.71 |
| morfeus | 6 | 561.50 | 13.63 |

Summed RSS double-counts shared pages and is not PSS or unique physical memory.
Sampling can miss brief peaks. These are one separately sampled run/configuration,
not repeated RSS distributions. [Raw memory traces and outputs](revision_policy/admitted_campaign/memory/)
and [memory_results.json](revision_policy/admitted_campaign/memory_results.json)
retain snapshots, PIDs, actual spacing and checked outputs.

## 13. Limitations and history

Only 56 conformers/11 ligands and the frozen finite estimator are admitted.
No continuum convergence bound was established; discretization error relative to
continuous occupied volume is outside this performance claim. No claim covers
Sterimol, pyramidalization, the combined CLI, quantum geometry generation,
ensemble reduction, models/databases, different radii/densities or other machines.
The larger corpora repeat the same geometries; this is a throughput test, not
10,000 independent chemistries. The small corpus emphasizes startup overhead.

The policy was frozen before its new application, but after prior exploratory
exposure. Exact predicates and rational operation models avoid selecting a
residual tolerance. This is an internally constructed independent reference,
not an assertion of external approval or a new external scientific audit.

The first campaign's volume and pyramidalization timings and initial gate labels
are preserved as withdrawn exploratory history. The older approximately 14×
study remains [historical](../../study_008/STUDY_008.md). Neither it nor the separate
5.46× native optimization ratio contributes to this fresh result.

## 14. Reproduction and evidence

From the repository root, preserve the recorded production commit/source hashes.
The frozen Python environment, requirements, build flags, package source hashes
and input manifest are linked above. No production code was changed.

```bash
# Integrity only, no timing rerun:
python3 docs/benchmarks/morfeus_current/verify.py

# Fresh policy evaluation and timings in a new sibling directory:
python3 docs/benchmarks/morfeus_current/reproduce.py /absolute/path/to/repo/docs/benchmarks/new_morfeus_run
```

The reproduction command refuses an existing directory, copies only frozen
inputs/scripts and builds both unchanged API drivers and the untimed probe. It
recreates the pinned environment, runs reference-only arithmetic checks, freezes
and verifies the same policy, recomputes all residuals and predicate/rounding
checks, and times only after admission. See its `--prepare-only` option for the
untimed preparation and scientific check without the timing campaign.

- [Input manifest](input_manifest.json), [software freeze](freeze.json), [environment](environment/).
- [Policy](EQUIVALENCE_POLICY.md), [policy freeze](revision_policy/policy_freeze.json), [reference revision](revision_policy/reference_revision.json).
- [Fresh raw comparisons](revision_policy/evaluation_2/raw/), [molecule residuals](revision_policy/evaluation_2/scientific_comparisons.csv), [metrics](revision_policy/evaluation_2/scientific_summary.json), [investigations](revision_policy/descriptor_investigations.json).
- [Pointwise evidence and exact scalar models](revision_policy/evaluation_2/volume_records.json); `reference-00.npz` through `reference-55.npz` retain all interval clearances and masks.
- [New timed campaign](revision_policy/admitted_campaign/), [all raw timing samples](revision_policy/admitted_campaign/raw_timings.json), [machine-readable results](benchmark_results.json).
- [Pre-policy archive](revision_policy/prior_publication.tar.gz) and [revision provenance](revision_policy/revision.json).

The top-level evidence manifest seals this revision and records its predecessor.
Its verification establishes artifact integrity and checks admission/measurement
ordering; it is not a new scientific calculation or a timing rerun.
