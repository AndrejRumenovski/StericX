raise SystemExit("Archived pre-policy campaign: publication claim withdrawn. See BENCHMARK.md and EQUIVALENCE_POLICY.md.")
"""Render measured results only; refuses incomplete or failed gates."""
import json
from pathlib import Path
from benchmark import HERE,read,write,sha


def fmt(x):return f"{x:.6g}"
def ms(x):return f"{x/1e9:.6f}"
def dist(d,scale=1):return f"{d['median']/scale:.4f} ± {d['MAD']/scale:.4f}"
def spread(d,scale=1):return f"{d['minimum']/scale:.4f}–{d['maximum']/scale:.4f}"
def select(results,n,w):return next(s for s in results["summaries"] if s["structures"]==n and s["workers"]==w)


def timing_table(results):
    lines=["| Files | Workers | Tool | Wall median ± MAD, s | Wall range, s | Files/s | CPU median, s | wait4 RSS median, MiB* |","|---:|---:|---|---:|---:|---:|---:|---:|"]
    for s in results["summaries"]:
        for tool,t in s["tools"].items():
            lines.append(f"| {s['structures']:,} | {s['workers']} | {tool} | {dist(t['wall_ns'],1e9)} | {spread(t['wall_ns'],1e9)} | {t['molecules_per_second']:,.1f} | {t['cpu_total_s']['median']:.4f} | {t['peak_rss_bytes']['median']/2**20:.2f} |")
    lines.append("\n*For multiprocessing, wait4 RSS is the largest individual process high-water mark, **not summed worker memory**. Full CPU/user/system/RSS MAD, range and samples are in JSON.")
    return "\n".join(lines)


def ratio_table(results):
    lines=["| Files | Equal workers/cores | Paired speedup median ± MAD | Full paired range | Ratio of wall medians |","|---:|---:|---:|---:|---:|"]
    for s in results["summaries"]:
        lines.append(f"| {s['structures']:,} | {s['workers']} | {dist(s['paired_speedup'])}× | {spread(s['paired_speedup'])}× | {s['ratio_of_wall_medians']:.4f}× |")
    return "\n".join(lines)


def main():
    primary=read(HERE/"benchmark_results.json")
    pyramid=read(HERE/"pyramidalization/benchmark_results.json")
    memory=read(HERE/"memory_results.json")
    assert primary["status"]==pyramid["status"]=="complete" and memory["complete"]
    assert primary["scientific_gate"]["passed"] and pyramid["scientific_gate"]["passed"]
    assert primary["measurement_gate"] and pyramid["measurement_gate"]
    assert read(HERE/"protocol_review.json")["passed"]
    assert read(HERE/"postflight_inputs.json")["passed"]
    one,six=select(primary,10000,1),select(primary,10000,6)
    s1,m1=one["tools"]["stericx"],one["tools"]["morfeus"]
    s6,m6=six["tools"]["stericx"],six["tools"]["morfeus"]
    scientific=read(HERE/"scientific_summary.json")
    freeze=read(HERE/"freeze.json")
    native=sha(HERE/"build/stericx")
    wrapper=sha(HERE/"build/volume-batch")
    numerical=["| Descriptor | N | MAE | RMSE | Max absolute | Median absolute | R² (identity) | Slope | Intercept |","|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for key,v in scientific["metrics"]["source_float64"].items():
        numerical.append("| "+" | ".join([key,str(v["N"]),*["undefined" if v[k] is None else fmt(v[k]) for k in ["MAE","RMSE","maximum_absolute_error","median_absolute_error","R_squared","slope","intercept"]]])+" |")
    memory_table=["| Tool | Threads/workers | Sampled peak simultaneous sum RSS, MiB | Largest sample spacing, ms |","|---|---:|---:|---:|"]
    for r in memory["results"]:memory_table.append(f"| {r['tool']} | {r['workers']} | {r['sampled_peak_sum_rss_bytes']/2**20:.2f} | {r['max_sample_spacing_s']*1000:.2f} |")
    throughput=["| Tool | Workers | Wall median, s | Files/s | Relative to ordinary morfeus 1 worker |","|---|---:|---:|---:|---:|"]
    for w in [1,2,4,6]:
        for tool,t in select(primary,10000,w)["tools"].items():throughput.append(f"| {tool} | {w} | {t['wall_ns']['median']/1e9:.4f} | {t['molecules_per_second']:,.1f} | {m1['wall_ns']['median']/t['wall_ns']['median']:.3f}× |")
    py1,py6=select(pyramid,10000,1),select(pyramid,10000,6)
    report=f'''# Current StericX versus morfeus: frozen end-to-end descriptor benchmarks

## 1. Executive result

On **10,000 files cyclically repeating 56 real conformers from 11 ligands**, the
current unmodified StericX buried-volume API was **{one['paired_speedup']['median']:.2f}× faster
on one CPU core** and **{six['paired_speedup']['median']:.2f}× faster at equal six-core parallelism**
than morfeus 0.8.0, using end-to-end command-line API drivers on this host.
These are medians of four fresh paired `morfeus wall / StericX wall` ratios.
No historical speed ratio enters the calculation.

The headline covers **nine buried-volume outputs, including all three donor-plane
calculations**, with exactly equal regional point populations and occupied counts.
Every floating-point residual is reconstructed exactly from those independent
counts using the native documented f32 rounding operations. No error tolerance
was invented or inferred from a historical maximum.

Separate pyramidalization timing gives **{py1['paired_speedup']['median']:.2f}× at one core**
and **{py6['paired_speedup']['median']:.2f}× at six cores** on the same 10,000-file repetition.
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

StericX commit: `{freeze['stericx_commit']}`; Rust 1.97.0
(`2d8144b78`, LLVM 22.1.6). Production source hashes, exact command and environment
are in [freeze.json](freeze.json). The pre-existing unrelated `docs/media/` files
were left untouched. No production source was changed or optimized for this study.

Build: `cargo build --locked --offline --release --bin stericx --target-dir …`.
Release uses opt-level 3, thin LTO, one codegen unit, no optional Cargo features,
no RUSTFLAGS, and rustc's default x86-64 target CPU. Enabled target features:
`fxsr`, `sse`, `sse2`; host capabilities are recorded separately.

- Shipping executable SHA-256: `{native}`.
- Timed volume driver SHA-256: `{wrapper}`.
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

{chr(10).join(numerical)}

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

{timing_table(primary)}

## 9. Single-core implementation/runtime comparison

For 10,000 files, StericX median wall time is **{s1['wall_ns']['median']/1e9:.4f} s**;
morfeus is **{m1['wall_ns']['median']/1e9:.4f} s**. Paired speedup is
**{one['paired_speedup']['median']:.4f}×**, MAD {one['paired_speedup']['MAD']:.4f},
range {spread(one['paired_speedup'])}×. Both are pinned to CPU 2, with one
native thread or one Python process and numerical-library threads limited to one.
This includes driver and runtime costs, so it is not a pure-language or algorithm-only attribution.

## 10. Equal-core parallel throughput

For 10,000 files and six physical cores, native six-thread median wall is
**{s6['wall_ns']['median']/1e9:.4f} s**; morfeus six-worker median is
**{m6['wall_ns']['median']/1e9:.4f} s**. Paired speedup is
**{six['paired_speedup']['median']:.4f}×**, MAD {six['paired_speedup']['MAD']:.4f},
range {spread(six['paired_speedup'])}×. This is an equal-hardware batch comparison,
including Python multiprocessing costs; it is not attributed wholly to Rust.

{ratio_table(primary)}

## 11. Throughput

The following table is **10,000 repeated input files**, not 10,000 independent
chemistries. The last column uses ordinary single-worker morfeus as a practical
baseline and must not be interpreted as per-core efficiency for multi-core rows.

{chr(10).join(throughput)}

StericX's six-thread versus ordinary morfeus one-worker wall-median ratio is
{m1['wall_ns']['median']/s6['wall_ns']['median']:.3f}×. This is explicitly an unequal-core practical
throughput comparison and is excluded from the equal-core README headline.

## 12. Memory

Timed wait4 peaks are in the raw summary above. For a single-process tool they
provide the kernel high-water RSS; for multiprocessing they **do not sum the tree**.
A separate 10,000-file execution per configuration samples concurrent parent,
worker and resource-tracker RSS, targeted every 10 ms. Sampling overhead does
not affect the headline timings. Peaks below are the largest **simultaneous**
sum; individual process peaks are never added together.

{chr(10).join(memory_table)}

Summed RSS counts shared pages in each process, so it is not unique physical
memory or PSS. Sampling can miss brief peaks; the actual maximum spacing is
shown. These are separate sampled runs, not repeated RSS distributions. Full
process/PID/RSS time series, peak snapshots, output checks and commands are in
[memory/](memory/) and [memory_results.json](memory_results.json).

## 13. Separate pyramidalization results and limitations

{timing_table(pyramid)}

{ratio_table(pyramid)}

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
python3 docs/benchmarks/morfeus_current/reproduce.py \\
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
'''
    (HERE/"BENCHMARK.md").write_text(report)
    primary.update(pyramidalization="pyramidalization/benchmark_results.json",memory="memory_results.json",scientific_summary="scientific_summary.json",scope="End-to-end public-API drivers; headline is 9 buried-volume outputs over three donor-plane grids; excludes Sterimol/combined; no historical ratio arithmetic",acceptance_gates=dict(scientific=True,fairness=True,measurement=True,documentation=True))
    write(HERE/"benchmark_results.json",primary)
    # Concise README edits happen only after all measurements and exact gates pass.
    root=HERE.parents[2]
    path=root/"README.md";text=path.read_text()
    old="| Historical, pre-remediation single-core buried-volume benchmark against `morfeus`; not a current-build comparison | **13.8×** · **1,546 input structures** · **one machine** |"
    new=f"| Current buried-volume API drivers vs `morfeus`, 10,000 files repeating 56 conformers; direct equal-core measurements | **{one['paired_speedup']['median']:.2f}× at 1 core** · **{six['paired_speedup']['median']:.2f}× at 6 cores** · [methodology](docs/benchmarks/morfeus_current/BENCHMARK.md) |"
    assert old in text or new in text;text=text.replace(old,new)
    section=f'''## Performance vs morfeus

On a Ryzen 5 5600G (Linux, six physical cores), current StericX and morfeus 0.8.0
processed the **same 10,000 XYZ files: repetitions of 56 conformers from 11 ligands**.
Thin command-line drivers call the unchanged APIs for nine equivalent buried-volume
outputs over three donor-plane grids, including startup, parsing and JSON output.
All regional occupied-point counts agree exactly; output differences are fully
accounted for by f32 rounding. Sterimol and combined-descriptor speedups are excluded.

| Buried-volume benchmark | morfeus median time (files/s) | StericX median time (files/s) | Paired speedup |
|---|---:|---:|---:|
| 1 core: 1 process / 1 thread | {m1['wall_ns']['median']/1e9:.3f} s ({m1['molecules_per_second']:,.0f}) | {s1['wall_ns']['median']/1e9:.3f} s ({s1['molecules_per_second']:,.0f}) | **{one['paired_speedup']['median']:.2f}×** |
| 6 cores: 6 processes / 6 threads | {m6['wall_ns']['median']/1e9:.3f} s ({m6['molecules_per_second']:,.0f}) | {s6['wall_ns']['median']/1e9:.3f} s ({s6['molecules_per_second']:,.0f}) | **{six['paired_speedup']['median']:.2f}×** |

The first row measures single-core implementation/runtime efficiency; the second
measures equal-core batch throughput. These are workload-specific API-driver results,
not timings of the full `stericx descriptors` command or universal speedups.
[Full methodology, accuracy, raw samples, memory, separate pyramidalization results,
and reproduction commands](docs/benchmarks/morfeus_current/BENCHMARK.md).
The older approximately 14× comparison remains historical and is superseded here.

'''
    anchor="## Quick Start\n";assert anchor in text
    if "## Performance vs morfeus\n" not in text:text=text.replace(anchor,section+anchor,1)
    text=text.replace("| **008** | Head-to-head speed benchmark vs `morfeus` | About 14× faster", "| **008** | Historical pre-remediation speed benchmark vs `morfeus` | About 14× faster")
    text=text.replace("**Current corrected build:** the 10,000-file descriptor batch", "**StericX-only optimization comparison (not morfeus):** the 10,000-file descriptor batch")
    path.write_text(text)


if __name__=="__main__":main()
