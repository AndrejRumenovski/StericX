"""Publish only the campaign launched after the frozen scientific policy passed."""
import json,hashlib,tarfile
from pathlib import Path
R=Path(__file__).resolve().parent;H=R.parent;ROOT=H.parents[2];C=R/'admitted_campaign'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def num(x):return 'undefined' if x is None else f'{x:.8g}'

def main():
 assert read(C/'completion.json')['passed']
 frozen=read(R/'policy_freeze.json');assert sha(H/'EQUIVALENCE_POLICY.md')==frozen['sha256']
 gate=read(R/'evaluation_2/scientific_gate.json');assert gate['passed'] and not gate['failures']
 results=read(C/'benchmark_results.json');assert results['measurement_gate']
 samples=read(C/'raw_timings.json');assert len(samples['samples'])==120 and samples['status']=='complete'
 assert frozen['frozen_utc']<read(R/'evaluation_1/execution_freeze.json')['started_utc']<gate['completed_utc']<samples['started_utc']
 summary=results['summaries'];by={(r['structures'],r['workers']):r for r in summary};s1=by[10000,1];s6=by[10000,6]
 one=s1['paired_speedup']['median'];six=s6['paired_speedup']['median']
 with tarfile.open(R/'prior_publication.tar.gz') as t:old=t.extractfile('BENCHMARK.md').read().decode()
 # Reuse only environment/input/convention prose; remove the old P admission paragraph.
 static=old[old.index('## 2. Hardware'):old.index('## 6. Scientific agreement')]
 static=static[:static.index('Pyramidalization uses')]
 intro=f'''# Current StericX versus morfeus: prospective finite-grid volume benchmark

## 1. Executive result

On **10,000 XYZ files cyclically repeating 56 conformers from 11 ligands**, StericX
was **{one:.2f}× faster on one CPU core** and **{six:.2f}× faster at equal six-core
parallelism** than morfeus 0.8.0 for the specified **nine-output, three-frame,
finite-lattice buried-volume workload**. These are medians of four alternating
paired wall-time ratios, measured after the new scientific gate passed.
End-to-end API drivers include startup, parsing, calculation and complete JSON
output. This is not a timing of the full shipping descriptor CLI.

The policy was frozen at `{frozen['frozen_utc']}` with SHA-256
`{frozen['sha256']}` before its first application. All **2,588,544** point occupancy
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

'''
 science='''## 6. Scientific equivalence and excluded descriptors

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
'''
 metrics=read(R/'evaluation_2/scientific_summary.json')['metrics']['source_float64']
 metrics.update(read(R/'admitted_scalar_comparisons.json')['metrics'])
 for k,v in metrics.items():science+=f"| {k} | {v['N']} | "+' | '.join(num(v[key]) for key in ['MAE','RMSE','median_absolute_error','maximum_absolute_error','R_squared','slope','intercept'])+f" | {'finite-volume gate passed' if k in gate['admitted_descriptors'] else 'unresolved/excluded'} |\n"
 science+='''
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
'''
 for row in summary:
  for tool,v in row['tools'].items():
   wall=v['wall_ns'];science+=f"| {row['structures']:,} | {row['workers']} | {tool} | {wall['median']/1e9:.4f} ± {wall['MAD']/1e9:.4f} | {wall['minimum']/1e9:.4f}–{wall['maximum']/1e9:.4f} | {v['molecules_per_second']:,.1f} | {v['cpu_total_s']['median']:.4f} | {v['peak_rss_bytes']['median']/2**20:.2f} |\n"
  assert row['tools']['stericx']['wall_ns']['maximum']<row['tools']['morfeus']['wall_ns']['minimum']
 science+='''
*For multiprocessing, wait4 reports a maximum individual-process high-water mark,
not total process-tree memory. CPU/user/system/RSS median, MAD, full range and
individual samples are in [benchmark_results.json](benchmark_results.json).

## 9. Single-core implementation/runtime efficiency

'''
 for workers,title in [(1,None),(6,'## 10. Equal-core parallel throughput')]:
  row=by[10000,workers];ratio=row['paired_speedup']
  if title:science+='\n'+title+'\n\n'
  science+=f"At {workers} {'core' if workers==1 else 'physical cores'}, StericX's median wall time is **{row['tools']['stericx']['wall_ns']['median']/1e9:.4f} s**, versus **{row['tools']['morfeus']['wall_ns']['median']/1e9:.4f} s** for morfeus. Paired speedup is **{ratio['median']:.4f}×**, MAD {ratio['MAD']:.4f}, full range {ratio['minimum']:.4f}–{ratio['maximum']:.4f}×.\n\n"
  science+=('This is the primary per-core comparison; it includes runtime and driver costs and does not isolate programming language alone.\n' if workers==1 else 'This compares six native threads with six Python workers on the same six physical cores, including multiprocessing overhead. It measures parallel throughput, not a per-core or Rust-only improvement.\n')
 science+='''
| Files | Equal workers/cores | Paired speedup median ± MAD | Full paired range | Ratio of wall medians |
|---:|---:|---:|---:|---:|
'''
 for row in summary:
  v=row['paired_speedup'];science+=f"| {row['structures']:,} | {row['workers']} | {v['median']:.4f} ± {v['MAD']:.4f}× | {v['minimum']:.4f}–{v['maximum']:.4f}× | {row['ratio_of_wall_medians']:.4f}× |\n"
 science+='''
## 11. Throughput

Throughput is file count divided by median wall time. These are files representing
repeated conformers, not independent chemistries.
The final column uses single-worker morfeus as an ordinary-use baseline. For
multi-core rows it is explicitly an unequal-core practical throughput ratio.

| Tool | Workers | Median wall, s | Molecules/s | Relative to ordinary morfeus 1 worker |
|---|---:|---:|---:|---:|
'''
 for w in [1,2,4,6]:
  for tool,v in by[10000,w]['tools'].items():science+=f"| {tool} | {w} | {v['wall_ns']['median']/1e9:.4f} | {v['molecules_per_second']:,.1f} | {s1['tools']['morfeus']['wall_ns']['median']/v['wall_ns']['median']:.3f}× |\n"
 science+='''
## 12. Memory

Separate 10,000-file runs sample the simultaneous sum of parent/worker/resource-
tracker RSS, targeting 10 ms intervals. Individual process peak values are never
added. These runs occur after speed measurements and do not supply headline times.

| Tool | Workers | Sampled peak simultaneous sum RSS, MiB | Largest sample spacing, ms |
|---|---:|---:|---:|
'''
 memory=read(C/'memory_results.json');assert memory['complete']
 for v in memory['results']:science+=f"| {v['tool']} | {v['workers']} | {v['sampled_peak_sum_rss_bytes']/2**20:.2f} | {v['max_sample_spacing_s']*1000:.2f} |\n"
 science+='''
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
'''
 (H/'BENCHMARK.md').write_text(intro+static+science)
 results.update(publication_status='admitted_after_frozen_policy',scope='Nine finite-lattice buried-volume outputs, three frames, 56 fixed conformers',raw_timings='revision_policy/admitted_campaign/raw_timings.json',policy_sha256=frozen['sha256'],scientific_evaluation='revision_policy/evaluation_2/scientific_gate.json',acceptance_gates={'scientific':True,'fairness':True,'measurement':True,'documentation':True},memory_results='revision_policy/admitted_campaign/memory_results.json',prior_timings_used=False)
 write(H/'benchmark_results.json',results)
 write(H/'PUBLICATION_STATUS.json',dict(status='admitted_finite_grid_volume_after_frozen_policy',prospective_scientific_admission=True,policy_sha256=frozen['sha256'],scope=results['scope'],excluded=gate['excluded'],prior_campaign_status='pre_policy_exploratory'))
 readme=ROOT/'README.md';text=readme.read_text();a=text.index('## Performance vs morfeus');z=text.index('## Quick Start',a)
 section=f'''## Performance vs morfeus

On a Ryzen 5 5600G (six physical cores, Linux), StericX and morfeus 0.8.0
processed the same **10,000 XYZ files repeating 56 conformers from 11 ligands**.
Both calculated nine buried-volume outputs over three donor-plane grids, using
the same finite-grid estimator, radii, center and hydrogen conventions.

| Buried-volume API-driver benchmark | morfeus median time (files/s) | StericX median time (files/s) | Paired speedup |
|---|---:|---:|---:|
'''
 for w in [1,6]:
  row=by[10000,w];m=row['tools']['morfeus'];s=row['tools']['stericx'];section+=f"| {w} {'core: 1 process / 1 thread' if w==1 else 'cores: 6 processes / 6 threads'} | {m['wall_ns']['median']/1e9:.3f} s ({m['molecules_per_second']:,.0f}) | {s['wall_ns']['median']/1e9:.3f} s ({s['molecules_per_second']:,.0f}) | **{row['paired_speedup']['median']:.2f}×** |\n"
 section+='''
The first row measures single-core implementation/runtime efficiency; the second
measures equal-core batch throughput. Both include startup, parsing and JSON output.
Before these fresh timings, the [frozen scientific policy](docs/benchmarks/morfeus_current/EQUIVALENCE_POLICY.md)
passed all 56 geometries: pointwise occupancy matched a certified interval reference,
and output differences were exactly accounted for by independent rounding models.
This establishes the fixed finite-grid calculation, **not continuum-volume accuracy**.
Sterimol, pyramidalization and combined-descriptor speedups remain excluded.

These workload-specific timings use thin drivers of unchanged APIs, not the full
shipping descriptor CLI. [Methodology, exclusions, raw evidence and reproduction](docs/benchmarks/morfeus_current/BENCHMARK.md).
The initial pre-policy claims are withdrawn and archived; the older approximately
14× study remains historical. Neither is combined with the separate native batch
optimization result.

'''
 text=text[:a]+section+text[z:]
 text='\n'.join(f'| Current finite-grid buried-volume API drivers vs `morfeus`; 10,000 files repeating 56 conformers | **{one:.2f}× at 1 core** · **{six:.2f}× at equal 6 cores** · [policy and measurements](docs/benchmarks/morfeus_current/BENCHMARK.md) |' if line.startswith(('| Current StericX–morfeus comparison |','| Earlier exploratory buried-volume API-driver measurements vs `morfeus`;')) else line for line in text.split('\n'))
 readme.write_text(text)
 print(json.dumps({'single_core':one,'equal_six_core':six,'policy_sha256':frozen['sha256']}))

if __name__=='__main__':main()
