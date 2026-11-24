# Accuracy, work and performance diagnostics

All new wall times were collected on the ordinary shared workstation. They describe
costs of these inputs and do not establish controlled general speedups. Numerical
admission uses independent intervals and unchanged tolerances, not timing. The frozen
production baseline and explicit new methods remain separate binaries/APIs.

## Original numerical requests and expanded budgets

All eighteen original failed requests were run in three alternating baseline/new
process pairs. The new method meets **18/18** targets; the frozen baseline meets
**0/18**. Repeating both with four times each request's original cell budget still
leaves the baseline at **0/18**, while the new method remains **18/18**. No tolerance
was relaxed. All original and expanded independent analytical reference intervals
remain enclosed. The largest new error bound among these requests is
3.197442310920451e-14 Å³ (the subnormal requests retain their own much smaller target).

The improvements establish analytical geometric proofs, so these cases need zero
adaptive subdivisions. The baseline spends work on conservative generic boundaries
or reaches its floating-point limit. A wall-time ratio would compare different
achieved accuracies and is therefore not an equal-accuracy speedup claim. Original
baseline process user CPU is 0.14–0.15 s for the eighteen requests; the new process
rounds below /usr/bin/time's 0.01-s display resolution. First-pair peak RSS is
7,844 versus 2,568 KiB, including input parsing and process startup.

A representative three-sphere chain allocates five times, requesting 296 bytes in
the new kernel. The allocation snapshot returns to its prior live-byte count.
These counts exclude input parsing and serialization, unlike process RSS.

## Regional and profile work

For twelve frozen molecular regional requests at R=3.5 Å, tolerance 0.05 Å³,
budget 1,048,576, all certificates meet target. They require 30,864–138,634 cells
and 29,673–135,197 allocation calls per request. The largest process-wide tracked
live heap across that serial batch is 15,955,688 bytes; RSS is 19,280 KiB. The first
request allocates 4,358,064 total requested bytes in 29,673 calls and twelve
reallocations. Requested allocation traffic is not live memory or RSS.

One four-entry molecular radius profile (2.5,3.5,4.5,3.5 Å) requires three distinct
integrations and reuses the fourth certificate. It makes 76,062 allocation calls,
35 reallocations, and requests 13,516,560 bytes cumulatively. Tracked peak live bytes
are 4,013,462; returned result storage retains 5,696 bytes at the snapshot. The
instrumented result equals the original certificate after excluding timing/allocation
fields. Distinct radii still have separate cell trees; shared preprocessing was only
0.000846% of median profile time in the independent profile campaign. There is no
evidence here for a substantial distinct-radius speedup from preprocessing alone.

The 405 experimental molecular frames all meet the same 0.05 Å³ target. Their full
native workflow, including other descriptors/statistics and JSON output, uses
85.27 s user CPU, 85.50 s elapsed and 52,624 KiB peak RSS. This workflow time must not
be assigned independently to every region. The mathematical API reports shared
runtime and work counters with per-region bounds/status.

## Parallel outer-batch diagnostic

Parallelism is over independent regional requests using a research Rayon adapter;
the individual regional API remains serial. All bounds, estimates, counters, statuses
and output order match exactly across thread counts. Three samples per count:

| Threads | Median elapsed seconds | Observed ratio to one thread | Ratio / threads |
|---:|---:|---:|---:|
| 1 | 0.63186 | 1.000 | 1.000 |
| 2 | 0.34491 | 1.832 | 0.916 |
| 4 | 0.21186 | 2.982 | 0.746 |
| 6 | 0.16547 | 3.819 | 0.636 |

These are exploratory batch measurements under shared load, not controlled
parallel-efficiency guarantees. They include scheduling/allocator effects and a
small, uneven twelve-request workload. No intra-region parallel algorithm is claimed.

## Instructions and memory correctness

Hardware `perf` instructions/cycles are unavailable (`perf_event_paranoid=4`). The
failed command is retained and no system permissions were changed. Callgrind 3.26.0
provides **simulated whole-process instruction counts**, including startup, parsing
and serialization:

| Diagnostic process | Callgrind instructions |
|---|---:|
| Frozen chain, original 1,023-cell budget | 3,308,352 |
| New analytical chain proof | 420,061 |
| One molecular regional request | 194,104,638 |
| Three-distinct-radius profile plus duplicate | 507,281,471 |

These counts describe different achieved accuracies for the old/new chain and are
not hardware performance counters. Native numerical output matches ordinary runs
exactly after removing named time/allocation diagnostics. Memcheck for regional and
profile workflows reports **zero errors and zero definitely/indirectly/possibly lost
bytes**, without added suppressions. The regional process retains 544 reachable
runtime bytes at exit, which is not silently labeled a leak-free zero-byte process.

The first Valgrind run used a library path containing a space; LD_PRELOAD split it,
so interception failed and memcheck emitted errors. That probe's memory/instruction
evidence is rejected and preserved. A whitespace-free symlink to the identical
library fixes loading; both workflows then pass. An initial comparison also omitted
`preparation_seconds` from its timing-field list; its sole mismatch was that field,
and the corrected checker retains every scientific field and work count.

Evidence: `performance_v1/`, `valgrind_v2/`, `profile_allocations_v1/` and earlier
method-specific run directories under `.stericx/experimental_validation/`. Controllers
are `profile_methods.py`, `profile_valgrind.py` and `verify_profile_outputs.py`.
The evidence supports the accepted exact geometric shortcuts and shared regional
calculation. Further allocator or tree optimization is not needed to make the
chemical conclusions trustworthy; better structures and experimental coverage are
currently more scientifically consequential.
