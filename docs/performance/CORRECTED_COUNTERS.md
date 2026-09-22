# Corrected baseline modeled counter and heap investigation

Cachegrind and DHAT completed the unchanged `conformers_56`, `screen_1000` and `descriptors_10000` workloads on the frozen corrected baseline, each at one Rayon thread pinned to CPU2. All six outputs match the corresponding native fingerprints exactly. No workload was reduced or extrapolated into a claimed full-workload result. The other seven workload types were not profiled by these tools in this investigation.

The executable is the uninstrumented symbol build from corrected commit `1ecfbd5be8b354711bdadd8a8d148417e03d8e7a`, bound to the frozen baseline source and builds receipt. The profiling-feature allocator/timers are absent. Native wall/CPU/RSS results remain in [CORRECTED_PROFILE.md](CORRECTED_PROFILE.md); these profiler runs neither measure native speedups nor compare C1.

Fresh hardware-event probes were denied (`perf_event_paranoid=4`). No privileges were changed. The retained local Valgrind installation was freshly checked: valgrind-3.26.0, both tool/preload paths and actual executable probes succeeded. Tool, source, input, binary, output, native-oracle and helper identities are in [the complete counter data](current_corrected_counters.json).

## Modeled instructions, data references, branches and caches

Cachegrind explicitly models 32 KiB, 8-way, 64-byte-line I1 and D1 caches and a 16 MiB, 16-way, 64-byte-line unified last-level cache, using freshly read CPU2 sysfs geometry. Intermediate L2 is not modeled. Matching cache geometry does not reproduce the full processor or its branch predictor. Counts below are modeled events; REP iterations differ from native retired-instruction counting. Data references mean reads plus writes; branches mean conditional plus indirect branches.

| Workload | Instructions | Data references | Branches | Branch misses | D1 misses | LL misses |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `conformers_56` | 381,443,863 | 111,861,488 | 49,729,574 | 2,735,753 | 710,155 | 21,734 |
| `screen_1000` | 2,776,643,885 | 731,078,310 | 230,234,019 | 7,152,974 | 2,550,022 | 151,563 |
| `descriptors_10000` | 68,077,660,702 | 19,956,355,170 | 8,871,083,151 | 485,283,798 | 124,816,733 | 289,244 |

Instruction attribution uses the function totals from `cg_annotate` at a 0.1% display threshold, including each function's source-file/inline attributions exactly once. These are instruction shares, not elapsed-time or Amdahl shares; optimized inlining also makes their categories differ from lexical diagnostic scopes.

| Workload | Function | Modeled instruction share |
| --- | --- | ---: |
| `conformers_56` | `<steric_x::geometry::buried_volume::BuriedVolumeCalculator>::compute_from_center` | 82.099% |
| `conformers_56` | `steric_x::geometry::buried_volume::integration_grid` | 10.661% |
| `conformers_56` | `steric_x::geometry::sterimol::params_from_projection` | 3.162% |
| `conformers_56` | `<core::str::iter::Split<core::str::IsWhitespace> as core::iter::traits::iterator::Iterator>::try_fold::<(), core::iter::traits::iterator::Iterator::find::check<&str, &mut core::str::IsNotEmpty>::{closure#0}, core::ops::control_flow::ControlFlow<&str>>` | 0.770% |
| `screen_1000` | `steric_x::geometry::sterimol::params_from_projection` | 38.230% |
| `screen_1000` | `<core::str::iter::Split<core::str::IsWhitespace> as core::iter::traits::iterator::Iterator>::try_fold::<(), core::iter::traits::iterator::Iterator::find::check<&str, &mut core::str::IsNotEmpty>::{closure#0}, core::ops::control_flow::ControlFlow<&str>>` | 9.599% |
| `screen_1000` | `steric_x::model::domain::beta_continued_fraction` | 6.165% |
| `screen_1000` | `core::slice::sort::stable::quicksort::quicksort::<f64, <[f64]>::sort_by<<f64>::total_cmp>::{closure#0}>` | 5.294% |
| `descriptors_10000` | `<steric_x::geometry::buried_volume::BuriedVolumeCalculator>::compute_from_center` | 82.153% |
| `descriptors_10000` | `steric_x::geometry::buried_volume::integration_grid` | 10.667% |
| `descriptors_10000` | `steric_x::geometry::sterimol::params_from_projection` | 3.165% |
| `descriptors_10000` | `<core::str::iter::Split<core::str::IsWhitespace> as core::iter::traits::iterator::Iterator>::try_fold::<(), core::iter::traits::iterator::Iterator::find::check<&str, &mut core::str::IsNotEmpty>::{closure#0}, core::ops::control_flow::ControlFlow<&str>>` | 0.771% |

## Requested heap and allocation sites

DHAT intercepts whole-process heap operations, including libc and shutdown. A reallocation adds a block and its entire new requested size; its stack remains attributed to the original site. Global heap peak is measured jointly, not by summing site-local peaks. Requested traffic and peak heap are not native RSS, allocation time or evidence of a leak. Lifetimes use modeled instructions. The local version-matched Cachegrind and DHAT manuals are hash-bound in the receipt.

| Workload | Requested bytes | Allocation blocks incl. realloc | Peak live heap bytes | Grid-attributed traffic |
| --- | ---: | ---: | ---: | ---: |
| `conformers_56` | 23,351,167 | 18,404 | 432,501 | 94.300% |
| `screen_1000` | 91,065,573 | 934,690 | 5,303,806 | 0.000% |
| `descriptors_10000` | 4,115,598,060 | 3,161,379 | 14,250,478 | 95.543% |

The JSON retains the ten largest individual allocation sites and their complete stacks, plus Rust-counter cross-checks with explicit scope differences. Grid-attributed traffic is summed over every DHAT program point whose stack includes `integration_grid`, not inferred from timing shares or only the displayed top sites.

These baseline results support candidate prioritization. Exact scientific gates and paired uninstrumented native measurements remain required for an optimization claim.
