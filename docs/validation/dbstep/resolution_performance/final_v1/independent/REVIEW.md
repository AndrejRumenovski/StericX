# Independent performance review

The controlled campaign supports performance preservation for the frozen workload
set. No native or Python configuration crosses the unchanged 10% investigation
trigger. This review does not independently admit a release: scientific, source
binding, packaging and release artifact gates remain separate requirements.

The audit checks all 504 native attempts in their exact expected order, including
56 warmups and 448 measured attempts (eight AB/BA pairs across 28 configurations).
It hashes each metrics file and raw stdout/stderr, recomputes scientific output
fingerprints with the established output parser, and independently recomputes
wall-time medians, minima, maxima, MAD and pair ratios. All samples and unusually
slow observations remain included. All 28 instrumented runs produce the same
scientific fingerprints as their corresponding native workload.

All four Python campaigns follow A/B/B/A order. Each contains eight blocks of
native/bridge/public methods, alternating method order and using 100 calls per
block: 96 block attempts total. All preflight results agree exactly, and all
32 timed native block outputs agree exactly. The bridge and public methods are
checked at preflight as specified by the unchanged protocol; this is not a claim
that every Python return value was independently rechecked during timing.

Six-thread 10,000-file throughput is 10,656.50 files/s baseline and 10,774.50
files/s candidate. Wall-time medians are 938.394 and 928.118 ms. The maximum
native side-median slowdown is 4.843% for a roughly 2 ms single-molecule process.
Python native/bridge/public median ratios are 0.997744/0.992293/0.995386. Identical
production executable and installed package bytes mean these timing differences
are host/sample variation, not evidence of algorithmic improvement.

All 14 native allocation pairs have identical allocation, reallocation and
deallocation counts, with no failed allocations. Whole-process requested bytes
differ by +2 to +12; baseline/candidate executable and output paths have different
lengths, so whole-process counters do not isolate the descriptor call. The six-
conformer warmed call does isolate descriptor/statistics work and has exactly
871 allocations, 73 reallocations, 871 deallocations and 2,454,897 requested
bytes on both sides, with zero live-byte delta. Process-wide peak requested
storage differs by one byte for that example. Six-thread native peak storage can
also vary with overlapping worker lifetimes.

The largest RSS median increase is 1.77% (database build, one thread). For the
six-thread 10,000-file workload, medians are 26.334 and 26.457 MiB. RSS is measured
per child process and is distinct from allocator requested storage.

The recorded interval is 261.128 seconds. Aggregate CPU counters show 82.41%
idle, 0.36% I/O wait and 17.23% other work. CPU pressure accumulates 0.214% of the
interval; I/O pressure accumulates 0.604% some/0.589% full; memory pressure does
not increase. All native attempts have zero major page faults. These counters
include benchmark work and are not a proof of total host isolation. The recorded
`ps %CPU` values are lifetime averages, not instantaneous interference measures.
No evidence here warrants excluding samples or repeating the campaign.

The first independent audit attempt incorrectly expected a `values` member in
native wall-time summaries. That reviewer exception is preserved in
`audit_attempt1.py` and `attempt_1_failure.txt`. The corrected audit reads the raw
samples and verifies the actual summary fields. The measurement campaign,
scientific outputs, summary values and acceptance criteria were never changed.

Reproduce this read-only audit from the repository root:

```sh
.venv/bin/python .stericx/dbstep_resolution/performance_independent_v1/audit.py
```

The machine-readable receipt is `review.json`; `attempt_2.log` records the
successful execution. This audit creates only its own derived report.
