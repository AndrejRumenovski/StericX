# Measurement-window review and prospective artifact handling

The interrupted window 2 and complete CLI attempts in windows 4 and 5 remain
retained. No failed host gate is silently promoted to a pass. Window 3's guard
refused to suspend an ancestor of the measurement process; subsequent attempts
use an independent transient user service and separate automatic-restore timer.

Window 4 completed all 630 CLI trials with exact scientific outputs and passed
its inner quiet endpoints. After its final journal write, the outer postflight
reported 0.168 I/O-wait core equivalents, above the unchanged 0.05 criterion.
The journal was 172,272,188 bytes, including 162,387,101 bytes of repeated command
arguments before indentation. Full arguments already existed in every trial's
raw metrics file. Compact references reduce the equivalent journal to 950,578
bytes; all 630 raw metric dictionaries were checked without changing a field.

Window 5 repeated all 630 exact-output trials using compact journals, in about
3.6 minutes. Its inner postflight still reported 0.067 I/O-wait core equivalents
and the outer postflight 0.113, while busy CPU was approximately 0.013 cores.
Both attempts remain available for inspecting variation, but are not the final
accepted controlled window. These observations are consistent with pending
artifact writes; they do not prove that every I/O event came from the benchmark.

Before the next attempt, the protocol explicitly drains benchmark artifacts:
`fsync` each regular file and directory in that measurement's output tree, then
apply the same ten-second quiet test. Symlinks and special files are skipped;
no cache is evicted, no input/scientific file is modified, and no global sync
is requested. Flush count, observed file bytes and elapsed time are retained.
This happens outside all measured intervals. Previous attempt artifacts are
also flushed before starting the new window.

The native process timer still measures the same invocation and regular-file
output scope. Durable-storage latency is not newly included in that timer.
Workload-generated I/O during the batch remains part of the observed workflow;
the postflight check now follows completion of the artifact writes. It is not
a timing filter or an outlier-removal rule. Every trial and failed gate remains.

The quiet thresholds, eight alternating pairs, warmups, CPU affinities, input
hashes, native binaries, exact scientific checks and reference tolerances remain
unchanged. All controller/helper hashes are bound to the new receipts. The
[artifact-flush check](artifact_flush_check.json) verifies two regular files,
a directory, symlink/special-file exclusion and unchanged input/content hashes.


## Accepted final windows

On 2026-11-17, window 6 completed all eight sequential stages: CLI A/B, kernel
A/B, kernel profile, CLI profile, Python A/B, Python/native layers, accuracy/cost
and competitive volume. Every required outer and inner quiet endpoint passed
without changing the 0.15 CPU-core/0.05 I/O-core limits. All scientific checks
passed within their declared scopes. Raw window 4/5 distributions remain, and
window 6 alone supplies final accepted whole-campaign timing claims.

Warp UI PID 6817 was stopped only after explicit user authorization, from an
independent transient user service; child shells/agent processes were not
signalled. Both an independent guard watchdog and a separate user-service timer
could restore that exact process identity. The window completed in about 1,799
seconds, with successful resume recorded. The earlier window 2 interruption and
Warp process disappearance remain documented; no unsupported crash cause is
inferred. Window 3's ancestor refusal prevented an unsafe launch.

A final coverage audit added the prospective twelve-case Sterimol protocol. Its
fixed accuracy gate preserved two external failures. All 36 cost trials were
completed in a separate approximately one-minute guarded window, both quiet
checks passed, and Warp resumed. Its restore timer was cancelled after successful
completion. Failed external methods receive no equal-accuracy speed ratio.

Pre-run-hashed main-window controllers were archived byte-for-byte after the
run and checked against the recorded hashes before any final reporting edits.
The archive receipt explicitly records this chronology. Raw service, timer and
guard receipts bind authorization scope, process identity and restoration.
Artifact flushing, like report analysis and plotting, is outside timed intervals.
