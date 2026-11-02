# Output-storage investigation

Specified before the second timing campaign, 2026-11-02 UTC.

The first complete 504-sample campaign triggered four latency investigations.
Median child CPU times changed little, whereas short-process wall times sometimes
varied by tens of milliseconds. Sixteen separate `strace -f -ttt -T` diagnostic
runs locate the dominant delay in the same 678-byte `write(1, ...)`: as much as
84.160 ms for baseline and 14.711 ms for candidate in this diagnostic sample.
Tracing is intrusive and these durations are not performance samples.

System I/O pressure was substantial. Other user workloads were also observed;
their processes were neither stopped nor modified. Consequently the original
wall-time results cannot establish an isolated engine regression or speedup.
The original data, all triggers and both favorable/unfavorable ratios remain.

Repeat the complete seven-workload, four-thread-count, eight-pair protocol using
regular output files under `/tmp`, which is tmpfs on this host. Keep executable
bytes, all frozen input paths/geometries, descriptor settings, affinities,
warmups, order and acceptance checks unchanged. The only intended experimental
change is the output filesystem. Exclude build/reference/download jobs during
timing. Record host pressure/process snapshots before and after.

This measures the engine with RAM-backed regular-file output, including normal
process startup and serialization. It is not a claim about sustained disk export
speed or every deployment filesystem. Do not replace or pool the original
storage-sensitive samples, discard outliers, or call a changed output location
a production-code optimization. Preserve raw output and metrics outside `/tmp`
after completion with a relocation receipt so reboot does not destroy evidence.
