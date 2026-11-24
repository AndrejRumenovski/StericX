# Remediation measurement adapter

This is a SUT adapter, never an independent reference. `geometry` preserves the
original request meaning (including inferred BV connectivity); checked API
failures are explicit `error` objects rather than numerical zero sentinels.
`geometry_topology` is an additional operation using the supplied three
`neighbors` for center/BV as well as pyramidalization. IEEE bit fields are
retained. Per-plane private observations retain their original layout.

Build with `scripts/replay_scientific_remediation.py`: it supplies
`src/observed_buried_volume.rs` as the captured current source followed by the
byte-identical historical private append. `topology_append.rs` adds only the
explicit-topology observation path and does not alter historical operations.
