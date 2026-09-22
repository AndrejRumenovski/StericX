# Documented-command check against corrected build v3

All 17 native commands extracted from the screening tutorial and model-format
producing section exited 0. The two parse/fit/screen workflows used separate new
output roots under `.stericx/scientific_remediation/demo/documented_commands_v2`.
The receipt preserves exact documented/executed arguments, archived document
hashes, native/build/input bindings, all logs and 17 generated output hashes.
Only the executable and destination root were substituted during execution.
Preparation was already verified by the immutable metadata-v2 input receipt.

The S001 tutorial fit selects `nbo_charge` alone (eight train, three blind), so
its screen consumes supplied electronics without recomputing unused sterics.
The full-record example selects `B5_x_nbo_charge` (ten train, one blind), and its
screen computes the declared means from all supplied conformers. Seven schema
cases also gave their expected outcomes: legacy 1/2 inspect as unknown; their
geometry substitutions fail; missing schema-3 aggregation, future schema and
malformed training geometry are rejected.

This verifies command/schema behavior for the specifically pinned v3 build, not
empirical predictive validity or a later snapshot's unchanged behavior. Interim
v2-build command results are retained separately at `documented_commands_v1`.
