# C1: ordered descriptor file parallelism

Status: accepted for descriptor file batches against corrected baseline commit
`1ecfbd5be8b354711bdadd8a8d148417e03d8e7a`. The original pre-edit proposal is
preserved in the candidate freeze. The prediction and acceptance criteria below
are that original plan, not retrospective estimates.

The frozen C1 build passes the complete scientific and engineering gates.
Forty native A/B configurations and five supplemental configurations retain all
970 launches. Paired 10,000-file speedups are 1.931×/3.693×/5.287× at 2/4/6
threads. The decision in `accepted_c1_v1/acceptance.json` explicitly accepts the
roughly 2% parser wall regression, small unrelated CPU costs, and higher batch
RSS/CPU in exchange for these large batch gains. Search wall measurements remain
too variable for a precise claim. See the evolving final performance report;
the acceptance applies to the frozen C1 executable, not unmeasured later edits.

Fresh evidence is `docs/performance/CORRECTED_PROFILE.md` and its linked
`current_corrected_measurements.json`, using `accurate_baseline_v1` and
`corrected_workloads_v2`. Descriptor file work occupies 99.161% of one-thread
and 99.267% of six-thread diagnostic wall time. Native 10,000-file medians at
1/2/4/6 threads are 4991.033/5133.013/4972.785/4958.143 ms, using about one CPU.
Using S=1/(1-f+f/6), ideal six-worker ceilings are approximately 5.758x/5.788x.
These are conditional ceilings. The proposal expects a realistic 2–4x at four
or six workers and neutral one-thread performance; neither is an observed result.

Only `descriptors_command` scheduling changes. An indexed Rayon slice iterator
computes one unchanged `descriptors_for_file_with_reference` call per original
input position. Each file retains its parser, authoritative topology, conformer
order, checked arithmetic, geometry kernels and reductions. Collecting into an
indexed vector preserves duplicates and original input order. The existing main
thread records outcomes and emits errors/results serially in that order.

The donor/reference batch-index guard remains before any file processing or pool
construction. Single-file invocations and explicit one-thread settings keep the
serial path. A private Rayon pool honors its existing thread-count configuration;
if pool construction fails, computation falls back to the original serial path
without adding an error or warning. A pool with one worker also uses the serial
calculation branch. No conformer or within-kernel parallel reduction is added.

Exactness follows from independent immutable per-file inputs, unchanged local
calculations, indexed collection and unchanged serial emitters. Completion times
and file-read order can differ; fixed immutable input files are the benchmark and
scientific contract. No shared numeric accumulator is introduced. Worker-local
geometry plus the temporary ordered outcome vector increase peak memory; measure
RSS/allocation changes rather than assume they are negligible.

Existing batch tests cover JSON/text/CSV, differently sized files, successes,
errors, duplicates, terminal all-error behavior and 1/2/4/6 threads. Meaningful
extensions cover both explicit-index guards and an SDF conformer ensemble in a
mixed ordered batch, including corrected finite means. Focused tests precede the
root-owned full scientific/CLI/reference/engineering gates. No native benchmark
or release build is part of this implementation step.

Acceptance requires exact scientific and observable CLI equivalence to the
corrected freeze, retained independent-reference classifications, complete
engineering checks, then alternating paired native measurements on all ten frozen
workloads at all four thread settings. The target is reproducible at least 10%
reduction in batch wall time beyond paired noise, without a material regression
in other workloads, one-thread behavior or memory cost. Preserve every sample,
error and attempt. No tolerance change, outlier exclusion or historical profile
substitution is permitted.

Before editing, the private candidate evidence directory captures this proposal,
all original `src` files, Cargo manifests, the existing batch test, full Git diff,
status, HEAD, corrected profile and baseline manifest with SHA-256 identities.
