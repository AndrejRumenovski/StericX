# C2 implementation against accepted C1b

This addendum activates the historical C2 proposal and test matrix against the
newly accepted C1b baseline. Those pinned historical documents remain unchanged.
No C2 performance result or acceptance is claimed here.

The baseline is the actual C1b build, manifest SHA-256
`def1aa57289814fd24e787dd95092cd0115a67091f250e4be9c49f0f4a5e663d`,
admitted by `accepted_c1b_v1/acceptance.json`. The fresh diagnostic reprofile
manifest has SHA-256
`602fb4e9c954ff0d5331c1aa282e5f0903433bcdcd668e16897965c330f7f0b5`.
All 20 measured screen runs make exactly 1,000 Student-t quantile calls. Their
median fractions of external wall time at 1/2/4/6 threads are
20.2458/20.1884/20.0176/20.0371 percent. The ideal complete-removal ceiling is
approximately 1.254x. Reusing 999 of the 1,000 identical-model values supports a
conditional 15–20% native screen wall reduction hypothesis, with small stack-cell
and branch costs. Acceptance requires a reproducible reduction of at least 10%
beyond paired noise and no unacceptable other-workload or memory regression.

Before production edits, the new ignored archive
`.stericx/profiling/scientifically_validated_optimization/candidate_c2_quantile_reuse/`
retains current source/tests/scripts, Git state/diff, dependency files, this
addendum, the old proposal/matrix, baseline admission/build/profile/workload and
controller identities. A tiny harness links the already built C1b library and
records the unchanged public interval method's `Option` status and exact f64
endpoint bits across valid, unavailable, nonfinite, overflow and unchecked
zero-parameter compatibility cases. The latter is existing direct-API behavior,
not a scientifically endorsed interval; CLI model validation rejects it.

The only production changes are an opaque borrowed interval evaluator in
`model/domain.rs`, its necessary library re-export, and one evaluator per screen
invocation. Its `OnceCell<Option<f64>>` stores only the multiplier. Construction
performs no validation or numerical work. Each candidate still calls the existing
leverage method first. The cell initializes only after that call returns `Some`,
including the existing unchecked direct-API edge. Both cached `None` and cached
finite multiplier bits persist for that borrowed immutable model only.

The existing public `prediction_interval` body stays byte-identical as the
differential reference. The new method duplicates its short arithmetic and finite
guards in the same order. It retains the same `uncertainty` profile scope label;
an additional diagnostic cache-initialization scope counts initialization without
native runtime counters. Quantile equations, leverage calculations, model
validation, confidence intervals, bootstrap uncertainty, applicability, candidate
order, filters, ranking and errors remain unchanged.

Focused tests compare status and endpoint bits to both the unchanged public
method and the pre-edit captured vectors. They distinguish uninitialized,
cached `None` and cached `Some`; test interleaved models and a fresh context after
model changes; and cover invalid-leverage, selected/unused nonfinite features,
signed zero, saturation and interval overflow followed by recovery. Existing and
targeted CLI tests retain model/error precedence and the independent Student-t
and bootstrap channels. Full scientific, reference, CLI and engineering gates
remain required before any native C2 timing.
