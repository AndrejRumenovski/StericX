# Fresh scientific observations and reference replay

This is preparation for the optimization task, not a post-optimization pass.
Production remains at `b515c4f1736a1519be19a3ecdf40c3c1327ac9b4`. Its scientific
limitations remain those of the published independent audit. See
[baseline status](BASELINE_STATUS.md) for the unresolved prerequisite.

## New observations of the current source

The live Rust observer replaces both historical frozen dependencies: its library
links current production code, and its private buried-volume observation module
contains the exact current source followed by the unchanged audit adapter.
There is no old scientific implementation hidden behind that private adapter.

The first complete observation covers eleven lanes: 844 geometry requests,
15 focused cases, four alignment cases, four rotation cases, eleven numerical
cases, twenty historical-frame cases, 31,721 Kraken conformers, 195 kinetic cases,
eleven aggregation cases, 33 model-domain cases, and a second full 31,721-conformer
stream exposing every buried-volume region. Original public observations,
failures, ordering and exact floating-point bits match the audited observations.
The additional Kraken stream exposes 95,139 frames and 1,141,668 quadrant/octant
values. Eight known connectivity failures prevent private frames; a ninth
conformer has a public descriptor failure while retaining its private frame data.

The strengthened CLI/Python v5 oracle captures 93 CLI commands and live Python
thermodynamic fixtures at each of 1, 2, 4 and 6 threads. Eighty-two CLI cases have
historical anchors; nine add batch contracts derived independently from frozen
single-file observations and the reviewed output/error contract; two add
candidate-deck exports. Historical scientific stdout,
stderr, status and generated records match at all four settings. Only the named
process timing/RSS fields and portable-model creation time are excluded from
scientific comparison; raw bytes are retained. Deck rows are checked against the
ranked JSON hits, including predictions, interval endpoints, uncertainty and
domain classifications, and the sidecar checksum/provenance must match.

## Fresh independent calculations

The reference harness runs the audit's unchanged independent equations and
Morfeus code against these new observations. It verifies pinned package/source
identities, stages immutable scientific inputs, and refuses to substitute an old
native executable. Reference results and historical differences are stored in new
directories. The audit's original outputs and claim classifications are unchanged.

| Campaign | Fresh evidence | Result against historical scientific values |
| --- | --- | --- |
| Geometry | 829 reference cases; 27,542 descriptor and 58,176 regional comparison rows; 238 invariance rows | All numerical comparisons exact; one exception traceback differs in formatting |
| Kraken | Full 31,721-conformer analysis and Morfeus/analytic references; 590 initial and 2,009 regional-outlier reference records | All individual numerical rows exact; one aggregate R² differs by one ULP |
| Thermodynamics and Eyring | Nine result artifacts, including 753 kinetic comparison records, Python weights, aggregation and packed ensembles | Exact values and failures |
| Models/statistics/screening | Nine independent scripts, 35 result artifacts, 72 reaction fits, 4,476 classifier predictions, 33 domain cases | Exact values and failures |

The geometry/Kraken exact comparison correctly returned a nonzero status for its
two differences. That result has **not** been changed to a clean pass:

1. `asymmetric_toy__translation_10000000.0` retains the same
   `ValueError("undefined direction")` and all other fields. The traceback quotes
   different source formatting and line numbers. The sealed reference source
   and its newly staged copy are byte-identical; the old captured traceback
   predates that formatting, as inferred from the text it quotes.
2. Full-primary B1 R² is `0.977389496618067` in the new replay versus
   `0.9773894966180671` historically: an absolute difference of
   `1.1102230246251565e-16`. The complete underlying descriptor CSV and
   Morfeus/analytic reference stream are byte-identical. On those identical
   31,713 valid rows, controlled BLAS reductions reproduce the new number with one
   thread and the
   historical number with 3/4/6/8/16 threads. An 80-digit Decimal calculation
   rounds to the historical f64 value. This isolates reference aggregation
   rounding; it does not establish the historical process's thread count.

No tolerance has been added. The exact mismatch receipts and a reproducible
bounded diagnosis remain beside the new reference results.

For a future candidate, historical replay must remain separate from an exact
comparison against this fresh independent baseline. Both complete campaigns must
have the same 84-result inventory (40 geometry/Kraken, nine kinetics and 35 model
artifacts), records, IDs, ordering, errors and schema. Source/build/executable,
input, capture, archived helper, package, Morfeus and loaded BLAS identities must
be verified. The reference wrapper pins BLAS/OMP to one thread; candidate admission
must additionally verify the loaded threadpool configuration. Candidate equations
must consume candidate observations and candidate-derived Kraken analysis.

Every reference result, residual, classification and aggregate must match this
fresh reference baseline exactly. In particular, candidate B1 R² must equal
`0.977389496618067`; either of two numbers is not an allowed range. The candidate's
historical-difference inventory must also match precisely the two diagnosed
locations, values and normalized traceback strings above, with no extra or
missing differences. Preserve the historical nonzero result. An environmental or
unexplained difference requires investigation, not a new tolerance or a generic
exception for tracebacks or one-ULP changes. This rule establishes independent
regression equivalence within the declared coverage, not global correctness.

## Oracle review and failed preparation attempts

Independent review found ways a future broken capture could falsely pass:
missing new deck exports in both builds, failed private-bin dumps in both builds,
an empty lane dictionary carrying stale full-coverage flags, and insufficient
binding of stored metadata to current manifest/executable bytes. The completed
captures contained the actual required outputs; these were guard defects, not
evidence of missing data in those captures.

The CLI v4 harness freezes its own source and each observed executable/Python
source, binds the actual manifest, rechecks complete raw observations, enforces
required nonempty artifacts and checks deck consistency. Adversarial tests reject
changed manifests, binaries, helpers, missing outputs, changed raw metrics,
one-bit packed differences and one-ULP scientific text differences. Old versions
and receipts are preserved. CLI v5 adds shuffled mixed-success/failure batches,
duplicates, all-failure batches and the multi-input donor-index precheck in all
three output formats. All 94 command/Python fingerprints agree across the four
thread settings; all historical and independently derived contracts pass.
The Rust v3 observer now requires actual nonempty lane sets, recomputes coverage
and historical fidelity from raw files, binds source/build/helper/executable
identities, and rejects malformed private frames or unreviewed dump failures.
The fixed full-Kraken lane requires finite regional values and positive grids.
Reviewed numerical nonfinite fixtures retain their original contract. All eleven
lanes of the v3 full replay pass exact historical fidelity, including every
available regional bin; v1/v2 evidence remains untouched.

A separate cross-version bridge verifies all 64,579 responses across the eleven
lanes: v1/v3 request bytes, stdout, stderr and exit status are identical. It also
verifies all 42 scientific sources against the immutable v3 snapshot and sealed
audit, and verifies the completed independent-reference result inventories and
their v1 observation links. Thus the stronger v3 capture is connected to the
already completed independent calculations without substituting an old kernel or
silently rerunning changed equations. These response counts include both public
and private replay of the same Kraken conformers; they are not 64,579 distinct
molecules.

Other retained preparation failures include the initial CLI count assumption
(33 model commands exist), an omitted `simulate` process timing field, a missing
Kraken percentage projection in model staging, and a subprocess guard correctly
blocking lazy platform metadata collection. None required a production change or
a changed scientific expectation.

## Engineering checks

The unchanged production baseline passes all 264 Rust tests, the separately run
published-screening regression, Clippy with warnings denied, rustdoc with warnings
denied and Rust formatting. The warning-denial rustdoc environment is additionally
recorded as an explicit command in `baseline_rustdoc_explicit_environment/`.
The first complete Python discovery run passes
54 tests, including the new CLI oracle tests. The final preparation suite passes
71 tests after the observer and batching checks are added. Python lint/format
also pass with preexisting, untracked `docs/media/` drafts excluded and untouched.
Study 011's verification command also passes.
These checks demonstrate software regression status, not scientific correctness.

## Scope and limitations

Independent re-evaluation is separate from baseline equality. Neither upgrades
an `INCORRECT`, `UNCERTAIN` or `OUT OF SCOPE` claim. Experimental generalization,
ensemble provenance, chemical donor justification and the original literature
limitations remain unresolved where the audit says so.

Focused/alignment/historical geometry supplements, Kraken focused rounding
diagnostics and formula-identity helper observations remain historical independent
evidence. Their current StericX observations are replayed where included above;
this reference pass does not recalculate every supplemental independent analysis.
Study 001/002/003 fits and Study 011 rankings are reanalyzed historical evidence,
not newly measured chemistry or newly fitted native models. Independent screening
equations use the final v2 audit campaign; older v1 captures are retained without
another separate reference calculation. Full-corpus regional observations compare
against current-baseline bytes; independent regional quadrature covers the audit's
existing campaigns, not every region of every conformer.

No post-optimization scientific claim is possible before a candidate exists and
passes these gates plus the required independent follow-up.

## Evidence locations

All paths below are relative to
`.stericx/profiling/scientifically_validated_optimization/`:

- `scientific_oracle/`, `scientific_observer_baseline/`,
  `scientific_observations_baseline/`: first complete current-source replay.
- `scientific_oracle_v3/`, `scientific_observer_baseline_v3/`,
  `scientific_observations_baseline_v3/`: final strengthened baseline replay with
  immutable source snapshots and complete coverage/fidelity gates.
- `scientific_oracle_bridge_v1_v3/`: exact cross-version stream/source/reference
  bridge and final admission self-check.
- `scientific_oracle_hardening_v3_tests/`: ten focused adversarial tests and lint
  receipts for the final observer harness.
- `cli_oracle_v4/`: first strengthened four-thread-setting CLI/Python captures.
- `cli_oracle_v5/`: those cases plus nine independently derived descriptor batch
  contracts, covering 93 CLI cases plus Python at all four settings.
- `independent_baseline_geometry_kraken/`: complete fresh reference campaign and
  its two unchanged exact-comparison failures.
- `independent_baseline_kinetics_check/`: fresh analytical thermodynamic/kinetic
  comparison.
- `independent_baseline_models_check_v3/`: fresh model/domain/reaction references;
  earlier staging failures are retained in separate siblings.
- `independent_reference_difference_diagnosis/`: source/data hashes, controlled
  BLAS reductions and Decimal calculation explaining the two differences.
- `baseline_engineering_checks/`: exact commands, statuses and raw logs.
- `baseline_rustdoc_explicit_environment/`: explicit warning-denial rustdoc command.
- `final_preparation_engineering_checks/`: 71-test Python discovery, lint/format
  results and captured helper/test source identities.

The executable/source/data freeze and measurement seals are documented in
[the fresh profile](CURRENT_AUDITED_PROFILE.md). All replay helpers refuse to
overwrite a completed output directory.
