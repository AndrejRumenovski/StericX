# Current optimization baseline: scientific remediation in progress

The performance task is **not complete**. No production optimization has been
accepted. A file-parallel experiment was withdrawn before benchmarking and its
source change restored; its patch and partial checks are preserved separately.

The current commit is `b515c4f1736a1519be19a3ecdf40c3c1327ac9b4`, which publishes
the independent scientific audit. Its 42 Cargo/Rust source files match the
audited implementation, and a fresh locked release build reproduces executable
SHA-256 `b577c49d3e98b4745213fce0e70e93a55fc61120c5abe6fdb0658547c43cab94`.
That commit is the historical audited baseline, not the corrected optimization
baseline. Remediation is now authorized and in progress in the working tree.

The [audit](../scientific_accuracy_audit/SCIENTIFIC_ACCURACY_AUDIT.md) supports
many scoped numerical claims, but does not establish the unrestricted passing
baseline requested by the new optimization task. For example:

| Required property | Current audited evidence |
| --- | --- |
| Valid symmetric geometries remain accepted | Valid PH3 occupied volume is rejected for zero quadrant asymmetry (G22). |
| Automatic axis selection is invariant to atom ordering | A tied-axis witness changes L/B5 by 2.57 Å (G06). |
| Requested thermodynamic temperature is honored | Parsed CREST populations ignore the requested temperature (C09). |
| Validation rejects invalid predictions | `evaluate` accepts NaN predictions with exit status zero (M69). |

Claim IDs refer to the individual propositions and limitations in
[CLAIMS.md](../scientific_accuracy_audit/CLAIMS.md); these are not four exhaustive
categories or a count of independent defects. Preserving known failures exactly
is an optimization regression check, not evidence that the failures are correct.

The user explicitly chose: **correct the demonstrated scientific failures first,
then freeze a new optimization baseline and continue the goal**. Follow the
[remediation plan](../scientific_remediation/REMEDIATION_PLAN.md). Independent
definitions remain fixed. Intentional corrections to the old scientific outputs
must be justified and recorded before a new optimization oracle is frozen.
Historical audit artifacts remain unchanged.

Preparation proceeds independently: the current audited system is frozen under
`.stericx/profiling/scientifically_validated_optimization/current_audited_b515c4f/`.
The directory name follows the requested task layout; its machine-readable
identity explicitly sets `scientifically_corrected_or_validated_baseline` to
`false`. All 74,298 sealed audit payload files, totaling 1,150,324,652 bytes, were
verified before the freeze. The original audit remains the immutable source of
the scientific observations and negative findings.

Fresh native measurements cover all ten requested workloads at 1/2/4/6 threads,
with one warmup and seven measured launches per configuration. These measurements
are a new performance baseline, not an optimization speedup. Separate diagnostic
builds measure attribution and allocations; they are not used for native speed
claims. Raw timings, resource measurements, output fingerprints and executable
hashes are retained under the new snapshot and its bound run directories.

Three new mechanisms close gaps in the historical regression workflow:

- `scripts/check_current_scientific_equivalence.py` links the selected live Rust
  source, including the private buried-volume observation module. It covers the
  full 31,721-conformer corpus and an additional stream exposing every regional
  bin, alongside geometry, kinetics, aggregation and model-domain witnesses.
- `scripts/check_scientific_cli_equivalence.py` freezes 93 CLI cases plus live
  Python weighting observations, preserving complete records, errors and
  candidate decks. Only explicitly named process metrics and the portable
  model's creation timestamp are excluded from exact comparison.
- Independent scientific reference replay is prepared separately from baseline
  equality, using the original equations, reference data and classifications.

Initial CLI harness failures are retained: preparation first detected an
incorrect assumed command count (33 model commands exist); the first full
fidelity run identified the previously omitted `simulate` timing field
`total_microseconds`. All other compared scientific values, artifacts and errors
matched. A new manifest version explicitly identifies this process timing field;
the earlier evidence was not rewritten. Neither event is a scientific tolerance
adjustment or a production optimization.

The requested [optimization report](SCIENTIFICALLY_EXACT_OPTIMIZATION.md) currently
contains an explicitly incomplete checkpoint. It must not be presented as a final
success report before the baseline prerequisite and both scientific and
performance acceptance gates have been resolved.
