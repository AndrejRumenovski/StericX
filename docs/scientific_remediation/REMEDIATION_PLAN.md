# Scientific remediation before optimization

The user authorized scientific corrections before restarting the performance
experiment. The starting source is audit publication commit
`b515c4f1736a1519be19a3ecdf40c3c1327ac9b4`. Its documented failures are observations
to correct, not numerical expectations to preserve in a later optimization oracle.
This plan is an in-progress record, not a declaration that the system is validated.

The original [audit](../scientific_accuracy_audit/SCIENTIFIC_ACCURACY_AUDIT.md),
reference programs, input structures, failed cases and frozen outputs remain
unchanged. The [ledger](remediation_ledger.json) preserves all 43 negative or
terminology propositions separately, with their original evidence and limitations.
Multiple propositions can refer to the same defect. Some deliberately test stronger
compatibility hypotheses than the implementation claims.

## Correction policy

For each implementation change, record the original witness, independent expected
behavior, mathematical/source justification and targeted validation before accepting
the change. A corrected scientific result can intentionally differ from the audited
build; that difference must be explained. Exact equality to the old build is not a
scientific gate during remediation.

Retain legitimate, explicit conventions unless independent evidence shows that the
implementation violates them. In particular, do not change the grid, radii or
coordinate convention merely to improve agreement with Morfeus or Kraken. Published
reference disagreement, empirical interval miscalibration and unavailable historical
thermochemistry cannot be repaired by silently changing reference expectations.
Correct overstated claims and retain unresolved evidence limitations.

## Work allocation and acceptance

| Area | Demonstrated issues to investigate and correct | Independent acceptance evidence |
| --- | --- | --- |
| Geometry | Symmetric-volume rejection; atom-order-dependent automatic axis; near-axis alignment; invalid/overflow/degenerate geometry; explicit topology versus proximity inference | Preserved synthetic and real conformers, analytic support formulas, matched-convention Morfeus, topology evidence, rigid transformations and permutations |
| Thermodynamics and kinetics | CREST temperature/populations; invalid and scale-sensitive weights; nonfinite aggregates; energy span; ee domain; absolute barrier versus barrier difference; reaction-temperature provenance | Published equations, official CREST/RDKit semantics, manually calculable ensembles, exact failing API/CLI inputs |
| Models and screening | Nonfinite evaluation; Student-t tail saturation; singular-domain geometry; screening aggregation mismatch; source classification boundaries; misleading uncertainty/trust terminology | Analytic Cauchy tails, independent linear algebra/statistics, explicit held-out-data limitations, source tables and consistent feature construction |
| Scientific descriptions | Kraken center/radius/grid compatibility; native versus Python precision; validation scope; formula versus structure identity; conditional versus complete-workflow validation | Original Kraken source/SI and independently frozen audit comparisons; historical numbers remain labeled historical |

Run all Rust/Python tests, Clippy, rustdoc and formatting, then the complete available
scientific corpus against unchanged independent definitions. Report every changed
result, still-failing case and unresolved claim. New tests supplement independent
evidence; they do not replace it.

Only after this recheck, freeze the corrected source, executable, datasets, outputs,
errors and reference results as a new optimization baseline. Profile it from scratch.
Every subsequent optimization must preserve that baseline's scientific and behavioral
contract, pass independent reference checks, and improve native runtime reproducibly.

## Earlier performance experiment

The unaccepted file-parallel descriptor experiment was withdrawn when the user
selected remediation. Its patch, source snapshot and partial verification record are
preserved under
`.stericx/profiling/scientifically_validated_optimization/candidate_01_file_parallel/`.
It was not benchmarked and has no accepted speedup claim. Production was restored
before remediation began. Earlier preparation documents are archived under
`history/pre_remediation_preparation/` in that profiling root, with hashes checked
against `preparation_checkpoint.json` before copying.
