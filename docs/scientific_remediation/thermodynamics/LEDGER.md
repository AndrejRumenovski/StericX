# Thermodynamics and kinetics remediation ledger

Prepared before production changes. The sealed audit at
`docs/scientific_accuracy_audit/` remains immutable. Its C01–C21 classifications
describe the audited implementation; this ledger records subsequent corrections
and does not replace the original findings. No performance optimization belongs
to this work.

## Evidence and conventions

The starting evidence is the complete C01–C21 claim set, the original
`kinetics/KINETICS_CONFORMERS.md`, all row-level kinetic/conformer results and
their frozen inputs/observations. `historical_evidence.json` binds those files by
SHA-256 before changes. Reference equations are evaluated independently of the
production routines, using exact SI constants and high-precision Decimal where
appropriate.

Primary sources read:

- [IUPAC Green Book](https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf),
  printed pp. 45–46 and 65–67; the audit retains the PDF and extracted text.
- [NIST CODATA constants](https://physics.nist.gov/cuu/Constants/Table/allascii.txt):
  exact SI definitions of kB, h and NA, with 4184 J per thermochemical kcal.
- [CREST keyword documentation](https://crest-lab.github.io/crest-docs/page/documentation/keywords.html):
  `--temp` sets the CREGEN population temperature. The frozen CREST v2.12
  `cregen.f90` obtains T from `env%tboltz` and sums individual rotamer
  populations into conformer-set populations. A conformer degeneracy is not
  permission to assume distinct rotamers have identical energies.
- [RDKit force-field API](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html):
  MMFF status 0 means converged; status 1 records failure to converge within the
  iteration limit. The retained optimizer status must remain visible.

For a first-order elementary process with transmission coefficient one,
`ln k = ln(kB T/h) - ΔG‡/(RT)` requires an **absolute activation free energy**.
Selectivity instead uses `ΔΔG‡ = G‡S - G‡R` and equal pathway prefactors. Positive
ΔΔG favors R. Absolute ee is `100 |pR-pS|`; ee-derived energy magnitudes discard
enantiomer identity. State populations are normalized nonnegative weights;
electronic or force-field energies are approximations to free energies.

## Decisions and retained witnesses

| Claims | Sealed witness and independent expectation | Planned correction or retained limit |
|---|---|---|
| C01–C04, C06 | `kinetics/results/kinetics.csv`; same-selectivity barrier pairs `(10,11)` and `(20,21)` have ratio 5.4076065654 but different absolute rates. Tiny nonzero selectivity is representable even when subtraction of rounded percentages cancels it. | Preserve sign/units. Use a log-rate representation and an explicit checked representable-rate API; document legacy f32 underflow/overflow. Compute tiny ee and minority populations without subtracting two nearly equal rounded numbers. Invalid temperature/energy remains invalid. |
| C05 | `kinetics/frozen_outputs/cli_runs.json`, `simulate --ddg 1 --temp 298.15`; a difference cannot identify either absolute barrier. | Remove the absolute-rate output from the difference-only CLI. State its signed difference convention and assumptions. No inferred absolute barrier or new convention. |
| C07–C08 | `kinetics/results/mmff_weights.csv`: optimizer status 1 is retained, negative status/nonfinite energy excluded. | Retain explicit nonconvergence policy and statuses; clarify that retained conformers are not all proven converged. Validate direct weighting inputs. Do not silently discard additional states. |
| C09–C11 | `kinetics/inputs/crest_default.log` at 298.15 K, requested 500 K: unit-degeneracy populations must be `(0.7323226909,0.2676773091)`. `crest_negative.log` contains negative populations. | Pass configured `--temp`; validate table populations/temperature and energy/degeneracy provenance. Reweight only when the available state information justifies it; never silently reuse a wrong-temperature table or invent equal-energy rotamers. Reject malformed population tables. Corrected cache identities must not reuse old erroneous populations. |
| C12–C14 | `parse_provided.csv` weights 1:3 gives analytic Sterimol `(5.2,1.7,4.97)` Å; missing weights use a uniform mean; bad conformers reject the row. | Keep explicit uniform fallback and complete-ensemble rejection. Supplied populations are not automatically claimed to be energy-derived. |
| C15–C17 | `kinetics/inputs/aggregation.json`, `parse_tiny_weights.csv`: scaling 1:3 by a positive finite common factor must preserve normalized means; finite large sums must not overflow; NaN descriptor inputs must fail. | Normalize after maximum-weight scaling with f64 accumulation; reject negative/nonfinite/all-zero weights and nonfinite scalar inputs. Preserve extrema and the identity-associated minimum-volume conformer. Geometry ownership delegates only `aggregate()` and a separate regression-test module to this work. |
| C18 | `parse_relative_offset.csv` supplies energies `[1,2]`; its span is 1, not 2. | Compute max-minus-min after existing count/domain validation. Do not require the supplied minimum to be zero. |
| C19–C20 | `kinetics/inputs/ee.json`: finite `|ee|<100` follows the log-ratio equation; exactly 100 has an infinite limit; greater than 100 is invalid. | Reject censored/out-of-domain/nonfinite measured values instead of clipping to an invented finite target. Missing values may remain explicitly missing in dataset normalization. Validate positive finite temperature. Preserve magnitude convention. |
| C21 | `kinetics/results/ni_hda_target_temperature.csv` and corrected primary SI: experimental temperature is 80 °C = 353.15 K. ID 2064 has a separate source inconsistency. | Correct future Ni-hDA normalization metadata to 353.15 K while retaining published target values and flagging unresolved source inconsistency. Root owns historical dataset/Study provenance updates; this work does not rewrite historical tables or targets. |

## Verification plan

Retain the original failing observations. Add focused Python and Rust tests for
the above witnesses, invalid values, zero/tiny/large weights, common energy
offsets, degeneracy handling, temperature changes, sign reversal, and numeric
range behavior. Derive reference values from equations or analytic descriptors,
not from candidate outputs. Run related existing suites and new tests. Write a
separate result report with corrected behavior, evidence paths and remaining
method limitations; do not relabel an untested correction as validated.

MMFF/CREST sampling completeness, electronic-energy populations as approximations
to free-energy populations, absolute enantiomer assignment from unsigned ee,
and Eyring applicability to a multistep catalytic mechanism remain limitations.
