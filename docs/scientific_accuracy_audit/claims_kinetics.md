# Conformer and kinetics claim inventory

Numerical evidence and failed cases: [KINETICS_CONFORMERS.md](kinetics/KINETICS_CONFORMERS.md). Each entry is classified separately. These entries include documented behavior and explicitly identified audit propositions. A tested stronger interpretation is not represented as a quotation or a promise that StericX actually makes.

## C01 — Eyring rate arithmetic implements the activation-free-energy equation.

1. **StericX claim / audited proposition:** Eyring rate arithmetic implements the activation-free-energy equation. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** k=κkBT/h exp(-ΔG‡/RT); κ=1; absolute ΔG in kcal/mol, T inK; first-order units.
4. **Implementation location:** `src/kinetics/eyring.rs::calculate_rate_constant`
5. **Independent validation method:** 75-digit Decimal equation on frozen f32 inputs over temperature/barrier grid.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Normal-rate subset maximum relative error 6.9231e-8; subnormal and extreme failures separately retained.
7. **Confidence:** High
8. **Limitations:** Equation correctness does not identify a physical barrier or establish transition-state assumptions.

## C02 — Physical constants and kcal/Hartree conversion are accurate.

1. **StericX claim / audited proposition:** Physical constants and kcal/Hartree conversion are accurate. 
2. **Primary/reference source:** NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** kB/h/NA exact SI; R=kBNA/4184; Hartree energy fromCODATA.
4. **Implementation location:** `src/kinetics/eyring.rs constants; scripts/stericx_quantum.py constants`
5. **Independent validation method:** Independent Decimal arithmetic from primary constants.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — R truncation relative difference ≈3.22e-10; Hartree conversion difference ≈2.03e-10 kcal mol^-1 Hartree^-1.
7. **Confidence:** High
8. **Limitations:** Numerical effects generally smaller than f32 rounding; physical energy-model errors not tested by constants.

## C03 — Zero ΔΔG gives equal R:S and reversing sign exchanges preference.

1. **StericX claim / audited proposition:** Zero ΔΔG gives equal R:S and reversing sign exchanges preference. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** ΔΔG=G‡S−G‡R, common reactant and equal pathway prefactors; zero→50:50.
4. **Implementation location:** `src/kinetics/eyring.rs::product_ratio`
5. **Independent validation method:** Independent two-pathway rates with common energy shift, signed grid and signed zero.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — R/S max error 3.7596e-6 percentage points; zero and sign symmetry observed.
7. **Confidence:** High
8. **Limitations:** This algebraic convention cannot infer the absolute enantiomer from unsigned experimental ee.

## C04 — Enantiomeric excess arithmetic follows the two-pathway populations.

1. **StericX claim / audited proposition:** Enantiomeric excess arithmetic follows the two-pathway populations. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Absolute ee=100|pR−pS|; signed ee, if requested, must retain configuration convention.
4. **Implementation location:** `src/kinetics/eyring.rs::calculate_enantiomeric_excess`
5. **Independent validation method:** Independent pathway probabilities, zero and small/large differences.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Max absolute error 7.5191e-6 percentage points; tiny ee and minority populations can round to zero.
7. **Confidence:** High
8. **Limitations:** Finite precision saturation is not experimentally complete selectivity.

## C05 — simulate can output an absolute rate from the same ΔΔG used for selectivity.

1. **StericX claim / audited proposition:** simulate can output an absolute rate from the same ΔΔG used for selectivity. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** A difference fixes a rate ratio, not either absolute activation barrier.
4. **Implementation location:** `src/commands/simulate.rs; README simulate example`
5. **Independent validation method:** Two barrier pairs(10,11) and(20,21) have same difference but different rates; frozen CLI observation.
6. **Result:** **INCORRECT** — Same ratio 5.4076065654; R rates 290543.5244 versus 0.01358815005 s^-1.
7. **Confidence:** High
8. **Limitations:** The standalone rate function is correct when supplied an actual absolute barrier; the interface is the problem.

## C06 — Rate and selectivity calculations have a defined finite numerical domain.

1. **StericX claim / audited proposition:** Rate and selectivity calculations have a defined finite numerical domain. Audit of the accepted numerical domain, not an assertion that documentation promises finite outputs for every finite input.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; IEEE754 f32 range
3. **Exact convention expected:** Finite input can map outside representable output range; this must be reported.
4. **Implementation location:** `src/kinetics/eyring.rs; src/commands/simulate.rs`
5. **Independent validation method:** Large positive/negative barriers, temperature extremes, NaN/infinity and zero-temperature inputs.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Wider grid includes 16 rate underflows and 15 nonfinite rates; CLI checks inputs only.
7. **Confidence:** High
8. **Limitations:** This is a bounded precision/validation claim, not evidence that ordinary barriers have chemically material errors.

## C07 — MMFF ensemble weights use the normalized Boltzmann equation.

1. **StericX claim / audited proposition:** MMFF ensemble weights use the normalized Boltzmann equation. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** wi∝exp(-(Ei−Emin)/(RT)); kcal/mol; positive temperature; unit degeneracy.
4. **Implementation location:** `scripts/prepare_data.py::embed_and_optimize`
5. **Independent validation method:** Controlled energy/status injection into frozen routine; independent Decimal weights for single/equal/two/three/offset/extreme/temperature/window cases.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum population error 1.9429e-16 in isolated weighting tests.
7. **Confidence:** High for arithmetic
8. **Limitations:** Mocked optimizer results do not validate search completeness or MMFF energies as thermodynamic free energies.

## C08 — Calling the retained MMFF ensemble fully optimized implies that its conformers converged.

1. **StericX claim / audited proposition:** Calling the retained MMFF ensemble fully optimized implies that its conformers converged. Audit of a possible stronger reading of “optimized”; the frozen code retains optimizer status and does not explicitly promise convergence of every retained state.
2. **Primary/reference source:** RDKit installed2026.3.4 MMFFOptimizeMoleculeConfs docstring and official https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html
3. **Exact convention expected:** Return status 0 means converged; finite nonzero status may still be unconverged.
4. **Implementation location:** `scripts/prepare_data.py::embed_and_optimize`
5. **Independent validation method:** Inject status −1/0/1 and nonfinite energies; inspect retained IDs/statuses.
6. **Result:** **INCORRECT** — Status 1 is retained and weighted; negative status and nonfinite energy are excluded.
7. **Confidence:** High
8. **Limitations:** Do not confuse retaining an unconverged state with a proven material descriptor error; magnitude depends on that geometry.

## C09 — CREST populations correspond to configured QuantumConfig.temperature_k.

1. **StericX claim / audited proposition:** CREST populations correspond to configured QuantumConfig.temperature_k. 
2. **Primary/reference source:** CREST official --temp documentation and frozen v2.12 cregen.f90: https://crest-lab.github.io/crest-docs/page/documentation/keywords.html; IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Population temperature must equal requested T; external table temperature cannot be assumed.
4. **Implementation location:** `scripts/stericx_quantum.py::crest_ensemble and _conformer_thermodynamics`
5. **Independent validation method:** Frozen realistic two-state summary at 298.15 K, requested 500 K; independently evaluate weights.
6. **Result:** **INCORRECT** — Population error 0.1116127958 at 500 K; command omits --temp and parsed table bypasses recomputation.
7. **Confidence:** High
8. **Limitations:** Direct frozen helper/command audit; no actual external CREST execution is claimed.

## C10 — CREST missing-table fallback correctly weights supplied electronic energies.

1. **StericX claim / audited proposition:** CREST missing-table fallback correctly weights supplied electronic energies. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** Convert Hartree differences to kcal/mol, subtract minimum, normalize atconfigured T.
4. **Implementation location:** `scripts/stericx_quantum.py::_conformer_thermodynamics`
5. **Independent validation method:** Frozen two-frame energy input at 298.15/500 K, independent Hartree conversion.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum population error≈1.03e-12, including coordinate-comment energy rounding.
7. **Confidence:** High for tested arithmetic
8. **Limitations:** Degeneracy defaults to one and electronic energies are not state free energies; missing/incorrect degeneracies change expected populations.

## C11 — Parsed CREST weights are valid nonnegative probabilities.

1. **StericX claim / audited proposition:** Parsed CREST weights are valid nonnegative probabilities. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Each population ≥0; positive finite normalization total.
4. **Implementation location:** `scripts/stericx_quantum.py::_conformer_thermodynamics`
5. **Independent validation method:** Synthetic malformed population table with weights −0.5,1.5.
6. **Result:** **INCORRECT** — Helper accepts negative population because only finiteness and sum>0 checked.
7. **Confidence:** High
8. **Limitations:** Malformed-table validation test; does not assert real CREST normally emits negative populations.

## C12 — Native supplied conformer weights correctly average Sterimol descriptors.

1. **StericX claim / audited proposition:** Native supplied conformer weights correctly average Sterimol descriptors. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Normalize supplied nonnegative weights then Σwi di; Å; preserve ensemble extrema/count.
4. **Implementation location:** `src/reaction.rs::{conformer_weights,record_from_ensemble}; src/commands/parse.rs`
5. **Independent validation method:** Native frozen CLI on an analytic XYZ pair; independent packed-record decoding.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Accepted cases max descriptor error 2.6703e-7 Å; 1:3 supplied weights give L=5.2, B1=1.7, B5=4.97.
7. **Confidence:** High for specified cases
8. **Limitations:** Does not establish the supplied weights are thermodynamically justified.

## C13 — Missing weights with supplied conformer energies produce Boltzmann populations.

1. **StericX claim / audited proposition:** Missing weights with supplied conformer energies produce Boltzmann populations. Stronger Boltzmann interpretation tested by the audit; the function itself explicitly implements a uniform fallback when weights are absent.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Energy-derived populations require an explicit equation and T; uniform averaging follows equal state weights.
4. **Implementation location:** `src/reaction.rs::conformer_weights`
5. **Independent validation method:** Native parse CSV with unequal energies 0;1 and missing weights.
6. **Result:** **INCORRECT** — Native path uses a uniform average, L=4.7 Å, irrespective of energies; differs from a Boltzmann expectation.
7. **Confidence:** High
8. **Limitations:** Uniform fallback is explicit implementation behavior; documentation must not describe all native averages as energy-derived.

## C14 — Missing or invalid native conformer geometry does not silently produce a partial ensemble.

1. **StericX claim / audited proposition:** Missing or invalid native conformer geometry does not silently produce a partial ensemble. Validation behavior tested by the audit, not a quoted promise of partial-ensemble renormalization.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf; nativeingestionspecifiedsemantics
3. **Exact convention expected:** Reject an incomplete ensemble or explicitly track excluded population; avoid silent partial results.
4. **Implementation location:** `src/commands/parse.rs`
5. **Independent validation method:** A missing file and a malformed XYZ as the second conformer.
6. **Result:** **SUPPORTED** — Both return errors and no packed record; the native path does not silently discard the bad conformer.
7. **Confidence:** High for cases
8. **Limitations:** Python preparation intentionally filters some failed energies; these are distinct stage semantics.

## C15 — Buried-volume ensemble weightedmeans implement normalized sums.

1. **StericX claim / audited proposition:** Buried-volume ensemble weightedmeans implement normalized sums. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Use the same populations for each scalar; treat extrema and conformer-identity reductions separately.
4. **Implementation location:** `src/geometry/buried_volume.rs::aggregate`
5. **Independent validation method:** One/two/three/nonunit/permuted analytic descriptors versus Decimal.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Maximum weighted-mean error 6.3578e-7 Å³.
7. **Confidence:** High for finite validatedinputs
8. **Limitations:** Tiny/huge weight normalization and NaN parameters are separate failures; no complete chemical ensemble claim.

## C16 — Positive finite weights can be normalized independent of common scale.

1. **StericX claim / audited proposition:** Positive finite weights can be normalized independent of common scale. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Multiplying all weightsby any positive finite factor should not change valid probabilities when representable.
4. **Implementation location:** `src/reaction.rs::conformer_weights; src/geometry/buried_volume.rs::aggregate`
5. **Independent validation method:** Compare weights 1:3 with 1e-10:3e-10; also test sum overflow.
6. **Result:** **INCORRECT** — A valid tiny-weight ratio is rejected as zero total; a large finite sum can overflow.
7. **Confidence:** High
8. **Limitations:** Ordinary normalized weights are unaffected; a correction should rescale before summation.

## C17 — Public buried-volume aggregation has a scientifically defined result for accepted inputs.

1. **StericX claim / audited proposition:** Public buried-volume aggregation has a scientifically defined result for accepted inputs. Public-API domain/validation test; no evidence that the ordinary validated CLI supplies NaN descriptors.
2. **Primary/reference source:** Weightedmean mathematicaldomain and publicAPI
3. **Exact convention expected:** Nonfinite components require an error or an explicit missing-data meaning.
4. **Implementation location:** `src/geometry/buried_volume.rs::aggregate`
5. **Independent validation method:** Direct frozen public API: NaN descriptor with weight 1.
6. **Result:** **INCORRECT** — Returns Ok with NaN vbur_boltz; exact nonfinite output preserved.
7. **Confidence:** High
8. **Limitations:** Direct API finding; normal CLI geometry validation may prevent this input.

## C18 — Stored conformer energy_span is maximum minus minimum.

1. **StericX claim / audited proposition:** Stored conformer energy_span is maximum minus minimum. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** Energy span=max(E)−min(E); alternatively require relative energies to have minimum zero.
4. **Implementation location:** `src/reaction.rs::conformer_energy_span; PackedReactionRecord metadata`
5. **Independent validation method:** Accepted native CSV energies [1,2], independently decode the packed record.
6. **Result:** **INCORRECT** — Stored span is 2 instead of 1; code uses max and does not require minimum zero.
7. **Confidence:** High
8. **Limitations:** When relative energies are correctly referenced to the minimum, max equals max-minus-min.

## C19 — ee-to-ΔΔG conversion reproduces the standard magnitude equation.

1. **StericX claim / audited proposition:** ee-to-ΔΔG conversion reproduces the standard magnitude equation. 
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** |ΔΔG|=RTln((100+|ee|)/(100−|ee|)); |ee|<100.
4. **Implementation location:** `scripts/prepare_data.py::ee_to_ddg`
5. **Independent validation method:** Independent equation for validpositive/negative ee and two temperatures.
6. **Result:** **VERIFIED WITH NUMERICAL LIMIT** — Valid interior inputs agree within floating-point precision.
7. **Confidence:** High
8. **Limitations:** Requires kinetic selectivity assumptions; absolute value intentionally discards enantiomer identity.

## C20 — 100% and out-of-range ee can be converted to a finite exact barrier.

1. **StericX claim / audited proposition:** 100% and out-of-range ee can be converted to a finite exact barrier. Mathematical interpretation of an implemented clipping behavior; no measured detection limit was supplied.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** 100% is an infinite mathematical limit; >100% is invalid; censored measurements require bounds.
4. **Implementation location:** `scripts/prepare_data.py::ee_to_ddg`
5. **Independent validation method:** Cases ±100/±101 versus the analytic limit/domain.
6. **Result:** **INCORRECT** — All are clipped to 99.999%, giving finite 7.231911... kcal/mol at 298.15 K without censoring metadata.
7. **Confidence:** High
8. **Limitations:** If clipping represents a detection limit, report the assumption and lower bound instead of an exact observation.

## C21 — Ni-hDA stored targets use their recorded 298.15 K temperature.

1. **StericX claim / audited proposition:** Ni-hDA stored targets use their recorded 298.15 K temperature. 
2. **Primary/reference source:** Corrected primary Ni-hDA SI Table S3; primary author CSV; NIST CODATA2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
3. **Exact convention expected:** Experimental 80 °C, use 353.15 K (allow original rounded R).
4. **Implementation location:** `data/reactions_raw.csv; scripts/prepare_data.py::normalize_public_sigman; Study 011 temperature description`
5. **Independent validation method:** Independent reverse conversion of all 11 labeled source pairs.
6. **Result:** **INCORRECT** — Ten pairs imply 353.113701 K instead of 298.15 K; maximum room-temperature target difference ≈0.333685 kcal/mol.
7. **Confidence:** High
8. **Limitations:** Do not silently recompute targets; preserve published values and correct provenance separately.

## C22 — Ni-hDA ID2064 target and ee are mutually consistent.

1. **StericX claim / audited proposition:** Ni-hDA ID2064 target and ee are mutually consistent. 
2. **Primary/reference source:** Corrected primary Ni-hDA Table S3 and frozen author CSV
3. **Exact convention expected:** Use the same temperature and selectivity equation for each row.
4. **Implementation location:** `data/official/ni_hda_kraken.csv row 2064`
5. **Independent validation method:** Independently convert 3% ee at 353.15 K and infer T from the published target.
6. **Result:** **UNCERTAIN** — Source 3% ee and target ≈0.028072 correspond to ≈2% ee at 353 K; implied T=235.37 K.
7. **Confidence:** High confidence in the inconsistency; uncertain which datum is correct
8. **Limitations:** The disagreement exists upstream too; resolution requires author/experimental records.

## C23 — Population-weighted steric descriptors establish experimental conformer populations.

1. **StericX claim / audited proposition:** Population-weighted steric descriptors establish experimental conformer populations. Experimental interpretation tested for evidentiary support, not a claim that the README explicitly asserts measured conformer populations.
2. **Primary/reference source:** IUPAC Green Book (2008/2012), pp.45–46,65–67: https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf
3. **Exact convention expected:** True equilibrium populations depend on appropriate free energies, state counts, solvent and temperature.
4. **Implementation location:** `scripts/prepare_data.py; scripts/stericx_quantum.py; README ensembleclaims`
5. **Independent validation method:** Assess energy sources, degeneracy assumptions and failed-conformer handling against statistical thermodynamics.
6. **Result:** **OUT OF SCOPE** — Numerical weights are verified conditionally; no independent experimental population measurements were supplied.
7. **Confidence:** High about scope
8. **Limitations:** MMFF/xTB electronic-energy weighting and finite conformer search remain model assumptions.
