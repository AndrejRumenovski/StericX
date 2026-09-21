# Kinetics, selectivity and conformer audit

The rate equation is numerically accurate within the measured `f32` range, but the `simulate` command does not have enough information to infer an absolute rate from a selectivity barrier difference. Conformer weighting also has a reproducible temperature-handling failure in the CREST population-table path. These are separate findings from descriptor accuracy.

## Frozen evidence and independent equations

`inputs_manifest.json` precedes all observations. `sut_manifest.json` hashes the original executable, the observation adapter, frozen Python observations, native CLI outputs and packed records before analysis. The adapter calls the frozen Rust library; the Python observation phase imports frozen preparation modules. Neither is an independent reference. Decimal calculations in `../scripts/kinetics_audit.py` are the independent reference and do not call those kernels.

Primary references are the [IUPAC Green Book, third edition, second printing](https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf), printed pp. 45–46 and 65–67; [NIST CODATA 2022 constants](https://physics.nist.gov/cuu/Constants/Table/allascii.txt); and [official CREST keyword documentation](https://crest-lab.github.io/crest-docs/page/documentation/keywords.html). Downloaded bytes, HTTP statuses and hashes are in `sources/manifest.json`. Green Book formula pages were rendered and visually checked. Attempts to retrieve the original [Eyring 1935 paper, DOI 10.1063/1.1749604](https://doi.org/10.1063/1.1749604) and Gold Book entry returned HTTP 403; their failed responses are retained and are not represented as read full texts. CREST 2.12's original `cregen.f90` is also frozen; line 143 takes the population temperature from `env%tboltz` and line 853 calls the Boltzmann routine.

For an elementary first-order process under transition-state assumptions,

\[
k=\kappa \frac{k_B T}{h}\exp[-\Delta G^\ddagger/(RT)].
\]

This audit uses \(\kappa=1\). Other reaction orders require standard-state/concentration factors. The molecular barrier, temperature, prefactor assumptions and rate units must be specified; the formula does not establish that a regression output is an activation free energy.

For two competing irreversible pathways with equal prefactors and a common reactant,

\[
\Delta\Delta G^\ddagger=G_S^\ddagger-G_R^\ddagger,\quad
k_R/k_S=\exp(\Delta\Delta G^\ddagger/RT),\quad
p_R=k_R/(k_R+k_S),\quad ee=|p_R-p_S|.
\]

The independent reference evaluates the two pathway weights after subtracting their minimum barrier. It does not copy StericX's major/minor logistic calculation. Equal barriers give 50:50 and zero ee; reversing the difference exchanges R and S. A dataset containing only absolute ee does not identify an absolute product configuration.

For conformer states with energies \(E_i\) and degeneracies \(g_i\),

\[
w_i=\frac{g_i e^{-(E_i-E_{\min})/(RT)}}{\sum_j g_j e^{-(E_j-E_{\min})/(RT)}},\qquad
\langle d\rangle=\sum_iw_i d_i.
\]

Using MMFF potential energies or xTB electronic energies in this expression is an approximation to thermodynamic populations. It omits state-dependent free-energy and solvent contributions unless those are supplied. Unweighted geometry averaging, population averaging, and averaging predicted selectivity through a nonlinear transformation are different operations.

Independent SI constants give `R = 0.001987204258640831739961759082 kcal mol^-1 K^-1`. The exact SI definitions of kB, h and NA and the thermochemical calorie conversion are used. The current NIST Hartree conversion gives `627.5094740628974597648183556 kcal/mol`; StericX's `627.5094740631` differs negligibly for these tests.

## Kinetic numerical results

195 requests include eight temperatures from 50 to 1000 K; positive/negative differences and barriers, zero, signed zero, very small differences, invalid inputs, and finite extreme temperatures. `results/kinetics.csv` retains every descriptor and high-precision reference string, rather than only aggregate statistics.

| Quantity | Observed numerical limit | Interpretation |
|---|---:|---|
| R or S percentage | maximum absolute error 3.7595557672e-6 percentage points | Normal `f32` rounding at tested valid inputs |
| Absolute ee | maximum absolute error 7.5191115343e-6 percentage points | Small nonzero differences can round to zero ee |
| Rate, 0–50 kcal/mol, 200–1000 K, normal nonzero `f32` results | maximum relative error 6.9231311738e-8, N=98 | Well reproduced equation, conditional on barrier meaning |
| Same range including subnormal results | maximum relative error 7.4135428360e-6, N=100 | Loss of relative precision near underflow |
| Wider finite grid | 16 underflowed rates and 15 nonfinite rates | Output range is limited; these are preserved, not omitted |

For the widest grid the maximum finite absolute rate error is `8.5554991360e26 s^-1` at an enormous rate; this is not a useful chemical accuracy measure on its own. The row-level reference and relative errors make the scale explicit. Invalid input temperatures/barriers return NaN in the library; the CLI rejects invalid input, but a finite extreme input can still produce an infinite or zero output. The observer records IEEE bits and explicit nonfinite labels so JSON null cannot hide these results.

### KIN-01: the same ΔΔG does not determine an absolute rate — INCORRECT

1. **Failing input:** native `simulate --ddg 1 --temp 298.15`, frozen in `frozen_outputs/cli_runs.json`.
2. **StericX result:** an absolute `rate_constant_s^-1` calculated by inserting 1 kcal/mol into the absolute-barrier equation, alongside selectivity calculated from that same 1 kcal/mol as a difference.
3. **Independent result:** barriers (10,11) kcal/mol give R/S rates `290543.5244386831` and `53728.67292149751 s^-1`; barriers (20,21) give `0.01358815004787473` and `0.002512784516333518 s^-1`. Both pairs have the same R:S ratio `5.407606565365828`.
4. **Definition:** rates require each absolute activation free energy; their difference supplies the ratio under shared-prefactor assumptions.
5. **Cause:** `src/commands/simulate.rs` sends one `ddg_kcal` value into both APIs.
6. **Magnitude:** the rate is unidentifiable from this input; it is not a small rounding discrepancy.
7. **Affected claims:** README's simulate output interpretation, command help, any downstream interpretation of this absolute rate.
8. **Proposed correction:** separate absolute barrier(s) from barrier difference, or omit the absolute rate when only a difference is provided. Do not alter production code until its intended interface is decided.

## Conformer results and failures

The independent checks cover one, two equal-energy, two known-energy, three, shifted-energy, extreme-separation, temperature-change, energy-window, failed and nonfinite-conformer cases. The MMFF tests inject predetermined optimizer results into the actual frozen weighting routine; they validate filtering and mathematics, not conformer search completeness or MMFF chemical accuracy. Their maximum population error against 70-digit Decimal weighting is `1.9428902931e-16`. Failed (`status<0`) and nonfinite energies are discarded; finite `status=1` results are retained even though RDKit labels status 1 as not converged. That limitation must accompany any claim of fully optimized ensembles.

Native `parse` tests use two transparent XYZ structures with analytic bond-axis dimensions `(3.7,1.7,3.47)` and `(5.7,1.7,5.47)` Å. Independent byte decoding confirms supplied weights 1:3 produce `(5.2,1.7,4.97)` to `2.10e-7 Å`; across the accepted native cases the maximum error is `2.67e-7 Å`. Missing or malformed conformers cause a recorded error, rather than silently renormalizing a partial native ensemble. Missing weights give a uniform mean even if unequal relative energies are present. This is a fallback convention, not an energy-derived Boltzmann calculation.

Direct buried-volume aggregation agrees within `6.36e-7 Å³` in the tested finite examples. Unnormalized weights 1:3 normalize correctly, as do a three-state example and a permutation. However, multiplying valid weights by a tiny positive common factor causes rejection because the `f32` sum is compared with `f32::EPSILON`; it is not scale-invariant normalization. Very large finite weights also overflow their sum. A direct public aggregate call accepts NaN descriptor input and returns a nonfinite result. These direct-API findings do not imply the CLI's validated geometry path normally produces NaNs.

### CONF-01: CREST table populations ignore configured temperature — INCORRECT

1. **Failing input:** `inputs/crest_default.log`, two states separated by 1 kcal/mol, table populations at 298.15 K; `QuantumConfig.temperature_k=500`.
2. **StericX result:** populations remain approximately `(0.8439354867,0.1560645133)`.
3. **Independent result:** `(0.7323226909,0.2676773091)` at 500 K. See full-precision rows in `results/crest_weights.csv`.
4. **Definition:** Boltzmann populations depend on the specified thermodynamic temperature; official CREST provides `--temp` (default 298.15 K).
5. **Cause:** the external CREST command does not pass `--temp`; `_conformer_thermodynamics` uses the parsed table without reading or checking its temperature. Only its fallback energy path uses the configured temperature.
6. **Magnitude:** population error `0.1116127958...` (11.16 percentage points) for this ordinary energy gap; descriptor errors depend on differences between conformers.
7. **Affected claims:** temperature-controlled quantum ensemble populations and any weighted descriptors using this path.
8. **Proposed correction:** pass and verify the requested CREST population temperature, or recompute from appropriately defined state energies and degeneracies. Preserve state-degeneracy semantics. No external CREST run is claimed here: the actual frozen parser/thermodynamics function and command construction were tested/inspected.

Further preserved cases: a malformed CREST table with negative populations is accepted by that helper; normalization only checks finiteness and positive total. The absent-table fallback reproduces energy-derived populations to about `1.1e-12` on the two-state test but assumes degeneracy one. The stored `energy_span` for accepted nonzero-offset relative energies `[1,2]` is 2, although their actual max-minus-min span is 1; either require a zero minimum or define this field as maximum relative energy.

## Experimental target temperature and ee conversion

All 11 nonempty published Ni-hDA ee/ΔΔG pairs were checked independently. Ten imply approximately `353.113701 K` under the exact modern R, consistent with 80 °C and a historically rounded gas constant. StericX stores `Temp_K=298.15` beside those targets. The maximum target discrepancy from recomputing at 298.15 K is approximately `0.333685 kcal/mol` for 91% ee. This is a temperature-provenance error, not evidence that the stored published energy values should be silently rewritten. `results/ni_hda_target_temperature.csv` preserves every row; the model audit supplies the corrected primary SI showing the 80 °C reaction.

ID 2064 is inconsistent separately: source ee=3%, ΔΔG≈0.028072 kcal/mol, implying 235.37 K; the ΔΔG is consistent with about 2% ee at 80 °C. This inconsistency is also present in the corrected SI table and must be resolved against experimental provenance, not automatically assigned to either implementation.

For valid `|ee|<100`, the preparation helper reproduces `RT log((100+|ee|)/(100-|ee|))` within floating-point error. At exactly 100% or an invalid 101%, it silently clips to 99.999%, producing a finite 7.231911... kcal/mol at 298.15 K. Exact 100% has an infinite mathematical limit; a finite measured detection limit can justify a lower bound, but none is encoded by this helper. Magnitude conversion discards sign by design and cannot establish which absolute enantiomer is favored.

## Reproduction and limits

From the repository root:

```sh
uv run --extra science python docs/scientific_accuracy_audit/scripts/kinetics_audit.py analyze
```

This verifies input/output hashes and recomputes references from frozen observations. `freeze` is deliberately non-overwriting. To regenerate observations independently, use a fresh audit directory with the same frozen source/binaries/environment, then run `freeze` and `analyze`. All negative findings above are retained in the raw JSON/JSONL, CSV and `.sigpack` artifacts. Reference sources, scripts and results are also covered by the final audit manifest.

Nothing here establishes that a force-field ensemble is complete, that a predicted ΔΔG generalizes to new chemistry, or that Eyring theory describes a particular multistep catalytic reaction. Those are methodological or experimental questions, separate from these numerical checks.
