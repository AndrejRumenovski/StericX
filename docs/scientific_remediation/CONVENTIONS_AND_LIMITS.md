# Conventions retained during scientific remediation

These distinctions explain why not every negative audit proposition requires a
kernel change. They are not declarations that remediation testing is complete.
The original claim IDs and full evidence remain in the
[audit inventory](../scientific_accuracy_audit/CLAIMS.md).

## Sterimol and geometry

For origin `o`, unit attachment axis `u`, atomic displacement `d_i = r_i - o`
and van der Waals radius `rho_i`, the geometric support equations are

```
L_raw = max_i (d_i · u + rho_i)
B5 = max_i (||d_i - (d_i · u)u|| + rho_i)
B1 = min_{||v||=1, v·u=0} max_i (d_i · v + rho_i)
```

StericX approximates B1 by 360 azimuthal directions. That finite angular scan
has a measured error and orientation sensitivity; it is not an exact continuous
minimization. A real attachment atom is excluded in bond-axis mode; all real
atoms contribute in virtual-dummy mode. Explicit hydrogens contribute. The raw
API L and the coordination CLI's `L_raw + 0.40 Å` are distinct quantities.
See the [independent derivation and historical-source access limits](../scientific_accuracy_audit/geometry/METHODS.md)
and [official Morfeus documentation](https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html).

The atomic table is Bondi-style with extensions/fallbacks. It is not identical
to every tool's table named “Bondi”: for example, native Br is 1.85 Å and B is
1.92 Å, versus the audited Morfeus values 1.83 Å and a 2.0 Å fallback. Matching
a fallback is not evidence of physical correctness. No radius is changed to
improve a parity plot. Unknown-element fallback radii are a method limitation.

Distance-based connectivity is an inference, not chemical proof. Cordero radii
do not justify StericX's 1.3 cutoff for every molecule. Explicit connectivity
and selected axes carry information that XYZ coordinates alone may not provide.
An ambiguous axis cannot be resolved scientifically by atom numbering.

## Buried volume

The continuum quantity is the volume inside the integration sphere covered by
the union of selected atomic spheres. `%Vbur = 100 V_occ/(4πR³/3)`.
Default sphere radius is 3.5 Å, radius scaling 1.17, density 0.01 Å³, with
hydrogen spheres excluded. Donor-bound hydrogens still belong to the frame.
The geometric center uses the negative sum of the three **unit** donor-bond
vectors at distance 2.28 Å.

The native Cartesian grid retains sphere-boundary points and includes atomic
boundary hits with `distance_squared <= radius_squared`. Zero coordinate planes
belong to positive regions. Quadrants and octants estimate their occupied
fractions separately and multiply by one quarter/eighth of analytic sphere
volume. On odd grids, those independently normalized estimators need not add
exactly to the separately estimated total. This known discretization property
is retained and must be quantified in a comparison; it is not a continuum
identity. Coarse-grid and rotation effects do not disappear because R² is high.

The [SambVca documentation](https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html)
and [Morfeus documentation](https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html)
support the conventional physical setup. They do not make every finite-grid
boundary or normalization choice identical across programs. The remediation's
permutation-symmetric reduction over donor-plane orientations must be validated
and reported separately from these unchanged per-plane conventions.

Original [Kraken DFT source](https://github.com/the-matter-lab/kraken/blob/4eaad505c1343e6083032b4a3fda47e004e19734/conf_selection_and_DFT/PL_dft_library_201027.py)
uses **raw** bond vectors for its center and density 0.001 Å³. Its Sterimol
hydrogen radius is 1.09 Å and angular scan has 3,600 directions. Neither the
2.28 Å distance nor use of the same exported coordinates establishes complete
protocol equality. Eight absolute fields do not cover all published Kraken
volumes, ratios and thermodynamic reductions.

## Ensembles, kinetics and models

Supplied normalized weights define a weighted average. Missing weights define
the documented native uniform fallback, even if energies are present. Calling
either path Boltzmann requires known energy units, temperature and state
membership. Retained finite MMFF states can be unconverged; their statuses remain
part of provenance. CREST conformer groups can contain multiple rotamer energies.
Electronic/force-field energies and finite searches are approximations to
thermodynamic populations, not verified equilibrium measurements.

With `ΔΔG‡ = G‡S - G‡R` and equal prefactors, positive differences favor R.
Absolute rates require absolute activation barriers. Exact 100% ee has an
infinite barrier-difference limit; a measured detection limit would instead need
censoring metadata. Unsigned experimental ee gives a magnitude without R/S
identity. Ni-hDA response metadata is 353.15 K; historical 298.15 K population
weights must not be relabeled as if computed at the reaction temperature.

The model vocabulary and BIC selection are design choices, not mechanistic
constraints. Fixed-feature LOO, nested-alpha LOO and fixed-feature permutation
tests have conditional scope. Historical model IDs and serialized field names
do not turn those diagnostics into full-selection validation. Ranking panels
reuse ligands, and historical holdout 723 is not scaffold-disjoint.

Joint bootstrap intervals describe coefficient/model uncertainty under the
recorded resampling procedure. Marginal coefficient bands have no guaranteed
joint coverage. Nominal Student-t prediction intervals assume an adequate linear
model and residual distribution; they do not incorporate unknown chemistry,
selection uncertainty or changed descriptor conventions. Domain distances and
range checks describe support in descriptor space, not calibrated reliability.

## Evidence that code changes cannot supply

Historical Kraken per-conformer thermochemistry and retention records remain
unavailable. Some published extrema differ from both native and independent
calculations on currently exported structures. The Ni-hDA source conflict for
ligand 2064 remains unresolved. The repository has no new prospective experimental
outcomes. These remain limitations after numerical and validation bugs are fixed.
No threshold is tuned to remove them, no disagreeing molecule is dropped, and no
historical unfavorable result is rescored as part of a documentation correction.
