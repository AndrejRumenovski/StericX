# StericX: An Independent Reproduction of Physical-Organic Ligand Descriptors and a Nickel homo-Diels–Alder Selectivity Model

**Author:** Andrej Rumenovski (Dwight D. Eisenhower High School)

**Status:** Reproduction / validation study. StericX is an independent
reimplementation and is not affiliated with, endorsed by, or produced by the
Sigman or Reisman groups or the Kraken authors. It reuses only publicly released
data and descriptor definitions, which are cited below.

**Keywords:** physical-organic descriptors · buried volume · Sterimol · Kraken ·
reproducibility · organophosphorus ligands

---

## Abstract

StericX is an independent Rust implementation of molecular geometric descriptors:
Sterimol, coordination-aware buried volume, and donor pyramidalization. This report
separates numerical agreement with reference software, reproduction of published
descriptors, and reaction-prediction performance. On eleven matched structures,
Sterimol agrees with `morfeus` at \(R^2 \ge 0.9999\), with a \(B_1\) RMSE of
0.0105 Å; the buried-volume reference check on 56 conformers from eleven ligands
gives \(R^2 = 1.000000\) to the reported precision. Changing the geometry/conformer pipeline from RDKit/MMFF to
CREST/GFN2-xTB to Kraken's published DFT structures raises agreement with
`vbur_max_delta_qvbur_min` from \(R^2 = 0.8626\) to 0.9254 to 0.9937 at a fixed
2.1 Å virtual-metal distance. Adopting Kraken's documented 2.28 Å distance gives
\(R^2 = 0.9986\), while a residual remains.

Across the 1,541 ligands with matched published values and available DFT geometry
(31,611 conformers), that descriptor gives \(R^2 = 0.9852\), median absolute
error 0.11 Å³, and 90th-percentile absolute error 0.71 Å³. Unweighted means of
per-descriptor \(R^2\) are 0.9925 for eight buried-volume comparisons, 0.9887 for
six Sterimol extrema, and 0.99998 for four pyramidalization extrema. These
summaries describe agreement with published values; they are not reaction-model
scores or proof of universal geometric accuracy. Scaling also exposed a donor-frame
bug affecting phosphines with bonded hydrogens, which was fixed without removing
previously validated ligands.

The published-descriptor Ni-hDA model is reproduced on ten training ligands, but
compact native features give fixed-feature leave-one-out \(Q^2 \approx 0.002\).
A separate retrospective ligand-ranking experiment gives top-1 recovery 0.158
versus 0.333 for exact random selection, including ten failed screens among 57
planned splits. Ni and Pd cross-coupling datasets from Newman-Stonebraker et al.
provide further retrospective descriptor and classifier comparisons. In the
recorded single-core, warm-cache benchmark, StericX runs about 14× faster than
`morfeus` for the buried-volume workload; the full paired agreement is
\(R^2 = 0.9985\). A ten-candidate forecast is frozen, but no experimental
outcomes are recorded and its absolute-selectivity target needs clearer treatment
before experimental use. The evidence supports a fast descriptor implementation
with documented limitations, not a validated method for choosing new catalysts.

---

## 1. Introduction

Data-driven catalyst design relies on quantitative molecular descriptors —
Sterimol steric parameters and buried volumes for sterics, and quantum-chemical
quantities for electronics. The Kraken platform (Gensch et al., 2022) tabulated
such descriptors for 1,558 organophosphorus ligands, and a subsequent study
(Cadge et al., 2025) used them to model Ni-hDA enantioselectivity.

Reproductions of computational chemistry results are valuable because published
descriptors depend on a specific, multi-stage pipeline (conformer search →
selection → DFT → property calculation), and small implementation choices can
shift derived quantities. StericX reimplements the geometric descriptors
independently and asks two questions: does the implementation agree with
reference tools, and where do differences from published values originate?

## 2. Methods

**Descriptor kernels and conventions.** Sterimol \(L, B_1, B_5\) are computed by
aligning an attachment vector to the \(z\)-axis and scanning the van der Waals
envelope. Buried volume uses a deterministic voxel grid around a virtual metal;
quadrant/octant occupancies yield `qvbur` and the anisotropy descriptor
`max_delta_qvbur`. The geometry-source comparison in §3.3 holds the sphere radius
at 3.5 Å, grid density at 0.01 Å³, centre distance at 2.1 Å, and Bondi radii scale
at 1.17. The subsequent convention comparison and full-library studies use
2.28 Å; the grid study explicitly varies density. These settings are therefore
not identical across every study. Coordination-axis Sterimol uses the same
virtual-metal distance and a +0.40 Å correction on \(L\).

**Reference tools and data.** `morfeus` supplies the matched-geometry reference
calculations. Published descriptors and Ni-hDA responses come from the public
Kraken table and the Ni-Catalyzed-hDA repository. Kraken DFT geometries are
retrieved through the MolSSI descriptor-library API. Matching API per-conformer
`vbur_max_delta_qvbur` minima to the published `vbur_max_delta_qvbur_min` supports
the data mapping. The matched set is 1,541 ligands; this is an availability-defined
subset of the published organophosphorus collection, not all possible phosphines.
The Ni-hDA response is `ddG_abs`, the magnitude of the enantioselectivity-derived
free-energy difference. It contains no label for which enantiomer is favored.

**Geometry sources.** CREST 2.12 / GFN2-xTB 6.4.0 generate the semi-empirical
conformer ensembles. The Kraken study consumes the authors' published DFT
geometries directly; it does not reproduce the complete conformer-search,
optimization, energy, or electronic-property pipeline. The reported reference
workflow uses PBE/6-31+G(d,p) with GD3BJ dispersion for geometry optimization and
PBE0/def2-TZVP single points.

**Metrics.** Descriptor \(R^2\) measures reproduced-versus-reference agreement;
RMSE, absolute errors, and slopes expose differences a rounded \(R^2\) can hide.
A mean across descriptor \(R^2\) values is a descriptive summary, not a pooled
accuracy or an equivalence test. Reaction fit, leave-one-out, reaction-held-out,
and scaffold-held-out results are labeled separately. Repeated ligand–reaction
rows and overlapping ranking splits are not independent observations.

## 3. Results

### 3.1 Sterimol fidelity against morfeus (11 structures)

| Parameter | \(R^2\) | RMSE |
|---|---:|---:|
| \(L\) | 1.000000 | 0.000000 Å |
| \(B_1\) | 0.999959 | 0.0105 Å |
| \(B_5\) | 1.000000 | 0.000001 Å |

The \(B_1\) residual is attributed to the angular search resolution (1° scan
versus a denser search). Values printed as zero are rounded, not proofs of exact
equality.

### 3.2 Buried-volume geometry kernel

On 56 matched conformers from eleven ligands, the StericX voxel kernel agrees
with `morfeus` at \(R^2 = 1.000000\) to the reported precision. The largest
mean relative error among the six per-conformer quantities is
\(8.09\times10^{-6}\)%; ensemble difference quantities have larger relative
errors, up to \(4.47\times10^{-5}\)%. This supports fidelity of the tested integration path
under matched inputs and conventions; it does not exclude errors in other
structures, frames, or derived descriptors.

### 3.3 Geometry and convention contributions to descriptor agreement

**Step 1 — isolate the geometry.** Holding the descriptor kernel and a fixed
2.1 Å geometric coordination centre, varying only the conformer geometry source
(11 Ni-hDA ligands):

| Geometry source (2.1 Å centre) | \(R^2\) vs published Kraken | Notes |
|---|---:|---|
| RDKit / MMFF94 | 0.8626 | inexpensive force-field conformers |
| CREST / GFN2-xTB | 0.9254 | semi-empirical ensemble (322 conformers) |
| Kraken's own DFT | 0.9937 | \(r = 0.9993\), RMSE 0.5682 Å³ (135 conformers) |

Agreement increases across these three geometry/conformer pipelines and reaches
\(R^2 = 0.9937\) on the reference DFT structures, with a remaining offset. This
shows that the input pipeline explains much of the earlier shortfall. Because
the ensembles differ in size and generation method, this is not a comparison of
electronic-structure accuracy alone.

**Step 2 — test the reference-metal distance.** Kraken's descriptor code
(`PL_dft_library_201027.py`) places the reference metal 2.28 Å from phosphorus,
rather than the 2.1 Å used above. Adopting that documented value reduces the
offset:

| Reference-metal distance | \(R^2\) | RMSE (Å³) | Slope |
|---|---:|---:|---:|
| 2.1 Å (geometry-isolating baseline) | 0.9937 | 0.5682 | 0.93 |
| 2.28 Å (Kraken's convention) | **0.9986** | **0.2725** | 0.98 |

At Kraken's distance, agreement on the published DFT geometries is
\(R^2 = 0.9986\) (Pearson \(r = 0.9998\); Fig. 1). RMSE remains 0.2725 Å³,
so matching the distance improves agreement without fully reproducing the
reference coordination centre or descriptor.

**Matched-library comparison.** At 2.28 Å, all 1,541 Kraken ligands with matched
published values and available DFT geometry (31,611 conformers) give
\(R^2 = 0.9852\), Pearson \(r = 0.9927\), and median absolute error 0.11 Å³
(Fig. 2), after the frame fix in §3.5. This broadens the structural scope beyond
the eleven Ni-hDA ligands without establishing coverage of all organophosphorus
chemical space. The absolute-error distribution has 90th/95th/99th percentiles
0.71/1.08/1.88 Å³. Removing the largest-residual 1% or 5% raises \(R^2\) to
0.9897 or 0.9936; these residual-trimmed sensitivity summaries are secondary.
The untrimmed 0.9852 value remains the headline result.

![Figure 1. Buried-volume descriptor on Kraken's DFT geometries, 11 Ni-hDA ligands, at Kraken's 2.28 Å convention.](study_004/kraken_dft_parity.png)

*Figure 1. Reproduced vs published `vbur_max_delta_qvbur_min` on the eleven
Ni-hDA ligands, Kraken DFT geometries, 2.28 Å convention (\(R^2 = 0.9986\)).*

![Figure 2. The same experiment across all 1,541 Kraken ligands.](study_004/kraken_dft_scaled_parity.png)

*Figure 2. The same kernel across 1,541 ligands / 31,611 DFT conformers
(\(R^2 = 0.9852\), median absolute error 0.11 Å³).*

**Eight buried-volume quantities.** Beyond the derived `max_delta_qvbur`, the
study compares eight buried-volume quantities against published values across
1,541 ligands, using the minimum of each quantity over its conformer ensemble
(`studies/study_004_vbur_family.py`, `study_004/STUDY_004_FAMILY.md`). The
unweighted mean \(R^2 = 0.9925\) summarizes the eight comparisons below (Fig. 3);
it does not cover all published aggregation conventions.

| Descriptor | Kraken property | \(R^2\) |
|---|---|---:|
| Buried volume | `vbur_vbur` | 0.9982 |
| Quadrant \(V_\mathrm{bur}\), min | `vbur_qvbur_min` | 0.9900 |
| Quadrant \(V_\mathrm{bur}\), max | `vbur_qvbur_max` | 0.9925 |
| Octant \(V_\mathrm{bur}\), min | `vbur_ovbur_min` | 0.9982 |
| Octant \(V_\mathrm{bur}\), max | `vbur_ovbur_max` | 0.9845 |
| Near hemisphere | `vbur_near_vbur` | 0.9975 |
| Far hemisphere | `vbur_far_vbur` | 0.9940 |
| Max Δ quadrant | `vbur_max_delta_qvbur` | 0.9852 |

The repeated `max_delta_qvbur` result (0.9852) is an internal consistency check
between study drivers sharing the descriptor implementation. Near/far hemisphere
agreement also supports the orientation convention used in these comparisons.

![Figure 3. Buried-volume descriptor family vs published Kraken values.](study_004/kraken_vbur_family_parity.png)

*Figure 3. StericX vs published Kraken value for each of the eight
buried-volume descriptors, 1,541 ligands, at each descriptor's minimum over the
conformer ensemble (mean \(R^2 = 0.9925\)).*

**A second descriptor class: Sterimol.** Sterimol \(L\), \(B_1\), \(B_5\) is the
other classical steric descriptor Kraken publishes, and reproducing it tests a
completely separate kernel. It also repeated the §3.2 lesson about conventions.
StericX's default Sterimol axis runs along a P–substituent bond, but Kraken
measures Sterimol along the **coordination axis** — a virtual metal 2.28 Å from
phosphorus on the lone pair, the *same* centre the buried volume uses, with the
historical +0.40 Å Verloop correction on \(L\). The 2.28 Å distance follows the
reference convention; agreement is not exact, with median absolute errors for
\(L\) of about 0.11 Å. With the axis matched (exposed as
`stericx descriptors --sterimol-axis coordination`), StericX reproduces Kraken's
published Sterimol across the 1,541 ligands at each conformer-ensemble extreme,
mean \(R^2 = 0.9887\) (`studies/study_004_sterimol.py`,
`study_004/STUDY_004_STERIMOL.md`, Fig. 4):

| | \(L\) | \(B_1\) | \(B_5\) |
|---|---:|---:|---:|
| min over conformers | 0.9864 | 0.9815 | 0.9927 |
| max over conformers | 0.9935 | 0.9825 | 0.9955 |

![Figure 4. Sterimol vs published Kraken values, coordination axis.](study_004/kraken_sterimol_parity.png)

*Figure 4. StericX vs published Kraken Sterimol \(L\)/\(B_1\)/\(B_5\), 1,541
ligands, coordination axis, per conformer-ensemble minimum and maximum (mean
\(R^2 = 0.9887\)).*

**A third descriptor class: pyramidalization.** Kraken also publishes two
geometric pyramidalization descriptors for the donor — `pyr_P` (Radhakrishnan's
dimensionless pyramidalization) and `pyr_alpha` (the mean out-of-plane angle) —
both defined by `morfeus`' `Pyramidalization` class. Reading that definition,
`pyr_P` reduces to the absolute scalar triple product of the donor's three unit
bond vectors, \(|\det[\hat{a}, \hat{b}, \hat{c}]|\) (with `morfeus`' \(2 - P\)
acute correction), and `pyr_alpha` to the mean signed out-of-plane angle. StericX
implements these definitions in Rust using `f32` coordinates and arithmetic.
The study narrative reports double-precision checks of the closed forms against
`morfeus`; those are not the numerical precision of the native Rust path. Run on
Kraken's DFT conformers, the native kernel agrees with published values across
1,541 ligands at
each conformer-ensemble extreme (mean \(R^2 = 0.99998\);
`studies/study_005_pyramidalization.py`, `study_005/STUDY_005.md`, Fig. 5):

| | min over conformers | max over conformers |
|---|---:|---:|
| `pyr_P` | 0.999983 | 0.999977 |
| `pyr_alpha` | 0.999979 | 0.999968 |

Pyramidalization depends on the three donor bond directions and does not use the
virtual-metal centre or integration sphere. Its small residuals (RMSE
\(\sim 2 \times 10^{-4}\) for `pyr_P`, \(\sim 0.03^\circ\) for `pyr_alpha`)
are consistent with the finite precision of the cached DFT SDF coordinates.
That explanation is plausible, but a high \(R^2\) alone does not isolate the
source of the remaining difference.

![Figure 5. Pyramidalization vs published Kraken values.](study_005/kraken_pyramidalization_parity.png)

*Figure 5. StericX (native Rust) vs published Kraken `pyr_P` and `pyr_alpha`,
1,541 ligands, per conformer-ensemble minimum and maximum (mean
\(R^2 = 0.99998\)).*

### 3.4 Ni-hDA enantioselectivity reproduction

Using the published Kraken descriptor `vbur_max_delta_qvbur_min`, an ordinary
least-squares model over ten training ligands reproduces the reported
relationship: training \(R^2 = 0.8193\), fixed-feature leave-one-out
\(Q^2 = 0.7521\), and LOO RMSE 0.3430 kcal/mol. The one historical holdout,
ligand 723, has absolute error 0.3730 kcal/mol. This known historical outcome is
not a new prospective test and cannot support a holdout \(R^2\). The
CREST-geometry buried-volume model gives fixed-feature LOO \(Q^2 = 0.5941\) and
ligand-723 absolute error 0.1107 kcal/mol.

Substituting the compact StericX feature workflow (geometric Sterimol with
externally supplied donor NBO charge) selects `B5_x_nbo_charge` and gives
fixed-feature LOO \(Q^2 \approx 0.002\), with RMSE 0.6882 kcal/mol. This result
does not support using those features to predict Ni-hDA selectivity on this
small dataset. Because the feature is held fixed during this LOO calculation,
it is not an unbiased evaluation of the full feature-selection procedure.
Separately reported nested ridge and lasso baselines are also unfavorable
(\(Q^2 = -0.130\) and \(-0.350\)). The feature-ablation result and the
published-descriptor reproduction answer different questions; neither is
prospective validation.

### 3.5 A frame-construction bug surfaced and fixed at scale

The original quadrant-frame heuristic used the donor's three nearest heavy
atoms. In primary and secondary phosphines (R–PH₂, R₂P–H), this discards bonded
hydrogens and can substitute distant non-bonded carbons. Six ligands returned
spurious zero `max_delta_qvbur` values in the initial study. Because the published
quantity takes the minimum over conformers, one invalid conformer can affect the
whole ligand's result. A nearest-heavy-atom rule also has no general guarantee
of identifying all bonded substituents in other geometries.

The fix uses covalent-radius bond detection, including hydrogens in the frame:
two atoms are treated as bonded when their separation lies within 1.3 times the
sum of their Cordero covalent radii. Hydrogens define the donor frame but are
excluded from occupied volume in these buried-volume comparisons. Frame guards
reject degenerate constructions. This remains a geometry-based bonding heuristic,
not a general electronic bonding assignment.

The recorded full-set \(R^2\) rose from 0.9649 to 0.9852 and the validated count
from 1,535 to 1,541, without discarding a previously validated ligand. Every
spurious zero in that comparison was removed, and the top five ligands' share of
squared error fell from 50% to 18%. The eleven-ligand Ni-hDA comparison remained
unchanged. Residuals after the fix are evaluated separately below.

### 3.6 Evidence for a coordination-centre contribution to residual bias

Residual error remains across the library. The 1,517 tertiary phosphines have
mean signed residual −0.010 Å³ but mean absolute residual 0.256 Å³; a small
average bias does not imply negligible individual error. The nine secondary and
fifteen primary phosphines show larger positive mean residuals, +0.781 and
+1.457 Å³, respectively, for `max_delta_qvbur`.

Study 006 compares four centre-coupled descriptors (buried volume and Sterimol
\(L, B_1, B_5\)) with two centre-free pyramidalization descriptors on the same
geometries. The mean absolute standardized residual slope is 1.54 residual
standard deviations per P–H bond for centre-coupled descriptors versus 0.04 for
pyramidalization. This pattern supports the hypothesis that the geometric
lone-pair centre contributes to the P–H-dependent bias relative to Kraken's
xTB localized-molecular-orbital centre.

The comparison does not uniquely isolate that cause: these descriptors also
differ in their mathematical definitions and sensitivity to geometry. It cannot
rule out every frame, conformer, numerical, or descriptor-specific effect. A
stronger test would substitute the reference centre on identical geometries and
recompute each affected descriptor. The present result is a mechanistic
hypothesis supported by an internal comparison, with a small P–H subgroup
(`study_006/STUDY_006.md`).

### 3.7 A second published study: Ni cross-coupling

Newman-Stonebraker et al. (*Science* **2021**, *374*, 301) classify phosphines as
active or inactive using a single-node decision tree on minimum percent buried
volume, %Vbur(min). In Study 007, StericX descriptors agree with published values
at \(R^2 = 0.9992\) and mean absolute error 0.144 percentage points across 479
ligand–reaction rows from six Ni reaction datasets. Those rows span 103 ligands;
they are not 479 independent ligand structures. Experimental yields are read
locally from the supplementary information and are not redistributed.

Fitting the paper's per-reaction yield cutoffs and class weighting approximately
recovers its classifiers. The unweighted mean accuracy/MCC round to 0.69/0.50
for both StericX and the paper, but individual thresholds and scores differ. For
Reaction I, for example, StericX accuracy/MCC are 0.765/0.588 versus the paper's
0.79/0.62. Reaction V's threshold is 50.79% versus the paper's 51.53%, whereas
the other five StericX thresholds lie near 32%. These are retrospective fitted
scores, not held-out predictive performance. Per-reaction MCC ranges from 0.36
to 0.59, with wide bootstrap intervals at n = 34–89.

Pooling all 479 rows gives a threshold of 32.77% (accuracy 0.683, MCC 0.492).
Leave-one-reaction-out evaluation gives MCC 0.41–0.54, with training thresholds
32.41–33.12%. This tests transfer to a held-out reaction within the studied
family; ligands recur across training and test reactions, so it does not test
transfer to entirely unseen ligands or scaffolds.

A separate geometry comparison uses eighteen matched free-ligand structures
from the reaction authors' supplementary information. These give
\(R^2 = 0.9735\) against published %Vbur(boltz), mean signed offset +0.12
percentage points, and sixteen values within the ranges of the Kraken conformer
ensembles. This checks descriptor agreement using a different geometry source;
it compares a representative structure to an ensemble-derived reference and
does not validate a reaction prediction. The study driver and detailed results
are in `studies/study_007_crosscoupling.py` and `study_007/STUDY_007.md`.

Study 009 evaluates six Pd reaction datasets from the **same** published paper,
including two datasets the paper assembled from other groups. It recovers the
bulky-active direction for all six, with descriptor \(R^2 = 0.9994\) across 267
ligand–reaction rows and mean fitted classifier MCC 0.64 versus the paper's
0.67. These extend the retrospective comparison; Ni and Pd here are not two
independent prospective validation studies.

### 3.8 Throughput and agreement in the recorded benchmark

Study 008 compares StericX and `morfeus-ml` 0.8.0 on 1,546 input structures
(one geometry per ligand) on an AMD Ryzen 5 5600G. Both run with one software
thread, read from a warm OS file cache, and are timed end-to-end; the fastest of
three repetitions is reported. StericX takes 1.398 s (1,106 structures/s) versus
19.237 s (80 structures/s), a 13.76× speedup. StericX also computes Sterimol and
pyramidalization in its timed pass, while this `morfeus` benchmark computes
buried volume. This is a specific workflow comparison, not a speed claim for
all Python tools or all descriptors.

Among all 1,534 paired outputs, agreement is \(R^2 = 0.998549\), mean absolute
error 0.037 percentage points, and maximum absolute difference 8.55 percentage
points. The benchmark labels sixteen differences greater than 0.5 percentage
points as frame outliers; excluding these yields \(R^2 = 0.999999\) on 1,518
pairs, with maximum absolute difference 0.422 percentage points. Frame
conventions are a plausible contributor, but the cutoff is residual-based and
does not independently establish the cause of each discrepancy. The full paired
result is the primary agreement statistic.

A larger recorded run over 31,721 conformers gives a 13.80× speedup and
\(R^2 = 0.999800\) across 31,599 paired results. Absolute timings and ratios
can change with hardware, builds, workload, and software versions; these are
archived results, not measurements of the current checkout. The native binary
avoids a Python runtime for descriptor calculations. Quantum geometry generation
still uses external programs.

### 3.9 Retrospective ligand ranking: a negative result

Study 011 evaluates all 57 predefined, scaffold-disjoint three-candidate panels
from the eleven labeled Ni-hDA ligands, training on the other eight. Top-1
recovery is **0.158 versus 0.333** for exact random selection; top-2 recovery is
0.465 versus 0.667. Ten screens fail during fitting, and the primary metrics
include their predefined penalties. Among the 47 successful rankings, mean
Spearman correlation is −0.340 and pooled RMSE is 1.219 kcal/mol.

This experiment also exposed a descriptor mismatch: `fit` uses ensemble-averaged
Sterimol values while the tested CSV screening path uses a representative
conformer. Their predictions differ by as much as 1.209761 kcal/mol. The frozen
experiment therefore evaluates that complete workflow, with model and descriptor
aggregation effects confounded. Repeated ligands and overlapping splits prevent
treating 57 panels as independent experiments. These findings do not support a
claim of improved ligand selection; they identify a concrete priority for the
next validation cycle. See `study_011/STUDY_011.md` and the locked design.

### 3.10 Grid sensitivity

Study 010 sweeps grid density for sixty ligands against a finer 0.001 Å³
numerical reference. At the default 0.01 Å³ density, mean/max absolute %Vbur
differences are 0.056/0.160 percentage points. This bounds observed grid
sensitivity for that sample; it is not an analytic error bound, proof of exact
integration, or an assessment of every derived quadrant descriptor.

## 4. Limitations and prospective status

- **Descriptor agreement and reactivity are separate.** Good agreement with
  `morfeus` or Kraken does not validate a new reaction model. The compact-feature
  Ni-hDA test and retrospective ranking study are unfavorable; the latter also
  contains a fit/screen conformer-aggregation mismatch.
- **Small and overlapping samples.** Ni-hDA has ten training ligands and one
  historical holdout. Cross-coupling reaction splits reuse ligands; the 57
  ranking panels overlap. These cannot be counted as independent experiments.
- **Geometry and conventions remain material.** The library comparison consumes
  published DFT structures, uses a geometric approximation to the coordination
  centre, and has finite grid and angular-scan errors. P–H-dependent bias supports
  a centre-related hypothesis but does not explain away all residuals. Direct
  substitution of the published centre remains an open check.
- **The ten-candidate forecast is unmeasured.** The repository retains a frozen
  deck from the published-descriptor OLS model, with nominal 95% prediction
  intervals and nine of ten candidates below the leverage warning threshold.
  These intervals depend on linear-model assumptions and have no measured
  prospective coverage. SHA-256 verifies artifact identity; a locally recorded
  date and hash alone do not establish independent public preregistration.
- **The frozen protocol needs a target-interpretation amendment.** The forecast
  fits `ddG_abs`, a nonnegative magnitude. It cannot predict which enantiomer is
  favored, and its negative OLS interval limits are not evidence for the opposite
  enantiomer. The original protocol labels the feature as %Vbur even though
  `vbur_max_delta_qvbur_min` is a volume contrast in Å³. It also claims that
  changing from ee to ΔΔG makes an otherwise uninformative interval falsifiable;
  a monotonic unit transformation cannot add information. Before adopting an
  experimental protocol, document the magnitude response, response-domain
  handling, source conditions, and scoring interpretation in a dated amendment
  while preserving the original frozen deck and rule. The original criterion
  (at least six of eight intervals covering and positive rank correlation) is
  a proposed decision rule, not proof of predictive utility or calibrated 95%
  coverage. No experimental execution or outcomes are documented here.
- **Applicability labels are diagnostics.** Range, distance, and leverage can flag
  extrapolation; labels such as `reliable` do not certify predictions. In Study
  011, the interpolation stratum has descriptive RMSE 1.144 kcal/mol and nominal
  95% interval coverage 0.648 across repeated predictions.
- **Full quantum reproduction is out of scope.** The published DFT geometry and
  electronic-property workflow is not recomputed. Geometry files, conformer
  choice, atom assignment, and descriptor settings must accompany any use of the
  reproduced descriptors.

## 5. Conclusion

StericX provides a native implementation of geometric ligand descriptors with
strong measured agreement against `morfeus` and published Kraken values. The
full matched-library result for `vbur_max_delta_qvbur_min` is
\(R^2 = 0.9852\) over 1,541 ligands, median absolute error 0.11 Å³, with
additional Sterimol, buried-volume-family, and pyramidalization comparisons.
The recorded buried-volume benchmark is about 14× faster than the tested
`morfeus` workflow on one CPU thread.

The scientific contribution is the independently implemented, auditable descriptor
workflow and its measured boundaries. Published reaction relationships can be
approximately reproduced, but the current retrospective ligand-ranking result
is below random and prospective experimental validation is absent. Matching
fit/screen descriptor aggregation, testing the coordination-centre hypothesis
directly, and designing a clearly specified prospective experiment are the next
steps supported by the evidence.

## References

1. Gensch, T.; dos Passos Gomes, G.; Friederich, P.; Peters, E.; Gaudin, T.;
   Pollice, R.; Jorner, K.; Nigam, A.; Lindner-D'Addario, M.; Sigman, M. S.;
   Aspuru-Guzik, A. A Comprehensive Discovery Platform for Organophosphorus
   Ligands for Catalysis. *J. Am. Chem. Soc.* **2022**, *144* (3), 1205–1217.
   DOI: 10.1021/jacs.1c09718.
2. Cadge, J. A.; Lozano, C.; Merriman, M. T.; Oblad, P.; Sigman, M. S.; Reisman,
   S. E. A Data Science-Guided Approach for the Development of Nickel-Catalyzed
   Homo-Diels–Alder Reactions. *J. Am. Chem. Soc.* **2025**, *147* (34),
   31175–31186. DOI: 10.1021/jacs.5c09948.
3. Jorner, K. and contributors. `morfeus`: molecular featurizer.
   [Official documentation and project attribution](https://digital-chemistry-laboratory.github.io/morfeus/#about);
   [source repository](https://github.com/digital-chemistry-laboratory/morfeus).
4. Newman-Stonebraker, S. H.; Smith, S. R.; Borowski, J. E.; Peters, E.; Gensch,
   T.; Johnson, H. C.; Sigman, M. S.; Doyle, A. G. Univariate Classification of
   Phosphine Ligation State and Reactivity in Cross-Coupling Catalysis. *Science*
   **2021**, *374* (6565), 301–308. DOI: 10.1126/science.abj4213. (§3.7.)

## Data and code availability

Source code, all study drivers, provenance, and per-run results are in the
StericX repository (https://github.com/AndrejRumenovski/StericX). The §3.3
geometry experiment is reproduced by `studies/study_004_reproduction.py`
(11 ligands) and `studies/study_004_scaled.py` (the full 1,541-ligand set), both
of which download Kraken's DFT geometries from the public MolSSI API;
`studies/study_004_vbur_family.py` reproduces the §3.3 full buried-volume family
comparison, `studies/study_004_sterimol.py` the §3.3 Sterimol comparison,
`studies/study_005_pyramidalization.py` the §3.3 pyramidalization comparison, and
`studies/study_004_frame_residual.py` the §3.5/§4 residual-by-donor-class analysis. A one-page visual summary is at `docs/results.html`, and `REPRODUCE.md`
gives a clone-to-results walkthrough.
