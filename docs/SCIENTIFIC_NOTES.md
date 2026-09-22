# Scientific interpretation notes

**Review updated: 2026-09-21.** These notes qualify the historical study narratives for
the current demonstration. Original study numbers, failed results, prediction files,
and preregistration rules are retained. This review did not repeat the quantum
calculations or collect new experimental outcomes.

The subsequent [independent audit](scientific_accuracy_audit/SCIENTIFIC_ACCURACY_AUDIT.md)
found demonstrated implementation failures and convention/terminology errors.
[Scientific remediation](scientific_remediation/REMEDIATION_PLAN.md) now precedes
any further performance optimization. Until the corrected build's recheck is
complete, historical validation and benchmark numbers remain historical evidence.

## Corrections established by the independent audit

- Original Kraken DFT places a geometric virtual center using the normalized sum
  of **raw bond displacements**. StericX uses the sum of **unit bond vectors**.
  The earlier electronic-LMO attribution was wrong; tertiary phosphine centers
  need not coincide. The [source/convention matrix](scientific_accuracy_audit/kraken/REPORT.md)
  also records differences in hydrogen radius, angular scan and grid density.
- The historical Study 005 values of 4.4e−16 and 2.8e−14 describe a Python float64
  algebra experiment. Actual native-versus-Morfeus maxima on the audit's 31,721
  geometries were 3.2010947048632943e−7 for P and 4.8856123654239525e−5 degrees
  for alpha. Larger published-extremum discrepancies remain unexplained by
  rounding alone.
- A molecular formula does not establish structural identity. Study 007's guard
  cannot distinguish constitutional isomers; this does not prove that its 18
  historical mappings were wrong.
- Ni-hDA experimental targets correspond to 80 °C (353.15 K), while historical
  prepared files and some model metadata recorded 298.15 K. Published targets and
  frozen predictions are preserved; temperature-dependent interpretations need
  corrected metadata. The source conflict for ligand 2064 remains unresolved.
- A fixed descriptor vocabulary does not prove a mechanistic model. Fixed-feature
  LOO/permutation diagnostics condition on the chosen model; they do not validate
  the complete feature-selection procedure. Historical ligand 723 is a held-out
  observation, not a scaffold-disjoint holdout.
- Marginal coefficient bands are not necessarily conservative. The Student-t
  prediction interval is conditional on linear-model/residual assumptions and
  excludes model-selection and chemistry-domain uncertainty. Available overlapping
  Ni-hDA panels had 92/141 nominal 95% interval coverage (65.25%); this is evidence
  against blanket calibration claims, not a target for tuning interval widths.

Each correction is linked to its independent evidence in
[CLAIMS.md](scientific_accuracy_audit/CLAIMS.md). The original audit is preserved;
new code and recheck results belong to the remediation record.

## What the evidence supports

StericX's strongest contribution is a fast, inspectable implementation of established
steric descriptors. Separate numerical fidelity, agreement with published descriptors,
and reaction prediction when presenting the results.

| Result | Scope and interpretation |
|---|---|
| Library 1:1 R² = 0.9852; RMSE = 0.4906 Å³ | Specifically `vbur_max_delta_qvbur_min`, on 1,541 eligible Kraken ligands / 31,611 DFT conformers. Of 1,566 published entries, 20 lacked DFT structures and five failed donor checks. [Study 004](study_004/STUDY_004_SCALED.md) |
| Sterimol B1 RMSE = 0.0105 Å | On the reference-validation geometries, B1 is approximate. Rounded R² = 1.000000 for other paths does not establish mathematical identity for every kernel and input. [Reproduction commands](../REPRODUCE.md) |
| Approximately 14× speedup | Study 008's 1,546 existing structures, one CPU core, Ryzen 5 5600G, morfeus-ml 0.8.0, best of three warm-cache runs. Geometry generation is excluded. It is not a guarantee for another machine or workload. [Metrics](study_008/speed_metrics.json) |
| Benchmark paired R² = 0.998549 | All 1,534 paired structures, including 16 differences above 0.5 percentage points. R² = 0.999999 is a secondary result excluding those differences. A residual threshold alone does not diagnose each discrepancy's cause. [Metrics](study_008/speed_metrics.json) |
| Ni and Pd classifier comparisons | Two reaction families from the same Newman-Stonebraker paper, separate from Ni-hDA. Repeated ligand–reaction rows are not independent ligands. Leave-one-reaction-out validation permits recurring ligands. [Study 007](study_007/STUDY_007.md), [Study 009](study_009/STUDY_009.md) |
| Native Ni-hDA LOO Q² ≈ 0.002 | The compact native descriptor model does not match the published-feature model. [Study 001](study_001/STUDY_001.md) |
| Native ranking top-1 recovery 0.158 versus 0.333 random | All 57 predefined panels, including ten failed fits. Folds overlap. Training uses ensemble-averaged Sterimol while screening uses one representative conformer; the recorded maximum prediction difference is 1.209761 kcal/mol. This confounds model-ranking quality with descriptor aggregation. [Study 011](study_011/STUDY_011.md) |

The next ranking benchmark should use the same descriptor definitions and conformer
aggregation for both training and candidates, with its split and evaluation protocol
fixed before scoring. That would be a new experiment; it must not silently replace
Study 011's unfavorable result.

## Selectivity magnitude and the frozen forecast

The source Ni-hDA target is `ddG_abs`. This is visible in
[`normalize_public_sigman`](../scripts/prepare_data.py) and the training response in
[`preregister_prediction.py`](../scripts/preregister_prediction.py).
It represents |ΔΔG‡|, with no assignment of which enantiomer is favored.

The preserved [Study 003 preregistration](study_003/PREREGISTRATION.md) contains two
interpretation errors:

- Its “%Vbur feature” is actually `vbur_max_delta_qvbur_min`, in **Å³**.
- It interprets negative interval bounds as evidence about the opposite enantiomer.
  An unconstrained OLS interval for a nonnegative magnitude can extend below zero;
  this is a model-boundary limitation, not a stereochemical prediction. The signed-ee
  columns cannot establish R/S identity.

For a physically meaningful nonnegative magnitude, the corresponding ee magnitude
is `|ee| = tanh(|ΔΔG‡| / (2RT))`. The frozen regression predictions and interval
bounds have not been clipped, transformed, or rescored during this review.

The recorded coverage-plus-rank criterion is a proposed test, not evidence of
calibrated 95% coverage. The repository contains no measured prospective outcomes.
Before using that forecast experimentally, a dated protocol amendment should resolve
the magnitude target, interval-boundary interpretation, scoring convention, and
reaction conditions without changing the original record or claiming a new untouched
preregistration. A checksum demonstrates file integrity; it does not by itself
establish experimental blinding.

## Correction to the portable example's metadata

Earlier `stericx fit` output hardcoded “positive values favor the R product” regardless
of the input dataset. New fits default to an unspecified sign convention and accept
`--response-sign-convention` to record the dataset's definition explicitly. The demo
and tutorial use:

```text
Magnitude |ddG| from ddG_abs; larger values mean greater enantioselectivity; no R/S assignment.
```

On 2026-09-16, only `inference.response.sign_convention` was corrected in the
supplemental [Study 001 portable model](study_001/stericx_portable_model.json).
Its previous file SHA-256 was
`58337a99d1062d650538486139d90ed48ede9782d7d8977cdbe12d28896143b3`.
The numerical model, coefficients, selected features, diagnostics, and predictions
are unchanged; the file digest changes because the annotation changed. The original
legacy model and frozen prediction files remain intact. The `created` timestamp in
the portable model still identifies the original fit, not this metadata correction.

## Residual mechanism remains an interpretation

The [residual analysis](study_004/STUDY_004_RESIDUAL.md) shows tertiary phosphines
have near-zero mean signed residual, but their mean absolute residual is **0.256 Å³**
and class R² is **0.9869**. They are not error-free. The 24 primary/secondary
phosphines have a larger positive mean bias that increases with P–H count.

The [cross-descriptor control](study_006/STUDY_006.md) supports a coordination-centre
contribution without identifying a unique cause. The independent audit's direct
raw-vector versus unit-vector calculations provide stronger evidence on the
largest buried-volume discrepancies. The historical electronic-LMO explanation,
restriction to P–H ligands and claim that the kernel was never a limitation were
incorrect. Remaining geometry/ensemble provenance questions are retained.
