# Ensemble descriptors: independent agreement and historical behavior

The fixed 11-ligand/56-conformer study compares five reductions: the first
recorded conformer, arithmetic mean, MMFF-energy weighted mean, minimum and
maximum. The [protocol](ENSEMBLE_STUDY_PROTOCOL.md) and input hashes were recorded
before fresh Morfeus outputs were collected. All conformers and comparisons are
retained; there were no failed reference calls.

## Direct independent comparison

Unmodified Morfeus 0.8.0 evaluates the same accepted coordinates, radii, explicit
centers and three plane-atom identities. The new Rust results are compared with
independently reduced Morfeus results for **14 fields × 5 reductions × 11 ligands
= 770 pairs**, in addition to **784 per-conformer field pairs**. None exceeded
the predeclared numerical triage bounds. This establishes agreement at those
finite-algorithm bounds; it does not establish continuum accuracy or validate
the inferred center against an electronic structure calculation.

Examples of maximum absolute error among the 11 arithmetic-mean comparisons:

| Descriptor | Maximum absolute error | Unit |
| --- | ---: | --- |
| L | 8.92e-7 | Å |
| B1 | 1.55032e-2 | Å |
| B5 | 1.86e-6 | Å |
| Buried volume | 6.14e-6 | Å³ |
| Percent buried volume | 4.08e-6 | percentage points |
| Maximum adjacent-quadrant difference | 2.26e-6 | Å³ |
| Pyramidalization P | 1.50e-8 | dimensionless |
| Pyramidalization alpha | 8.15e-7 | degree |

B1 has a larger residual because the default native angular scan and Morfeus's
denser scan differ. Its threshold was derived from the angular support-function
bound before observation. Averaging does not remove that approximation, and no
angular settings or scientific constants were fitted to these data.

All MAE, RMSE, median and maximum errors, R², slope and intercept results are in
[summary.json](ensemble_study_v1/summary.json), with every
[conformer residual](ensemble_study_v1/conformer_residuals.csv) and
[ensemble residual](ensemble_study_v1/ensemble_residuals.csv). Constant-reference
R²/slopes are unavailable, not assigned perfect scores. R² is `1-SSE/SST`, not
squared correlation. The reference outputs retain all per-plane volumes and
errors for further inspection.

## Historical endpoint comparison

The fixed historical target is Kraken's published **minimum** adjacent-quadrant
volume difference. The original historical geometries and energies have not
been recovered. The table compares different estimates from the retained MMFF
workflow at the explicit current 2.28 Å center convention, with no fitted
correction or selection of a production default. Units are Å³.

| Native estimate | N | MAE | RMSE | Median absolute error | Maximum absolute error | R² |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| First conformer | 11 | 3.45639 | 4.98228 | 2.67486 | 13.29608 | 0.51936 |
| Arithmetic mean | 11 | 3.66304 | 4.42543 | 3.06670 | 9.38075 | 0.62079 |
| MMFF-energy weighted mean | 11 | 3.71739 | 5.09406 | 3.06670 | 13.45051 | 0.49755 |
| Minimum | 11 | 2.53752 | 3.01874 | 2.09306 | 5.82131 | 0.82355 |
| Maximum | 11 | 5.02252 | 6.64918 | 3.94465 | 16.61801 | 0.14395 |

The arithmetic mean reduces RMSE relative to the first conformer here, but
increases MAE. Energy weighting increases both. The minimum performs best for
this **minimum-defined** endpoint, so the result cannot demonstrate universal
superiority of an ensemble method. Geometry, center, sampling, population and
reduction differences are entangled. Slopes/intercepts and all 55 molecule-level
comparisons are retained in the summary and
[historical residual table](ensemble_study_v1/historical_residuals.csv).

This is exploratory historical workflow agreement, not identical-geometry
descriptor validation or an experimental prediction benchmark. No features,
models, corrections or hyperparameters were fitted, and no held-out prediction
claim is made. The fixed historical cohort is small and previously studied;
the existing 723/2064 provenance and prediction limitations remain. The MMFF
weights are conditional on retained states at unit multiplicity and omit
thermochemical entropy. No historical energies were invented.

## Reproduction

```bash
.stericx/paton/reference/venv/bin/python scripts/study_conformer_ensembles.py \
  --output docs/paton_recommendations/conformers/ensemble_study_v1
.stericx/paton/reference/venv/bin/python -m unittest discover \
  -s tests -p test_conformer_study.py -v
```

Use a fresh output directory. The reference environment contains Morfeus 0.8.0,
NumPy 2.5.4 and SciPy 1.18.1 under Python 3.12.13; exact package/source hashes,
protocol, input manifest and command are in the
[receipt](ensemble_study_v1/receipt.json). The orchestration calls the installed
reference packages and the already independently tested Decimal population
reference. It does not implement a second production descriptor engine.
`step5_study.log` also preserves the cclib import SyntaxWarning from that
reference environment. Closed-form metric tests independently distinguish R²
from squared correlation and retain missing/constant-target cases.
