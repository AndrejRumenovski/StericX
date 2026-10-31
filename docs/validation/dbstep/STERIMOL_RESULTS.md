# Sterimol comparison with DBSTEP 2.2.0

All **61 cases** completed: every one of the 56 retained MMFF conformers and
five declared synthetic witnesses. The two reference lanes produce 366 field
comparisons. No reference call failed and no unexplained discrepancy exceeded
the [pre-execution acceptance criteria](STERIMOL_PROTOCOL.md). This is a
conditional numerical result under the recorded definitions and angular bounds;
it does not make different defaults interchangeable.

## Supplied common radii

The unmodified DBSTEP support function receives the same explicit radii, with
coordinates independently aligned by DBSTEP. Native coordination L's documented
+0.40 Å correction is removed for the raw-L comparison. L's common-definition
subset has 58 cases; all three backward witnesses remain in the complete table.

| Field | N | MAE (Å) | RMSE (Å) | Median absolute error (Å) | Maximum absolute error (Å) | R² | Slope | Intercept (Å) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Raw L, common axial definition | 58 | 5.79690e-7 | 6.71479e-7 | 5.27117e-7 | 1.09579e-6 | 0.999999999999870 | 0.999999874 | 7.01727e-7 |
| B1 | 61 | 0.00537834 | 0.00740139 | 0.00416663 | 0.01891206 | 0.999881792 | 0.999782618 | 0.00121818 |
| B5 | 61 | 1.70992e-7 | 2.07650e-7 | 1.45911e-7 | 4.26773e-7 | 0.999999999999987 | 1.000000011 | -7.21368e-8 |

L/B5 agreement is at the f32 arithmetic scale for these inputs. B1 differences
reflect independently phased finite one-degree scans; the predeclared bound is
derived from the support function, not fitted to these residuals. The observed
maximum 0.01891 Å must remain visible. Agreement within a coarse angular bound
does not imply microangstrom accuracy for B1.

## Native reference defaults

The official `Dbstep.from_coords` API uses its unchanged Bondi table. Its H
radius is 1.09 Å, versus the existing native 1.20 Å. The maximum radius-change
bound is recorded independently for each case. No radius table is patched.

| Field | N | MAE (Å) | RMSE (Å) | Median absolute error (Å) | Maximum absolute error (Å) | R² | Slope | Intercept (Å) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Raw L, common axial definition | 58 | 0.105349 | 0.107409 | 0.109999 | 0.110000 | 0.996617 | 1.007964 | 0.027714 |
| B1 | 61 | 0.086321 | 0.093185 | 0.096778 | 0.122766 | 0.980308 | 1.023908 | -0.000087 |
| B5 | 61 | 0.099075 | 0.103920 | 0.110000 | 0.110000 | 0.996754 | 1.013506 | -0.011693 |

High R² coexists with approximately 0.11 Å systematic differences. It would be
incorrect to call these defaults numerically identical. Their residuals remain
within the independently defined radius-change plus arithmetic/angular bounds.

## Definition and invariance witnesses

Adding a carbon behind the reference plane changes DBSTEP's absolute-extent L
to 6.70000005 Å, while native signed raw L remains 3.15000018 Å. Morfeus signed
raw L is 3.15000010 Å. The earliest divergence is the absolute value in DBSTEP's
axial envelope, not atom mapping or an integration error. The three backward,
rotated-backward and permuted-backward geometries account for six labelled
convention differences across the two reference lanes. All remain in the raw
residuals and all-case metrics.

The equilateral projected synthetic geometry has analytical B1 = 2.4 Å,
B5 = 3.1 Å and signed raw L = 2.8 Å. The native result is 2.40000010,
3.10000014 and 2.80000004 Å; the common-radius DBSTEP result is 2.40000004,
3.10000007 and 2.79999995 Å. The atom permutation preserves native values.
Rigid rotation preserves L/B5 within floating precision but changes native B1
from 2.39999986 to 2.40608644 Å through scan phase. DBSTEP also exhibits a finite
scan-phase change. Exact rotational invariance is therefore not claimed for
either default angular discretization.

## Evidence and reproduction

* [Residuals](sterimol_v1/residuals.csv): every field, lane, comparison value,
  Morfeus check, arithmetic/angular/radius bound and classification.
* [Case evidence](sterimol_v1/case_evidence.json): source coordinates, atom maps,
  both aligned coordinate arrays, both radii arrays, centers, axis indices,
  envelope values, settings, stdout and stderr.
* [Summary](sterimol_v1/summary.json): all-case and explicitly common-definition
  metrics, including slopes/intercepts and all six convention records.
* [Receipt](sterimol_v1/receipt.json) and [manifest](sterimol_v1/manifest.json):
  pre-execution protocol, input and executable hashes, synthetic source files,
  exact commands and all artifact hashes.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .stericx/paton/reference/venv/bin/python scripts/validate_dbstep_sterimol.py \
  --engine target/debug/examples/ensemble_report \
  --output docs/validation/dbstep/sterimol_v1
```

Choose a new output directory for replay. The native example uses the Rust
scientific engine; the reference functions and constants are unmodified.
No production calculation, default, constant or tolerance was changed to obtain
this result. Buried-volume and density comparisons are separate gates.
