# Buried-volume comparison with DBSTEP 2.2.0

All 69 fixed cases completed with three captured native orientations each.
The predeclared matched-lattice gate **fails** for the R=2.5 Å variant of
ligand 401 conformer 0: all three frames differ at the spherical boundary.
The other 68 cases / 204 frames pass the unchanged arithmetic limits.
No production numerical code, constant or tolerance was changed.

Native total volume is the occupied fraction times analytic sphere volume.
DBSTEP occupied-grid volume is count times nominal cell volume; its public
`occ_vol` covers the whole molecular grid, not just the buried region.
The comparison uses common fractions and explicitly derived analytic-normalized
volumes. All three quantities and their denominators are retained.

| Comparison | Field | N | MAE | RMSE | Max absolute | Median absolute | R² | Slope | Intercept |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Matched, all cases | buried_volume | 69 | 0.000350569257 | 0.00289509207 | 0.0240484302 | 1.61137295e-06 | 0.999999965 | 1.00006365 | -0.00419828677 |
| Matched, all cases | percent_buried_volume | 69 | 0.000534023764 | 0.00442343277 | 0.0367437887 | 1.42143358e-06 | 0.999999698 | 1.00006888 | -0.00284111138 |
| Matched, passing subset | buried_volume | 68 | 2.07130213e-06 | 2.75162181e-06 | 6.80246085e-06 | 1.59654261e-06 | 1 | 0.999999939 | 2.47348259e-06 |
| Matched, passing subset | percent_buried_volume | 68 | 1.52722041e-06 | 1.93549824e-06 | 5.08230176e-06 | 1.34715988e-06 | 1 | 0.999999997 | 3.2840066e-07 |
| Independent DBSTEP grid 0.05 Å | buried_volume | 66 | 0.119306304 | 0.252202763 | 1.19544539 | 0.0517064065 | 0.999734451 | 1.00563853 | -0.273974372 |
| Independent DBSTEP grid 0.05 Å | percent_buried_volume | 66 | 0.0680567147 | 0.14188015 | 0.665639334 | 0.0287905659 | 0.999685093 | 1.00656558 | -0.184478373 |

Volume errors are Å³; percentage errors are percentage points. Passing-subset
metrics do not erase the failed case. Independent-grid metrics exclude the
three H-included cases because DBSTEP Bondi H=1.09 Å and native H=1.20 Å define
different surfaces; those cases remain in the complete residual table.
All three DBSTEP spacings (0.20, 0.10, 0.05 Å) are retained. Native integration
density is a nominal point volume, not a spacing. Different lattices have no
fitted PASS threshold and convergence need not be monotonic.

At R=2.5 Å, native single-precision sphere inclusion gives 5,565 points;
DBSTEP recomputes distance in double precision and keeps 5,533 of those points.
The mean difference is -0.02404843 Å³ / -0.03674379 percentage points.
Point-level diagnosis and Morfeus triangulation are a separate forensic step.
The strict gate remains false even when a numerical difference is explained.

The first run used the unsupported DBSTEP Python keyword `radius` for the two
nondefault radii. Its warnings, incomplete comparisons and failures are frozen
in [volume_v1](volume_v1/summary.json), with the exact original controller in
[controllers/validate_dbstep_volume_v1.py](controllers/validate_dbstep_volume_v1.py).
The corrected controller uses `r` and asserts the applied radius. It also
matches the native operation order when labelling nominal cell spacing.
Neither correction changes the independent acceptance criteria.

Reproduce in a new output directory:
```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .stericx/paton/reference/venv/bin/python scripts/validate_dbstep_volume.py \
  --output docs/validation/dbstep/volume_replay
```

Evidence: [protocol](BURIED_VOLUME_PROTOCOL.md), [summary](volume_v2/summary.json),
[per-frame residuals](volume_v2/matched_frame_residuals.csv),
[per-molecule residuals](volume_v2/matched_summary_residuals.csv),
[independent lattices](volume_v2/independent_lattice_residuals.csv),
[raw reference settings/grids](volume_v2/reference_raw.json),
[receipt](volume_v2/receipt.json), [artifact hashes](volume_v2/manifest.json).
Source coordinates, atom mapping, accepted radii, frame transforms and actual
native points are preserved in the source-bound observer capture.
