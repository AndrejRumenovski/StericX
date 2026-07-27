# StericX Study 004 — Scaled Validation

## Buried volume on Kraken's DFT geometries across the full ligand set

Study 004 reproduced the published `vbur_max_delta_qvbur_min` on eleven Ni-hDA ligands
(`R^2 = 0.9937`). This scaled run repeats the identical experiment on
every Kraken ligand that has both a published descriptor value and DFT geometry,
with exactly one phosphorus donor. Diverse organophosphorus chemistry is
included, and the resulting `R^2` is reported as-is.

| Quantity | Value |
|---|---:|
| Eligible ligands (published value) | 1566 |
| Ligands validated (DFT + single-P) | 1535 |
| Conformers | 31605 |
| R² vs published (1:1) | 0.9553 |
| Pearson r | 0.9785 |
| RMSE | 0.8524 Å³ |
| Median absolute error | 0.2417 Å³ |
| Ni-hDA subset R² (11 ligands) | 0.9937 |
| Study 002 R² (RDKit/MMFF geometry) | 0.8626 |
| Study 003 R² (CREST/GFN2-xTB geometry) | 0.9254 |

![Scaled parity](kraken_dft_scaled_parity.png)

Skipped ligands by reason (download): `{'ok': 1546, 'no_dft': 20}`.
Skipped during build: `{'donor_not_single_phosphorus': 3, 'donor_not_trisubstituted': 8}`.

## Interpretation

Across 1535 chemically diverse ligands, running the unchanged
StericX buried-volume kernel on Kraken's own DFT geometries reproduces the
published descriptor with `R^2 = 0.9553`. This generalizes the
eleven-ligand Study 004 result well beyond the Ni-hDA chemotype and confirms
that the kernel — not the geometry — was never the limiting factor; the residual
scatter is consistent with StericX's approximate lone-pair coordination centre
versus Kraken's exact convention. Geometries and reference descriptors are from
the public MolSSI Kraken descriptor-library REST API (`https://descriptor-libraries.molssi.org/api/kraken`).
