# StericX Study 004

## Buried volume on Kraken's own DFT geometries

The StericX voxel kernel was run, unchanged, on Kraken's public DFT-optimized
conformer geometries (PBE/6-31+G(d,p), GD3BJ), downloaded from the MolSSI
descriptor-library REST API, using Kraken's documented 2.28 Å reference-metal
distance. `vbur_max_delta_qvbur_min` is the minimum over the ensemble of each
conformer's `max_delta_qvbur`, requiring geometries only.

| Quantity | Value |
|---|---:|
| Ligands | 11 |
| Conformers | 135 |
| R² vs published (1:1) | 0.9986 |
| Pearson r | 0.9998 |
| RMSE | 0.2725 Å³ |
| Study 002 R² (RDKit/MMFF geometry) | 0.8626 |
| Study 003 R² (CREST/GFN2-xTB geometry) | 0.9254 |

![Buried volume on Kraken DFT geometry](kraken_dft_parity.png)

## Per-ligand comparison

| Kraken ID | Published | StericX on DFT | Absolute error (Å³) |
|---|---:|---:|---:|
| 401 | 2.3962 | 2.4128 | 0.0166 |
| 498 | 5.5820 | 5.4317 | 0.1503 |
| 723 | 19.6818 | 19.2789 | 0.4029 |
| 724 | 27.0649 | 26.5988 | 0.4661 |
| 785 | 6.7199 | 6.6089 | 0.1110 |
| 1057 | 16.5828 | 16.1201 | 0.4627 |
| 1058 | 15.5317 | 15.3508 | 0.1808 |
| 2062 | 6.9239 | 7.1218 | 0.1979 |
| 2063 | 6.1959 | 6.1427 | 0.0533 |
| 2064 | 6.4647 | 6.5739 | 0.1092 |
| 2067 | 9.8409 | 9.5229 | 0.3180 |

## Interpretation

Swapping the geometry source alone (holding StericX's 2.1 Å geometric-centre
convention) already raised agreement with the published descriptor from Study
002's 0.8626 to R² = 0.9937, isolating the earlier
shortfall to conformer geometry generation rather than the voxel kernel.

The remaining offset was then resolved by adopting Kraken's own documented
reference-metal distance of 2.28 Å (versus the 2.1 Å used to isolate geometry):

| Reference-metal distance | R² | RMSE (Å³) |
|---|---:|---:|
| 2.1 Å (geometry baseline) | 0.9937 | 0.5682 |
| 2.28 Å (Kraken's documented convention) | 0.9986 | 0.2725 |

At Kraken's convention the kernel reproduces the published buried-volume
descriptor on identical DFT geometries to R² = 0.9986
(Pearson r = 0.9998), confirming the residual was a
coordination-centre convention difference, not the structures or the kernel.

## Provenance

Geometries and reference descriptors were downloaded from the public MolSSI
Kraken descriptor library REST API (`https://descriptor-libraries.molssi.org/api/kraken`). The API's per-conformer
`vbur_max_delta_qvbur` minimum matches the published `vbur_max_delta_qvbur_min` value,
confirming the retrieved geometries correspond to the published dataset. StericX
is an independent reproduction; cite the original Kraken work (see the README).
