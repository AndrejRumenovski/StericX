# StericX Study 004 — Scaled Validation

## Buried volume on Kraken's DFT geometries across the full ligand set

Study 004 reproduced the published `vbur_max_delta_qvbur_min` on eleven Ni-hDA ligands
(`R^2 = 0.9986`). This scaled run repeats the identical experiment,
using Kraken's documented 2.28 Å reference-metal distance, on every Kraken
ligand that has both a published descriptor value and DFT geometry, with exactly
one phosphorus donor. Diverse organophosphorus chemistry is included, and the
resulting `R^2` is reported as-is.

| Quantity | Value |
|---|---:|
| Eligible ligands (published value) | 1566 |
| Ligands validated (DFT + single-P) | 1541 |
| Conformers | 31611 |
| R² vs published (1:1) | 0.9852 |
| Pearson r | 0.9927 |
| RMSE | 0.4906 Å³ |
| Mean absolute error | 0.2705 Å³ |
| Median absolute error | 0.1115 Å³ |
| R² at 2.1 Å baseline (before convention fix) | 0.9553 |
| Median abs error at 2.1 Å baseline | 0.2417 Å³ |
| Ni-hDA subset R² (11 ligands, 2.28 Å) | 0.9986 |
| Study 002 R² (RDKit/MMFF geometry) | 0.8626 |
| Study 003 R² (CREST/GFN2-xTB geometry) | 0.9254 |

![Scaled parity](kraken_dft_scaled_parity.png)

Skipped ligands by reason (download): `{'ok': 1546, 'no_dft': 20}`.
Skipped during build: `{'donor_not_single_phosphorus': 3, 'donor_has_no_bonded_heavy_neighbour': 1, 'donor_not_trivalent': 1}`.

## Robust statistics for a heavy-tailed error distribution

The absolute residual is strongly right-skewed: the median (0.1115 Å³) is far below the mean (0.2705 Å³), so a handful of ligands carry a large share of
the squared error that the headline 1:1 `R²` penalises quadratically. Reporting
the residual percentiles and a trimmed `R²` is the correct way to summarise such
a distribution — it is **not** outlier removal, because the headline `R²` above
remains the full-set value over every validated ligand.

| Robust quantity | Value |
|---|---:|
| Absolute residual, 90th percentile | 0.7080 Å³ |
| Absolute residual, 95th percentile | 1.0781 Å³ |
| Absolute residual, 99th percentile | 1.8784 Å³ |
| R² excluding the worst 1 % of ligands | 0.9897 |
| R² excluding the worst 5 % of ligands | 0.9936 |

For ninety percent of the set the descriptor is reproduced to within
0.71 Å³; the full-set `R²` of 0.9852
rises to 0.9897 once the thin 1 % tail is set
aside.

## Interpretation

Across 1541 chemically diverse ligands, running the StericX
buried-volume kernel on Kraken's own DFT geometries reproduces the published
descriptor with `R^2 = 0.9852` (median absolute error
0.1115 Å³). This generalizes the eleven-ligand Study 004
result well beyond the Ni-hDA chemotype and confirms that the kernel — not the
geometry — was never the limiting factor. Adopting Kraken's documented 2.28 Å
reference-metal distance (from the 2.1 Å distance used to isolate geometry) more
than halves the typical error on this set (median absolute error
0.2417 → 0.1115 Å³), resolving the
systematic offset.

Scaling to the full set also exposed a genuine kernel bug that the eleven
trisubstituted Ni-hDA ligands could never trigger. The quadrant frame had
identified a donor's three substituents as its three nearest *heavy* atoms.
That rule is exact for a trisubstituted phosphine — no non-bonded atom can lie
closer than a real P-X bond — but it silently mis-framed **primary and
secondary phosphines** (R-PH₂, R₂P-H), whose bonded hydrogens it discarded in
favour of distant non-bonded carbons. The result was either a spurious
`max_delta_qvbur = 0` (the coordination sphere caught only the donor atom) or a
gross overestimate. The frame now selects the donor's covalently bonded atoms —
hydrogens included — by covalent-radius bond detection, which is identical to
the old behaviour for trisubstituted donors and correct for the rest. This
removed every spurious zero, tightened the residual tail, and raised the
full-set `R²` from 0.9649 to 0.9852 **without discarding a single
ligand** — the validated count in fact rose as correctly three-coordinate
phosphines that the old heavy-atom filter had rejected were admitted. Data are
from the public MolSSI Kraken descriptor-library REST API (`https://descriptor-libraries.molssi.org/api/kraken`).
