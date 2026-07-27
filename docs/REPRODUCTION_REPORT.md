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

StericX is a from-scratch Rust engine for physical-organic molecular
featurization. It computes Sterimol (\(L, B_1, B_5\)) and coordination-aware
buried-volume descriptors from Cartesian coordinates. This report evaluates
whether an independent implementation can reproduce two published results: (i)
the Kraken buried-volume descriptor `vbur_max_delta_qvbur_min`, and (ii) the
ten-ligand nickel-catalyzed homo-Diels–Alder (Ni-hDA) enantioselectivity
relationship. Against `morfeus`, StericX reproduces Sterimol parameters to
\(R^2 \ge 0.9999\) and buried-volume geometry to numerical precision
(\(R^2 = 1.000000\)). A two-step controlled experiment then isolates the source
of a residual descriptor gap. Changing only the conformer geometry source, the
native descriptor rises from \(R^2 = 0.8626\) (RDKit/MMFF) to \(0.9254\)
(CREST/GFN2-xTB) to \(0.9937\) on Kraken's own DFT geometries, localizing the
gap to geometry rather than the kernel. Adopting Kraken's documented 2.28 Å
reference-metal distance (from the 2.1 Å used to isolate geometry) then resolves
the remaining offset, reaching \(R^2 = 0.9986\) (Pearson \(r = 0.9998\)). The
result generalizes: across all 1,535 Kraken ligands with a published value and
DFT geometry (31,605 conformers), the unchanged kernel reproduces the descriptor
with \(R^2 = 0.9649\) and a median absolute error of 0.11 Å³. The compact StericX
descriptors do **not** replace the published coordination-aware descriptor for
the small Ni-hDA family; that negative result is reported rather than hidden.

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

**Descriptor kernels.** Sterimol \(L, B_1, B_5\) are computed by aligning the
attachment vector to the \(z\)-axis and scanning the van der Waals envelope.
Buried volume uses a deterministic voxel grid with a virtual metal centre placed
2.1 Å from the donor; quadrant/octant occupancies yield `qvbur` and the anisotropy
descriptor `max_delta_qvbur`. All studies use identical geometric settings
(sphere radius 3.5 Å, grid density 0.01 Å³, centre distance 2.1 Å, radii scale
1.17) so that only the input geometry varies between them.

**Reference tools and data.** Sterimol and buried-volume references are computed
with `morfeus`. Published descriptors and the Ni-hDA dataset are the public
Kraken table and the Ni-Catalyzed-hDA repository. Kraken's DFT geometries were
retrieved from the public MolSSI descriptor-library REST API; the API's
per-conformer `vbur_max_delta_qvbur` minimum matches the published
`vbur_max_delta_qvbur_min`, confirming provenance.

**Quantum geometries.** CREST 2.12 / GFN2-xTB 6.4.0 (checksum-pinned) generate
conformer ensembles for the xTB-geometry study. Kraken's reference geometries
were optimized at PBE/6-31+G(d,p) with GD3BJ dispersion (Gaussian), with
PBE0/def2-TZVP single points — a level StericX does not recompute and instead
consumes directly.

## 3. Results

### 3.1 Sterimol fidelity against morfeus (11 structures)

| Parameter | \(R^2\) | RMSE |
|---|---:|---:|
| \(L\) | 1.000000 | 0.000000 Å |
| \(B_1\) | 0.999959 | 0.0105 Å |
| \(B_5\) | 1.000000 | 0.000001 Å |

The \(B_1\) residual is a known angular-discretization difference (1° scan vs a
denser search).

### 3.2 Buried-volume geometry kernel

On identical structures the StericX voxel kernel matches `morfeus` buried volumes
to \(R^2 = 1.000000\) (worst mean relative error \(8\times10^{-6}\)%). The kernel
is therefore not the source of any disagreement with published values.

### 3.3 Localizing and resolving the descriptor gap

**Step 1 — isolate the geometry.** Holding the descriptor kernel and a fixed
2.1 Å geometric coordination centre, varying only the conformer geometry source
(11 Ni-hDA ligands):

| Geometry source (2.1 Å centre) | \(R^2\) vs published Kraken | Notes |
|---|---:|---|
| RDKit / MMFF94 | 0.8626 | inexpensive force-field conformers |
| CREST / GFN2-xTB | 0.9254 | semi-empirical ensemble (322 conformers) |
| Kraken's own DFT | 0.9937 | \(r = 0.9993\), RMSE 0.5682 Å³ (135 conformers) |

Agreement rises monotonically with geometry quality and reaches \(R^2 = 0.9937\)
on the reference DFT structures, with a near-constant offset — localizing the
earlier shortfall to conformer geometry generation, not the kernel.

**Step 2 — resolve the residual.** Kraken's descriptor code
(`PL_dft_library_201027.py`) places the reference metal 2.28 Å from phosphorus,
not the 2.1 Å used above. Adopting that documented value closes the offset:

| Reference-metal distance | \(R^2\) | RMSE (Å³) | Slope |
|---|---:|---:|---:|
| 2.1 Å (geometry-isolating baseline) | 0.9937 | 0.5682 | 0.93 |
| 2.28 Å (Kraken's convention) | **0.9986** | **0.2725** | 0.98 |

At Kraken's convention the kernel reproduces the published descriptor on
identical DFT geometries to \(R^2 = 0.9986\) (Pearson \(r = 0.9998\); Fig. 1),
confirming the residual was a coordination-centre convention difference.

**Generalization to the full library.** Repeating the experiment at Kraken's
2.28 Å convention across every Kraken ligand with a published value and DFT
geometry — 1,535 ligands, 31,605 conformers spanning the full organophosphorus
chemical space — gives \(R^2 = 0.9649\), Pearson \(r = 0.9823\), and a median
absolute error of 0.11 Å³ (Fig. 2). The wider spread than the eleven Ni-hDA
ligands is expected across such diverse chemistry; the large-sample agreement
confirms the conclusion is not specific to the Ni-hDA chemotype.

![Figure 1. Buried-volume descriptor on Kraken's DFT geometries, 11 Ni-hDA ligands, at Kraken's 2.28 Å convention.](study_004/kraken_dft_parity.png)

*Figure 1. Reproduced vs published `vbur_max_delta_qvbur_min` on the eleven
Ni-hDA ligands, Kraken DFT geometries, 2.28 Å convention (\(R^2 = 0.9986\)).*

![Figure 2. The same experiment across all 1,535 Kraken ligands.](study_004/kraken_dft_scaled_parity.png)

*Figure 2. The same kernel across 1,535 ligands / 31,605 DFT conformers
(\(R^2 = 0.9649\), median absolute error 0.11 Å³).*

### 3.4 Ni-hDA enantioselectivity reproduction

Using the preregistered Kraken descriptor `vbur_max_delta_qvbur_min`, an ordinary
least-squares model over ten training ligands reproduces the published
relationship (training \(R^2 = 0.8193\), leave-one-out \(Q^2 = 0.7521\),
LOO RMSE 0.3430 kcal/mol; historical-blind ligand 723 MAE 0.3730 kcal/mol). The
CREST-geometry buried-volume model gives fixed-feature LOO \(Q^2 = 0.5941\) and a
historical ligand-723 error of 0.1107 kcal/mol.

## 4. Limitations (reported, not hidden)

- **Compact native descriptors underperform.** StericX's own Sterimol/NBO feature
  set does not replace the published coordination-aware descriptor for this
  reaction family (native-descriptor LOO \(Q^2 \approx 0.002\)). This is an
  intentional ablation.
- **Small sample.** The Ni-hDA model has ten training ligands; leave-one-out
  metrics are correspondingly unstable, and improved descriptor fidelity in §3.3
  did not raise held-out kinetic \(Q^2\).
- **Residual tail on the full set.** After matching Kraken's 2.28 Å convention
  the median absolute error is 0.11 Å³, but a minority of the 1,535 ligands
  scatter further, where the geometrically inferred lone-pair centre diverges
  most from Kraken's exact convention.
- **No prospective validation.** A frozen ten-candidate deck exists with
  predictions recorded; its experimental outcomes are unmeasured, so no
  predictive-success claim is made.
- **DFT not recomputed.** §3.3 consumes Kraken's published DFT geometries rather
  than regenerating them (Gaussian + NBO7 are proprietary; full re-optimization
  is out of scope).

## 5. Conclusion

An independent reimplementation reproduces published Sterimol and buried-volume
descriptors to reference precision. A two-step controlled experiment first
localizes the descriptor-value gap to conformer geometry generation
(\(R^2\): 0.86 → 0.93 → 0.99 as geometry quality increases) and then resolves the
remaining offset by adopting Kraken's documented 2.28 Å coordination-centre
distance (\(R^2 = 0.9986\), Pearson \(r = 0.9998\)). The conclusion holds across
the full 1,535-ligand library (\(R^2 = 0.9649\), median error 0.11 Å³). Finally,
the compact native descriptors are shown, honestly, not to substitute for the
published coordination-aware descriptor on a small reaction family. All passed
and failed gates are retained.

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
3. Luchini, G.; Paton, R. S. et al. `morfeus`: molecular featurizer.
   https://github.com/digital-chemistry-laboratory/morfeus.

## Data and code availability

Source code, all study drivers, provenance, and per-run results are in the
StericX repository (https://github.com/AndrejRumenovski/StericX). The §3.3
geometry experiment is reproduced by `study_kraken_dft_reproduction.py`
(11 ligands) and `study_kraken_dft_scaled.py` (the full 1,535-ligand set), both
of which download Kraken's DFT geometries from the public MolSSI API. A
one-page visual summary is at `docs/results.html`, and `REPRODUCE.md` gives a
clone-to-results walkthrough.
