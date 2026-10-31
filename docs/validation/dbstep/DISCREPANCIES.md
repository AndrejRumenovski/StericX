# DBSTEP discrepancy investigation

No production numerical correction is justified by this campaign. The strict
buried-volume gate remains failed for one nondefault-radius case. Explained
finite-grid differences are not relabelled exact numerical agreement.

## Classification ledger

| Observation | Earliest divergence | Classification | Independent check |
|---|---|---|---|
| Default Sterimol offsets, up to roughly 0.11 Å | Bondi H radius: DBSTEP 1.09, native 1.20 Å | Reference convention difference | Supplied common radii and Morfeus remove this surface difference |
| Three backward witnesses, both lanes | DBSTEP `abs(z)+r` versus native signed `z+r` L | Reference convention difference | Morfeus signed L agrees with native; frozen rotated/permuted witnesses preserve the distinction |
| B1 phase differences, common-radius maximum 0.0189121 Å | Finite angular grids and frame phase | Numerical integration difference | Morfeus 3600-direction support, analytic symmetric witness, independently derived angular bound |
| Default hard-sphere comparisons, independent lattices | Different grid placement/spacing and frame rotations | Numerical integration difference | Every case retained at three DBSTEP spacings and triangulated with unmodified Morfeus |
| H-included volume cases | Different H radii, before integration | Reference convention difference | Common-radius matched-grid lane passes; distinct-default results retained |
| R=2.5 Å, ligand 401 conformer 0 | f32 versus f64 sampling-sphere inclusion | Numerical integration difference | Exact rational point distances, independent IEEE model, native point observations, DBSTEP and Morfeus |
| First volume controller's unsupported `radius` keyword | Reference call settings | Implementation error in study controller | Corrected to official `r`; all applied radii asserted; failed run/source retained |

There is no geometry mismatch in the paired calls: accepted native f32 source
coordinates, original atom indices, centers, frame transforms, radii, H policy
and source hashes are frozen for each case. This does not establish identity
with historical Kraken geometries. Unresolved historical provenance and
predictive limitations remain unchanged.

## Sphere-boundary witness

The native R=2.5 Å lattice has 5,565 points. Of these, 32 have exact squared
distance greater than 6.25 Å² by 1.08372e-7 to 2.16744e-7 Å². Their native f32
squared-distance evaluation rounds to the inclusive boundary. DBSTEP's f64
sphere test removes them, leaving 5,533 points. Decimal 80-digit arithmetic on
the exactly represented f32 coordinates establishes the direction of rounding.

For planes 15, 13 and 10, native occupied counts are respectively 1,396, 1,398
and 1,397. DBSTEP's atom-occupancy calculation on the entire captured lattice
finds exactly those counts. An independent IEEE f32 distance model agrees at
every point, and reconstructs each native absolute volume bitwise. Six occupied
boundary points in each frame are excluded by DBSTEP's subsequent sphere test,
giving 1,390, 1,392 and 1,391. Ninety-six separate native boundary-point calls
confirm all 32 point decisions for each frame. No cancellation of opposing
atom-occupancy differences is hidden by the total counts.

Morfeus's independently generated f64 projection lattice has 5,575 points and
1,399 occupied points in each supplied aligned frame, giving 16.42409612 Å³.
Native mean is 16.43008804 Å³; DBSTEP on the captured/reclipped lattice is
16.45413647 Å³. None is a continuum oracle. Changing the existing native grid
to agree with one reference would change established default arithmetic and
would need a separate scientific correction proposal, not a benchmark fit.

## Evidence

[Boundary analysis](discrepancies_v1/boundary_analysis.json) records every point,
exact distance, all counts and checks. [Morfeus calls](discrepancies_v1/morfeus_all_cases.json)
cover all 69 cases, with individual planes and no removed failures. The native
single-point requests/responses, reference logs, source/build/environment hashes
and manifest accompany them. The [Sterimol evidence](sterimol_v1/case_evidence.json)
contains complete atom maps, both alignments, both radius arrays, signed and
absolute envelopes, and fresh Morfeus results for all 61 cases. All 366
Sterimol residuals remain in their original table.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .stericx/paton/reference/venv/bin/python scripts/investigate_dbstep_discrepancies.py \
  --output docs/validation/dbstep/discrepancies_replay
```

The forensic checks pass; this means the stated explanations are reproduced.
It does not change the original predeclared arithmetic gate or claim exact
rotational invariance, continuous-volume accuracy, or historical reproduction.
