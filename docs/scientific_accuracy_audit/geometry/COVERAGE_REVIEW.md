# Final geometry coverage review

Reviewed 21 September 2026. “Covered” below means an experiment or explicit
evidence limitation is reported, not that the scientific claim passed. Production
scientific algorithms remain unchanged. The independent references do not call
StericX kernels; the observer is explicitly the system under investigation.

| User requirement | Authoritative evidence | Result and remaining limit |
|---|---|---|
| Freeze structures, versions and outputs before comparison | Parent initial/observer manifests; `input_freeze.json`, `supplement_input_freeze.json`, `sut_freeze.json`; focused/numerical/historical/alignment phase captures | 806 initial + 23 supplemental molecular cases and 15 raw occupancy cases; 15 focused API cases, 35 CLI calls, 11 numerical cases; 20 actual historical conformers; 4 full/minimized alignment cases + 4 exact rotations. Every captured SUT output is retained. |
| Separate claim inventory with all eight requested fields | [claims_geometry.md](../claims_geometry.md), [claims.json](claims.json) | 39 separately classified claims, including source-access limits and later falsification of near-axis exactness. |
| Sterimol definitions, L correction, radii, H, axes and units | [METHODS.md](METHODS.md), [SOURCE_ACCESS.md](sources/SOURCE_ACCESS.md), frozen official Morfeus source | Direct projection/support equations; analytic B1; raw versus +0.40 Å L kept distinct. Original Verloop full text was not inspected; no complete historical-definition verification is claimed. |
| Representative Sterimol set and per-molecule metrics | [descriptor_comparisons.csv](descriptor_comparisons.csv), [metrics.json](metrics.json), [REPORT.md](REPORT.md) | Twelve graph-derived P/N molecules, two contact variants, small substituents and synthetic edge cases. Absolute/relative errors, MAE/RMSE/max/median, slope/intercept and identity R² retained. Base maxima are explicitly restricted to that population. |
| Every public BV output and individual quadrants/octants | [bin_comparisons.csv](bin_comparisons.csv), [descriptor_comparisons.csv](descriptor_comparisons.csv) | 58,176 individual regional comparisons; all three frames, all public extrema, near/far and total/% included. Holding the center fixed validates an integration comparison, not chemical center inference. |
| BV radius, center, scaling, H, lattice and inequalities | [METHODS.md](METHODS.md), `supplement_inputs.json`, `reference_results.jsonl` | Defaults and controls are explicit. Sphere/atom boundaries and sign regions are tested. Morfeus and StericX differ at odd grids; both are nonadditive there. |
| Coarse/default/finer/very fine convergence | [convergence.csv](convergence.csv), [convergence.png](convergence.png), focused analytic lens | Four molecules at densities 0.1, 0.01, 0.001, 0.0001 Å³ plus 0.008 parity probe; independent analytic lens quantifies default error 0.268072 Å³, 0.149266 percentage points. No universal floor or monotonic convergence guarantee. |
| At least 100 random rotations plus X/Y/Z and translation | [invariance.csv](invariance.csv), frozen initial inputs | Six representatives × 103 transforms = 618. Each has 100 seeded random rotations and separate X/Y/Z rotations, with random translations. All descriptors have max/mean/SD summaries. BV invariance in this campaign is not a global guarantee. |
| Atom permutations | Same files; focused tied-axis CLI cases | 120 mapped permutations preserve explicit axes. Six actual CLI tie permutations expose a 2.57 Å L/B5 axis-selection effect. |
| Float-adjacent sphere/occupancy and region boundaries | 15 raw cases; 3 sphere-grid cases; `numerical/` atom-center witnesses | Exactly on, next representable inside/outside, signed zero and smallest subnormal signs retained. Atom centers outside the sphere can still contribute through sphere intersection. |
| Donor-bound H, primary/secondary, asymmetry, planarity, collinearity, close contacts | Base/chemical-edge/focused inputs; [connectivity.json](connectivity.json) | Standard graph neighbors match; PH3 eligibility, symmetry rejection, planar fallback ambiguity, collinear zero sentinel and coincident/close cases remain documented failures/limits. |
| Historical nearest-heavy bug, including actual case reconstruction | [historical/](historical/), [focused/reference_results.json](focused/reference_results.json) | Controlled tertiary/no-bug and H/nonbonded-heavy cases plus six actual old zero ligands and all 20 primary conformers. SDF bond graphs establish intended neighbors independently; old outputs only identify cases. |
| Donor and bond inference, P/N, thresholds and multiple donors | [connectivity.json](connectivity.json), focused CLI captures | Sole element/explicit donor, inferred neighbors and heavy-reference requirements recorded. 1.3× covalent-radius cutoff is a heuristic, not universal chemical bond evidence. Charge/lone-pair/coordination validity remains unestablished. |
| Pyramidalization equations, ordering, sign/range, degeneracy | [METHODS.md](METHODS.md), focused tetrahedral/orthogonal/acute cases, permutations | Ordinary determinant/angle agreement quantified. Tetrahedral P≈0.7698 contradicts P=1 wording; acute correction can exceed 1 and alpha can be negative. Nearly collinear finite alpha≈30° becomes a zero sentinel. Original full text remains unavailable. |
| Morfeus equivalence with minimized discrepancies | Full tables, [focused/](focused/), [alignment/](alignment/) | B1 four-atom witness; symmetric PH3 and full-cover witnesses; odd-grid/radius cases; real near-axis L/B5 cases minimized to donor plus active atom. Reference limitations are retained, including Morfeus's own B5 scan and regional nonadditivity. |
| f32/f64, overflow, underflow, NaN, infinity, signed zero, cancellation | [numerical/](numerical/), precision cases, raw boundary probes, [alignment/rotation/](alignment/rotation/) | f64 references isolate rounding from definitions. Large translations erase coordinates; finite 1e38 Å input overflows B5; short/degenerate axes produce sentinels. Ordinary-coordinate near-Z alignment produces ~0.002 Å errors and is independently falsified by an exact rotation. |
| Quantify scientific scale without claiming predictive validity | [scientific_scale.json](scientific_scale.json), report scale discussion | IQR and rounded nearest-spacing comparisons are labeled descriptive scales, not experimental resolution. Close-candidate distinctions can be affected; no predictive/experimental validity inferred. |
| Preserve negative findings and harness errors | [failures.csv](failures.csv), all focused/numerical/alignment outputs; archived initial reference attempts | 65 main comparison failure records; later cases separately retained. Original incorrect Morfeus-key adapter and inferred-center attempt remain archived and labeled harness errors, not SUT failures. |
| Reproducible commands, local links and readable plots | Report command block; all linked local paths checked; both PNG figures visually reviewed | Render-only convergence x-axis inversion bug corrected; frozen data untouched. Source/input generation refuses replacement of frozen content; reference/report regeneration remains possible. |

## Master-report corrections supplied to the lead auditor

- Do not give the restricted base L/B5 maxima as the maximum over real audited
  geometries. Full-case errors reach **0.001722320168697 Å for L** and
  **0.002769611443506 Å for B5** in the controlled alignment cases. The full
  Kraken analytic campaign may establish larger corpus-wide maxima separately.
- The **0.0261794 Å B1** minimal witness is not the largest observed error:
  full conformer 963/48394 has **0.041657831296984 Å** error; its exactly rotated
  copy retains **0.040254976377550 Å** angular-scan error.
- The random campaign's L/B5 deviations do not cover near-Z failure branches.
  Exact rotations of the targeted full conformers change SUT L/B5 by
  **0.001722335815430/0.002768516540527 Å**, while independent values are invariant.
- Odd-grid nonadditivity occurs in **both** implementations. At density 0.1,
  triphenylphosphine's SUT octant sum exceeds total by **6.5763741 Å³**;
  Morfeus's region sum exceeds total by **3.4082601 Å³**. A continuum boundary
  convention does not make either finite-grid normalization physically exact.
- Published Kraken/Morfeus density **0.001 Å³** differs from compared StericX
  default **0.01 Å³**. Bondi/Cordero full-text access limits join the already
  acknowledged Verloop/Radhakrishnan limits; official software and available
  publisher abstracts support only the specifically identified facts.

Original-source access and broad chemical validity are scientifically unresolved
claims, not missing results hidden behind successful numerical tests.
