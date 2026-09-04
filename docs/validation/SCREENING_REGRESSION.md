# Screening regression validation

What `tests/screening_regression.rs` does and does not guarantee about
`stericx screen`.

The suite exists so that a future change cannot silently alter a screening
result the repository has already published. Every assertion is anchored to an
artifact committed to this repository, so it runs in CI with no network access
and no third-party supporting information.

Run it alone with:

```bash
cargo test --test screening_regression
```

CI runs it as its own named step, `Guard published screening results`, on Linux,
macOS and Windows — separate from the broad `cargo test --all-targets` run, so a
regression against the published studies is identifiable from the job list and
survives any future narrowing of the general test step.

## Which studies can be regression cases, and why

| Study | Model form | Usable as a `screen` regression case? |
| --- | --- | --- |
| **001** — Ni-catalysed homo-Diels-Alder | Fitted linear model over the physical-organic feature space | **Yes**, end to end |
| **007** — Ni cross-coupling | `%Vbur(min)` single-node threshold classifier | **No** — see below |
| **009** — Pd cross-coupling | `%Vbur(min)` single-node threshold classifier | **No** — see below |

Study 001 is the only published study whose fitted model can be driven through
`screen`. Its model (in both published schema versions), its library, its frozen
prediction and its evaluation are all committed, so the full pathway is
exercised against real published values.

**Studies 007 and 009 cannot be run through `screen` for a structural reason,
not a missing-data one.** They are threshold classifiers over `%Vbur(min)`, and
the screening feature space (`stericx.physical_organic.v1`) contains only an
intercept, the three Sterimol terms, donor electronics, two interaction terms
and IR frequency. It has **no buried-volume descriptor**, so these classifiers
cannot be expressed as `screen` models at all.

That is asserted rather than assumed. `the_screening_feature_space_still_cannot_express_a_buried_volume_classifier`
fails the moment a buried-volume term is added to the feature space, which is
the point at which Studies 007 and 009 should be promoted into the end-to-end
section rather than left where they are.

Separately, their per-ligand experimental yields live in copyrighted AAAS
supporting information (*Science* **2021**, *374*, 301) that this repository
deliberately does not redistribute. Even with a compatible feature space, a
per-ligand reproduction could not run in CI from committed data alone. Nothing
copyrighted was added to make this suite work.

## What is tested

### Study 001 — end-to-end, against published values

| # | Area | What is pinned |
| --- | --- | --- |
| 1 | **Model loading** | Both published artifacts load and report their schema version (1 and 2); the portable document validates with zero errors and warnings; the model still selects `B5_x_nbo_charge` over 10 training ligands in 9 groups |
| 2 | **Descriptor mapping** | The model's stored feature space drives which library columns are required (`sterimol_b5`, `nbo_charge` — both factors of the product term, and nothing else); each hit reports the descriptors it consumed; a library lacking a required descriptor is refused by name rather than defaulted |
| 3 | **Feature scaling** | The standardized coordinate is recomputed independently from the model's own recorded mean and scale and must reproduce the reported nearest-training distance; the recorded scale is a real training standard deviation, not 1 |
| 4 | **Prediction** | Every one of the 11 published ligands agrees with `RegressXPredictor` — the kernel `fit`, `evaluate` and the study driver all use — within the documented tolerance |
| 5 | **Ranking** | Order follows the direction the model records (`maximize` → descending); the top-ranked ligand is pinned; the schema-1 artifact still refuses to rank without an explicit direction, and ranks identically to the schema-2 document when given one |
| 6 | **Applicability** | Every verdict agrees with the model's own published training range; the library is entirely `interpolation`, matching the frozen artifact's `inside_training_range`; the neighbour threshold is recomputed from the published training points |
| 7 | **Uncertainty serialization** | The portable model carries a well-formed bootstrap ensemble; screening consumes it without refitting; every interval brackets its point estimate; the method name still disclaims what it does not cover; the schema-1 artifact reports no interval rather than inventing one |
| 8 | **Candidate filtering** | `--top` selects after inference so nothing is un-screened; `--exclude-tested` removes before inference and reconciles excluded + remaining against the library size; an unresolved identifier is reported |
| 9 | **Determinism** | Text, CSV, JSON and diverse-selection runs are byte-identical across repeats; an exported deck and its metadata sidecar are reproducible |

### Cross-pathway consistency

The generic `screen` pathway is checked against the study-specific one:

- Every published ligand agrees with the `RegressXPredictor` engine kernel.
- The single frozen, pre-registered blind prediction (`SIG-NIHDA-723`,
  1.2330213 kcal/mol) is reproduced, and its applicability verdict matches.
- The published evaluation MAE still follows arithmetically from the published
  prediction and experimental value, and `evaluate`'s prediction agrees with the
  frozen CSV.

**Tolerance: 1e-4 kcal/mol**, the same figure `stericx evaluate` already applies
to a frozen prediction. The gap is the f32 engine kernel against the f64
accumulation the screen loop uses; in practice the two agree to about 1e-6. The
tolerance is a documented contract, not a fudge factor — see
`tests/screening_regression.rs`.

### Studies 007 and 009 — at the level the committed artifacts support

- Their published aggregate metrics — StericX's own outputs, committed — are
  pinned exactly: ligand counts, descriptor-fidelity `n` and R², pooled
  transferability, and the presence of a threshold and paper comparison for
  every reaction. A descriptor-kernel change that moved them cannot land
  without the study being re-run and the numbers deliberately updated.
- Every published classifier threshold must remain inside the distribution of
  `percent_buried_volume` in the committed 1,541-ligand descriptor database, and
  that database must still hold 1,541 rows.

## What is **not** tested

Stated plainly, because the gaps matter more than the coverage:

- **Studies 007 and 009 are not reproduced.** Their classifiers are not re-fit,
  their per-reaction accuracy and MCC are not recomputed, and no per-ligand
  prediction is checked. The suite pins their *recorded outputs* and a
  distribution-level property of the descriptor they threshold. It does not
  verify that re-running the studies would produce those numbers again — that
  needs the copyrighted supporting information.
- **Descriptor kernels are not recomputed from geometry.** The Kraken DFT
  geometries live in a gitignored fetch cache, so the descriptor database is
  taken as committed input. Sterimol, buried-volume and pyramidalization
  kernels are covered by their own unit tests and by Studies 002–006, not here.
- **No claim about predictive accuracy.** This suite pins *reproducibility*, not
  correctness. Study 001's model has a leave-one-out Q² of about 0.002 — it does
  not model the response — and nothing here contradicts that. A test passing
  means the number has not moved, never that the number is good.
- **The blind holdout is a single ligand.** Study 001 has exactly one frozen
  prospective prediction. It is pinned, but one point is an error estimate, not
  a validation set.
- **Statistical diagnostics are out of scope.** Bootstrap coefficient intervals,
  permutation p-values, ridge/LASSO baselines and group-LOO are pinned by
  `tests/model_training_api.rs` against a golden report, not by this suite.
- **Timing and throughput are not covered.** The ~14× speed claim is Study 008's
  benchmark and is not re-measured in CI, where runner variance would make any
  threshold either meaningless or flaky.

## Keeping it honest

The suite was mutation-tested when written. Perturbing a committed artifact
makes it fail rather than pass quietly:

| Mutation | Result |
| --- | --- |
| Frozen blind prediction shifted by 1e-3 | 2 tests fail |
| A model weight perturbed by 0.1 % | 12 tests fail |
| Study 007 descriptor R² shifted by 1e-9 | 1 test fails |

If one of these tests fails after a deliberate change, the fix is to re-run the
affected study, confirm the new number is intended, and update the pin in the
same commit that changes the behaviour — never to loosen a tolerance so the
existing pin passes.
