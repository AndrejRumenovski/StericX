# A 5–7 minute StericX demonstration

Lead with this: **“StericX is a Rust tool for computing and searching steric descriptors,
with a reproducible workflow for examining reaction models. Its strongest evidence is
descriptor reproduction; the small native Ni-hDA model still needs scientific work.”**

## Rehearse before the meeting

From the repository root:

```bash
cargo build --release
bash scripts/demo.sh
```

The demo runs offline using the release binary and checked-in inputs. It requires Bash
and ordinary shell utilities; it does not download data, install packages, invoke Python,
or regenerate quantum calculations or published study artifacts. Build beforehand; a
first Cargo build may require dependency downloads. Rebuild after changing Rust code.

Each run creates a fresh ignored `.stericx/demo-*` directory. It saves command arguments
in `commands.txt` and numbered text logs. No existing results are overwritten.
The script resolves the repository from its own
location, so an absolute invocation works from any directory, including paths with spaces.

For the live presentation, pause between steps:

```bash
bash scripts/demo.sh --pause
```

Use the logs from a successful rehearsal as a fallback if the terminal environment changes.
The calculations are short; spend the meeting time interpreting them.

## The walkthrough

| Time | Show | Explain |
|---|---|---|
| 0:00–0:45 | Opening statement; [scaled descriptor validation](study_004/STUDY_004_SCALED.md) | The full-set **1:1 R² = 0.9852** concerns `max_delta_qvbur_min` on matched Kraken DFT geometries: 1,541 ligands and 31,611 conformers. It is not a reaction-prediction R² or an accuracy claim for every descriptor. |
| 0:45–1:45 | Step 1: ligand 723 descriptors | The CLI identifies the P donor and returns Sterimol, buried volume, quadrant asymmetry, and pyramidalization. This input is one representative geometry. State the convention: bond-axis Sterimol; 3.5 Å sphere, 2.28 Å reference-metal distance, and Bondi radii scaled by 1.17. **B1 is an approximation.** |
| 1:45–2:30 | Step 2: compare 723 and 724 | Similar lengths can coexist with different widths and buried volumes. The σ column expresses their difference relative to the shipped library's spread. These are single geometries; the library contains per-ligand conformer aggregates. |
| 2:30–3:30 | Step 3: constrained search | Of 1,541 library ligands, **271** satisfy 30–35% buried volume and B5 ≤ 8 Å. Show the five displayed hits. This query orders matches by buried volume; its empty distance column correctly indicates there is no similarity query. These are geometric matches, not claims of reaction success. |
| 3:30–4:45 | Steps 4–5: validate and inspect | Document validation checks model consistency. Scientific inspection exposes the compact model's selected `B5_x_nbo_charge` term, ten training rows, **R² ≈ 0.3625**, **LOO Q² ≈ 0.0020**, and **group-LOO Q² ≈ 0.0013**. The native compact model has poor predictive support. Keep this separate from the published-feature reproduction in [Study 001](study_001/STUDY_001.md). |
| 4:45–6:00 | [Study 011](study_011/STUDY_011.md), or the optional screening steps below | Across all 57 predefined retrospective screens, top-1 recovery was **0.158 versus 0.333 random**; **10 fits failed**. Training used ensemble-averaged Sterimol while the screening CSV used representative conformers, producing a recorded prediction mismatch. Show this limitation as part of the result. |
| 6:00–7:00 | Discussion | Ask which reaction family, descriptor convention, and scaffold holdout would make a next benchmark chemically useful. An appropriate next experiment would freeze predictions using consistent descriptor aggregation before collecting untouched outcomes. |

The Ni-hDA source response is **`ddG_abs`**, so interpret this demonstration as selectivity
magnitude; see [scientific notes](SCIENTIFIC_NOTES.md) for the response-metadata correction.
It does **not** assign the favored R/S stereoisomer. An implied ee is also a magnitude
interpreted at the recorded temperature. No prospectively measured validation has been
demonstrated by these retrospective studies.

## Optional: fit and export a retrospective example

Rehearse the extended path with:

```bash
bash scripts/demo.sh --screen
# Live version:
bash scripts/demo.sh --screen --pause
```

This adds the predefined `S001` split from the [reaction-screening tutorial](REACTION_SCREENING.md):
eight observations train the model, and three are held out. The fit uses 1,000 bootstrap
replicates, 500 permutations, and seed `20260904`, and explicitly records the magnitude-only
response convention using `--response-sign-convention`. It then inspects the model, excludes
the eight training ligands, applies deterministic diversity selection, and exports two
illustrative candidates with a provenance sidecar. All outputs stay in the new run folder.

The expected S001 diagnostics are training **R² ≈ 0.6929**, row-wise **LOO Q² ≈ 0.4625**,
and **group-LOO Q² ≈ −192.5603**. The scaffold-validation failure rules out presenting
this example as a model ready to guide an experimental campaign. Applicability labels
such as `interpolation` or `reliable` describe support in descriptor space; they do not
overturn failed predictive validation.

Open these files if Professor Sigman wants to see what a reviewer can audit:

- `fit-report.json`: validation diagnostics and fit details.
- `model.json`: saved coefficients, transformations, reaction context, and provenance.
- `frozen-predictions.csv`: the fit's withheld-row predictions.
- `09-screen.txt`: rankings, uncertainty, applicability, exclusion, and selection details.
- `illustrative-candidate-deck.csv` and `illustrative-candidate-deck.meta.json`: candidate
  export and its metadata, including hashes.

The source library contains known experimental response columns; screening ignores them,
and the exported deck excludes experimental responses. This is a retrospective software
demonstration. Also keep the screen-versus-fit conformer-aggregation mismatch visible:
the frozen fit predictions and the screen predictions are not guaranteed to agree.
For a future prospective study, use a target-free candidate library and match descriptor
conventions and aggregation between training and screening before freezing predictions.
