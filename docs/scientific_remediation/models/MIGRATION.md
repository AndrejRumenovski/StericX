# Explicit descriptor aggregation migration

The example below creates a separate corrected-input/model destination; it does
not overwrite historical studies. It is a reproduction recipe, not an executed
scientific-result receipt. Use the corrected binary from the final reviewed
build and a new destination suffix for each run.

```bash
python3 scripts/prepare_remediation_inputs.py \
  --output .stericx/scientific_remediation/ni_hda_demo_v1

target/debug/stericx parse \
  --csv .stericx/scientific_remediation/ni_hda_demo_v1/reactions.csv \
  --xyz-dir data \
  --output .stericx/scientific_remediation/ni_hda_demo_v1/reactions.sigpack

target/debug/stericx fit \
  --data .stericx/scientific_remediation/ni_hda_demo_v1/reactions.sigpack \
  --metadata .stericx/scientific_remediation/ni_hda_demo_v1/reactions.csv \
  --output .stericx/scientific_remediation/ni_hda_demo_v1/model.json \
  --predictions .stericx/scientific_remediation/ni_hda_demo_v1/predictions.csv \
  --portable-model .stericx/scientific_remediation/ni_hda_demo_v1/portable.json \
  --descriptor-aggregation supplied_weight_mean \
  --response-temp-k 353.15 \
  --response-sign-convention 'Published absolute ddG; no R/S assignment' \
  --optimize maximize

target/debug/stericx screen \
  .stericx/scientific_remediation/ni_hda_demo_v1/portable.json \
  --library .stericx/scientific_remediation/ni_hda_demo_v1/reactions.csv \
  --temperature 353.15 --format json
```

This screen is a method-consistency demonstration on the supplied rows, not a
new holdout or experimental validation. The input helper preserves targets,
split labels and supplied weights; it does not re-equilibrate old populations.

For unrelated packed data with unknown population provenance, omit the flag or
use `--descriptor-aggregation supplied_record_values`; screen using matching
precomputed descriptor columns. A sigpack cannot prove how its descriptor values
were obtained. For demonstrably single-geometry training records, declare
`single_geometry`; provide both `Attach_Atom_Idx` and
`Primary_Bond_Vector_Idx` in reaction-style CSVs to preserve the intended bond
axis. With no row indices, the explicit CLI donor/axis convention applies and
must match training. More than one conformer is refused in this mode.

For a declared supplied-weight mean, geometry CSVs must include
`Conformer_XYZ_Paths`, `Conformer_Boltzmann_Weights`, `Attach_Atom_Idx` and
`Primary_Bond_Vector_Idx` (a single `Ligand_XYZ_Path` can represent one weighted
conformer). Each path must contain one geometry. Supplied weights must be finite,
nonnegative, positive in total and match the conformer count. All-zero, negative,
NaN, missing and mismatched weights are refused; screening never substitutes
uniform weights for this contract.

Schema 1/2 models remain readable. If aggregation is absent they retain
`unknown`; geometry fallback is disabled. They may consume supplied descriptors,
whose provenance the caller must establish. Wrapping a fit report with the
portable API records supplied values, not an inferred conformer method. Refit
with an explicit justified method to migrate geometry-driven workflows. The old
`*_boltz` feature names remain for compatibility and do not establish
thermodynamic population provenance. JSON `descriptors[].source` and CSV
`descriptor_sources` distinguish supplied from computed values when mixed.
