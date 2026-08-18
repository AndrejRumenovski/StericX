# StericX Portable Model Format

A StericX model document is a single JSON file that carries a fitted model
together with everything needed to score a new structure on another machine —
without the training data, the training code, or the study driver that produced
it.

The format is implemented by [`src/model/portable.rs`](../src/model/portable.rs)
as typed Serde structures. Documents are read only through
`PortableModel::from_json`, which validates before returning; there is no path
that hands a caller an unvalidated model.

## Versions

| `schema_version` | Meaning |
|---:|---|
| 1 | The original `stericx fit` artifact. Readable, but **not portable**: it cannot state what it predicts, what data trained it, or who produced it. |
| 2 | Version 1 plus the `inference`, `provenance`, and `created` sections. |

Version 2 is a **strict superset** of version 1: every version 1 key keeps its
name, position, and meaning at the document root, and three sections are added
after them. Two consequences follow.

- A reader written against version 1 — including `stericx evaluate` and the
  Python study drivers — deserializes a version 2 document unchanged.
- A version 1 document loads as a `PortableModel` whose added sections are
  absent. `is_portable()` returns `false` and `inference()` / `provenance()`
  return `MissingSection`. Nothing is invented to fill the gap.

A document whose `schema_version` exceeds the version this build understands is
rejected outright, before any other field is examined, so a future format is
never partially interpreted.

## Document layout

```jsonc
{
  // ---- schema version 1 fields, unchanged --------------------------------
  "schema_version": 2,
  "model": "mechanistically_constrained_ols",
  "training_count": 10,
  "training_group_count": 9,
  "feature_names": ["intercept", "L_boltz", "...", "ir_frequency"],
  "selected_feature_indices": [6],
  "selected_features": ["B5_x_nbo_charge"],
  "weights": [1.863738, 0.0, 0.0, 0.0, 0.0, 0.0, -0.11433673, 0.0],
  "standardized_means":  [/* per feature column */],
  "standardized_scales": [/* per feature column */],
  "training":                {"count": 10, "r2": 0.3625, "mae": 0.4484, "rmse": 0.5500},
  "fixed_feature_loo":       {"count": 10, "r2": 0.0020, "mae": 0.5583, "rmse": 0.6882},
  "fixed_feature_group_loo": {"count": 10, "r2": 0.0013, "mae": 0.5586, "rmse": 0.6884},
  "ridge_baseline":  {"model": "ridge", "regularization": 0.01, "weights": [/*…*/],
                      "training": {/*…*/}, "nested_loo": {/*…*/}},
  "lasso_baseline":  {"model": "lasso", "...": "..."},
  "coefficient_intervals": [{"feature": "intercept", "estimate": 1.86,
                             "lower_95": -0.32, "upper_95": 2.90}],
  "response_permutation_p_value": 0.0640,
  "correlation_matrix": [[/* 8 x 8 */]],
  "variance_inflation_factors": [null, 1.0, null, "..."],
  "applicability_domain": [{"feature": "B5_x_nbo_charge",
                            "minimum": 5.0535, "maximum": 14.4759}],
  "notes": ["Descriptors were standardized using training rows only.", "..."],

  // ---- schema version 2 additions ----------------------------------------
  "inference":  { "...": "..." },
  "provenance": { "...": "..." },
  "created":    { "...": "..." }
}
```

### `inference`

Everything required to turn a structure into a predicted number.

| Field | Meaning |
|---|---|
| `response.name` | Machine-readable response identifier, e.g. `ddg_double_dagger`. |
| `response.units` | Physical units, e.g. `kcal/mol`. |
| `response.description` | What the quantity is. |
| `response.sign_convention` | Which enantiomer a positive value favors. Without this a consumer can invert the selectivity. |
| `response.temperature_k` | Temperature the response refers to, or `null`. |
| `response.optimization` | Which direction counts as better: `maximize`, `minimize`, `maximize_magnitude`, or `unspecified`. |
| `feature_space.definition` | Stable identifier for the feature construction, e.g. `stericx.physical_organic.v1`. |
| `feature_space.feature_names` | All eight column names, in model order. |
| `feature_space.transformations` | How each column is built (see below). |
| `intercept` | Raw-scale intercept. |
| `terms[]` | One entry per selected descriptor. |

Each entry of `terms` carries `feature_index`, `feature_name`, the raw-scale
`coefficient`, the `training_mean` and `training_standard_deviation` used when
fitting, and the `training_minimum` / `training_maximum` that define the
applicability domain for that descriptor.

`transformations` is a tagged union, so the feature vector can be rebuilt
without the StericX source:

```jsonc
{"kind": "constant"}
{"kind": "descriptor",  "descriptor": "sterimol_b5"}
{"kind": "interaction", "factors": ["sterimol_b5", "nbo_charge"]}
```

Descriptor names refer to packed-record quantities: `sterimol_l`,
`sterimol_b1`, `sterimol_b5`, `nbo_charge`, `ir_frequency`.

Scoring a structure is therefore:

```text
prediction = intercept + Σ  coefficient_t × value(transformations[index_t])
                         t
```

The coefficients are on the **raw descriptor scale**, so no standardization is
applied at inference time. `training_mean` and `training_standard_deviation`
are recorded because they are needed to interpret the fit and to re-derive the
standardized coefficients, not because inference needs them.

`optimization` exists because a prediction is only a number: whether a larger
one is preferable is a property of the chemistry. For a signed selectivity
response, `+1.5` and `−1.5` describe equal selectivity for *opposite*
enantiomers, so neither "larger is better" nor "smaller is better" is
universally correct — `maximize_magnitude` covers the case where either product
is acceptable. The field is optional and defaults to `unspecified`; documents
written before it existed read back that way. Consumers must treat
`unspecified` as *no direction stated* rather than choosing one:
`stericx screen` refuses to rank such a model until the caller says which way
is better.

### `provenance`

| Field | Meaning |
|---|---|
| `model_id` | Stable identifier. Derived deterministically from the model and its training-data digests when not supplied. |
| `stericx_version` | Crate version that produced the fit. |
| `record_format` | Record layout the descriptors were read from, e.g. `sigpack_v1`. |
| `training.record_count` | Training rows. Must agree with `training_count`. |
| `training.group_count` | Distinct scaffold groups. Must agree with `training_group_count`. |
| `training.dataset_digests[]` | One digest per training input: `artifact`, `algorithm`, `digest`, `byte_count`. Never empty. |
| `training.fit_options` | `max_terms`, `bootstrap_samples`, `permutation_samples`, and `seed` — the configuration needed to re-derive the fit. |
| `reaction.*` | Chemistry context: `reaction_family`, `catalyst_metal`, `ligand_class`, `source_url`, `notes`. |

The digest carries its own `algorithm` because the strength of the guarantee
varies. `stericx fit` currently records `fnv1a64`, which detects accidental
substitution but is **not** a cryptographic hash and proves nothing against a
deliberate one. A pipeline that computes SHA-256 can record `sha256` instead;
consumers should read the algorithm rather than assume it.

Every `reaction` field is optional and serializes as an explicit `null` when
unknown. Nothing is ever defaulted to a plausible value.
`PortableModel::missing_provenance()` lists what is absent, and `stericx fit`
prints the same list as `portable_model_missing_provenance`, so a gap is
visible rather than silent. A consumer that requires complete context can
refuse a model whose list is non-empty.

### `created`

| Field | Meaning |
|---|---|
| `created_utc` | RFC 3339 UTC timestamp. |
| `produced_by` | Producing command, e.g. `stericx fit`. |

## Validation

`PortableModel::validate` runs on every read and every write. It rejects:

- a `schema_version` above the supported maximum, or below 1;
- a version 2 document missing `inference`, `provenance`, or `created`;
- a missing required field in any section (Serde reports the field by name);
- non-finite weights; a selected index outside the descriptor columns, or
  pointing at the intercept; a duplicated selected index; a selected name that
  disagrees with its column;
- a standardization scale that is zero, negative, or non-finite, which would
  divide by zero during inference;
- an inverted or non-finite applicability range;
- an empty `model_id`, `stericx_version`, or `created_utc`;
- an empty `dataset_digests` list, or a digest that is not hexadecimal;
- provenance training counts that disagree with the fit report.

### Cross-checking

`inference` restates values that also appear in the flattened fit report:
coefficients against `weights`, standardization against `standardized_means`
and `standardized_scales`, and ranges against `applicability_domain`.
Validation requires the two copies to agree within a relative tolerance of
1e-12. That converts the redundancy into an integrity check: a document whose
`inference` block was edited without editing `weights` is rejected as
malformed, rather than scored with whichever copy a reader happened to consult.

## Determinism

Serialization is deterministic. Field order follows the Rust declarations, no
unordered maps appear anywhere in the format, and re-serializing a document
reproduces it byte for byte.

Exact float round-tripping requires `serde_json`'s `float_roundtrip` feature,
which is enabled in `Cargo.toml`. The default parser is not correctly rounded
and shifts some seventeen-digit `f64` values by one unit in the last place,
which would rewrite published validation statistics simply by reading a model
and writing it again. `tests/portable_model_format.rs` guards this.

## Producing a document

`stericx fit` writes its usual version 1 artifact to `--output` and, when asked,
an additional version 2 document:

```bash
./target/release/stericx fit \
  --data data/reactions.sigpack \
  --metadata data/reactions_raw.csv \
  --output docs/study_001/stericx_model.json \
  --predictions docs/study_001/stericx_frozen_predictions.csv \
  --portable-model docs/study_001/stericx_portable_model.json \
  --reaction-family "Ni-catalyzed homo-Diels-Alder" \
  --catalyst-metal Ni \
  --source-url "https://github.com/SigmanGroup/Ni-Catalyzed-hDA" \
  --response-temp-k 298.15
```

`--output` is unchanged by this flag, so existing study artifacts stay
byte-identical.

## Inspecting and validating a document

```bash
stericx model inspect model.json     # scientific summary
stericx model validate model.json    # every problem, with field paths
```

Both read the document in stages — JSON syntax, then `schema_version`, then the
fixed-width numeric arrays, then the typed fields, then the semantic checks —
so a reader is told which field is wrong instead of receiving a decoder message
about the whole file. `validate` reports all findings rather than stopping at
the first, exits non-zero on any error, and accepts `--strict` to fail on
warnings too. Both support `--format json`.

Findings carry a stable `code`, a `location` field path, and a `severity`:

| Severity | Meaning |
|---|---|
| `error` | The model cannot be trusted for inference. |
| `warning` | Usable, but records less than it should — a legacy document, or unrecorded chemistry context. |

Issue codes include `unsupported_schema_version`, `missing_schema_version`,
`invalid_json`, `missing_field`, `dimension_mismatch`, `non_numeric_value`,
`non_finite_value`, `invalid_scale`, `inverted_range`, `no_descriptors`,
`descriptor_name_mismatch`, `duplicate_feature_index`, `invalid_feature_index`,
`missing_section`, `missing_dataset_digest`, `malformed_digest`,
`incomplete_digest`, `training_count_mismatch`, `training_group_mismatch`,
`term_count_mismatch`, `inference_disagrees_with_fit`, `legacy_schema`, and
`unrecorded_context`.

## Compatibility policy

- Adding an optional field to an existing section does not change the version.
- Adding a required field, renaming one, or changing the meaning of one
  requires a new `schema_version`, and this build must then reject it until it
  is taught to read it.
- Version 1 documents remain readable. They are reported as legacy rather than
  upgraded in place, because the metadata a version 2 document requires cannot
  be recovered from a version 1 file — and guessing it is exactly what this
  format exists to prevent.
