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
varies. `stericx fit` records **`sha256`**. Documents written before that change
record `fnv1a64`, which detects accidental substitution but is **not** a
cryptographic hash and proves nothing against a deliberate one; they still load
unchanged. Consumers should read the algorithm rather than assume it — that is
exactly what the field is for, and a screening report prints the tag alongside
the digest so an FNV-1a value is never mistaken for a cryptographic guarantee.

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

## Applicability domain

`training_geometry` records what the training set actually covers, so a screened
ligand can be placed relative to it. It is optional — models written before it
existed load without it and are reported as `unknown` rather than assumed to be
in domain.

| Field | Meaning |
|---|---|
| `feature_indices`, `means`, `scales` | The standardized frame the fit solved in. |
| `xtx_inverse` | `(X'X)⁻¹` in that frame; the source of leverage. |
| `observations`, `parameters` | `n` and `p`, counting the intercept. |
| `residual_standard_error` | `s = √(RSS/(n−p))`. |
| `warning_leverage` | `h* = 3p/n`. |
| `standardized_training_points` | One row per training observation, one column per selected descriptor. |
| `neighbor_calibration` | The training set's own nearest-neighbour spacing. |
| `training_labels` | Optional. Identifier of each training observation, positionally aligned with `standardized_training_points`, so a screened ligand's nearest training neighbour can be named. **Identifiers only** — no experimental response is recorded, because a model document must never carry a blinded target. Absent from documents written before this field existed; consumers must degrade to reporting the distance without a name. |

### The three measures

**Per-descriptor range check.** Each selected descriptor is compared against the
`training_minimum` / `training_maximum` the model records. A departure is
reported as `normalized_exceedance = overshoot / (maximum − minimum)`, so "half a
training range beyond the edge" is 0.5 regardless of units. When a training
range has zero width — every training value identical — the exceedance is
infinite rather than a fraction of a range that does not exist.

**Nearest-neighbour distance.** Euclidean distance in the standardized
descriptor space to the closest training observation. This is what
`standardized_training_points` exists for; a range check alone cannot see a hole
in the middle of the training box.

**Mahalanobis distance.** Not a second covariance estimate. The standardized
columns are centred on the training means by construction, so `Z'Z` is block
diagonal and leverage decomposes exactly as

```text
h = 1/n + z'S⁻¹z ,  S = (n−1)·Cov   ⟹   Mahalanobis = √((n−1)(h − 1/n))
```

recovered from `xtx_inverse`, which guarantees it can never disagree with the
leverage the same model reports. It is **declined**, with a stated reason, when
it would be unreliable rather than merely inconvenient:

- `n < k + 2`, where the sample covariance of `k` descriptors is not estimable;
- `xtx_inverse[0][0] ≠ 1/n` or a non-zero intercept row, meaning the stored
  matrix is not the centred one the decomposition assumes;
- a negative squared distance from numerical error.

### `uncertainty` — the bootstrap ensemble

Schema 2 optionally carries the bootstrap replicates behind `coefficient_intervals`,
so a saved model can produce an uncertainty estimate without the training data
and without refitting.

| Field | Meaning |
| --- | --- |
| `method` | How the replicates were generated. |
| `replicate_count` | Replicates actually usable; a resample with a singular design is skipped, so this can be below `requested_samples`. |
| `requested_samples` | What the fit was asked for. |
| `seed` | Resampling seed, so the ensemble is reproducible. |
| `columns` | Feature name per coefficient position, intercept first. |
| `column_indices` | Index into the eight-feature model vector for each column. |
| `replicates` | `[replicate][column]` coefficients on the raw feature scale. |

A prediction under replicate `b` is `sum_j replicates[b][j] * x[column_indices[j]]`,
with the intercept column contributing 1. The reported interval is the empirical
2.5 and 97.5 percentiles of those predictions.

Storing the replicates rather than only the marginal percentiles is what makes a
**joint** interval possible: the per-coefficient intervals discard the
correlation between the intercept and the slopes, and recombining them by
interval arithmetic is conservative rather than correct.

The resulting interval is a confidence interval for the **fitted mean response**.
It is not a prediction interval and is not named like one: it excludes residual
scatter, and it carries no information about extrapolation. A candidate outside
the training range can have a narrower band than one inside it, because width
tracks the coefficient spread rather than domain membership.

The ensemble dominates the document size — one row per replicate per
coefficient. `stericx fit --omit-bootstrap-ensemble` writes the document without
it; screening then reports no bootstrap interval rather than inventing one.

This section is absent from schema 1 entirely, and the ensemble is never written
into the flattened fit report, so `--output` stays byte-identical.

### The one threshold, and where it comes from

Only the nearest-neighbour measure needs a boundary, and it is **measured, not
chosen**. For every training point, take the distance to the nearest *other*
training point; the boundary is the maximum of those distances:

```text
threshold = max over training points of ( distance to nearest other training point )
```

A candidate is inside the sampled region when it is no farther from the training
set than the training set's sparsest point is from its own neighbour. This has
no free parameter — no `0.5σ`, no percentile, nothing to tune.

Being set by the loosest part of the training set makes it **permissive**, and
that is a real limitation: one isolated training observation widens the boundary
for everything. `mean`, `standard_deviation`, and `median` of the same distance
distribution are serialized alongside it, and `stericx screen --domain-rule`
selects which of them bounds the domain:

| `--domain-rule` | Boundary | Notes |
| --- | --- | --- |
| `max-neighbor` (default) | `maximum` | No distributional assumption. Permissive, as above. |
| `mean-plus-sd` | `mean + σ` | ~84th percentile if the distances are roughly normal. |
| `mean-plus-2sd` | `mean + 2σ` | ~97.7th percentile under the same assumption. |

Every option is a statistic **of the training set's own neighbour-distance
distribution**, computed at fit time and serialized into the model. None is a
severity grade, and none introduces a cutoff chosen to make results look a
particular way; the two `mean + kσ` rules do add an assumption the maximum does
not — that the distances are roughly symmetric — which is why the parameter-free
maximum remains the default. `max-neighbor` is not always the loosest: a tight
training set with one far outlier can put `mean + 2σ` above the maximum, so the
rules are ordered by the statistic they name, not by severity.

The applied rule and the boundary it produced are reported in every screen
(`domain_rule`, `domain_rule_description`, `domain_threshold`), so a stricter run
is never mistaken for a default one. Changing the rule moves only the boundary:
the measured `nearest_training_distance` is a property of the candidate and the
training set, and is identical under every rule. The `rule` field on the stored
calibration carries the default derivation in words.

### Verdicts

No severity grades are invented. Each verdict is a stated combination of the two
checks above:

| Verdict | Definition |
|---|---|
| `interpolation` | Every descriptor inside its training range, and nearest-neighbour distance ≤ threshold. |
| `sparse_interpolation` | Every descriptor inside its range, but distance > threshold — a gap the training set did not sample. |
| `extrapolation` | At least one descriptor outside its training range. |
| `unknown` | No training geometry, or no calibration to compare against. Not a claim of either state. |

Applicability is computed from descriptors alone. `assess_applicability` takes
no prediction and there is no way to pass one, so a favourable-looking number
can never make a ligand appear more in-domain than it is.

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
  --bootstrap 2000 --permutations 2000 \
  --portable-model docs/study_001/stericx_portable_model.json \
  --model-id mechanistically_constrained_ols \
  --reaction-family "Ni-catalyzed homo-Diels-Alder" \
  --catalyst-metal Ni \
  --ligand-class "monodentate phosphorus(III)" \
  --source-url "https://raw.githubusercontent.com/SigmanGroup/Ni-Catalyzed-hDA/main/data/kraken.csv" \
  --response-temp-k 298.15 \
  --optimize maximize
```

That is the exact command that produced the checked-in
`docs/study_001/stericx_portable_model.json`. The replicate counts are not the
defaults — the Study 001 artifact was fitted with 2,000 bootstrap and 2,000
permutation replicates, and the defaults reproduce every other field but shift
the coefficient intervals and the permutation p-value.

The checked-in `stericx_model.json` predates `training_geometry.neighbor_calibration`
and `standardized_training_points`, which this build records. Re-running the command
against that path re-emits it with those two fields added; every pre-existing value,
including the bootstrap intervals and the permutation p-value, stays byte-identical.
It was deliberately left as published, so the version 1 artifact still screens with
applicability verdicts of `unknown` while the version 2 document reports
`interpolation`.

`--output` is unchanged by this flag, so existing study artifacts stay
byte-identical. Study 001 therefore publishes both documents side by side: the
version 1 `stericx_model.json` that the released studies and the Python drivers
already read, and the version 2 `stericx_portable_model.json` that carries the
response definition, the provenance, the applicability calibration, and the
optimization direction `screen` needs to rank without an explicit flag.

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
