# StericX v0.3.0 release audit

**Audit date:** 2026-09-04

**Release:** `v0.3.0`

**Decision:** PASS — eligible for tagging only after the release commit's clean-worktree
checks and GitHub Actions run complete successfully.

No scientific threshold, tolerance, split, or failed observation was removed or weakened
for this release. Study 011's below-random ranking result and all ten failed fits remain in
the release artifacts and in the changelog.

## Environment

The local audit used:

```text
rustc 1.97.0 (2d8144b78 2026-07-07)
cargo 1.97.0 (c980f4866 2026-06-30)
uv 0.11.32
Python 3.12.13
ruff 0.16.0
gh 2.96.0
```

The authoritative pre-tag repetition is run from a clean Git worktree. The primary
workspace also contains unrelated, untracked `docs/media/*.py` files; an initial
`uv run ruff check .` correctly reported 13 findings in those untracked files. They are not
part of the release, were not modified or deleted, and no Ruff exclusion was added. The
unchanged command is required to pass in the clean release worktree.

### Resolved clean-clone blocker

The first clean-worktree run stopped before release because Study 011 verification could
not find `data/reactions.sigpack`. The 704-byte file was checksum-locked as a study input but
globally ignored as a generated artifact, so the study reproduced only in a workspace that
had previously run `stericx parse`.

The gate was not weakened or skipped. The exact derived sigpack is now allow-listed and
tracked. Its SHA-256 is
`6dd2d1190127ad88b3c883d025d24cbd365099060d943cdce498f4f226d406ae`, matching the
Study 011 design lock and run manifest. Re-running `stericx parse` during the smoke test
reproduces it byte-for-byte from the already-committed reaction CSV and geometries. No split,
prediction, ranking, metric, or validation threshold changed.

## Required gates

| Area | Exact command | Result |
|---|---|---|
| Rust formatting | `cargo fmt --all --check` | PASS |
| Rust lint | `cargo clippy --all-targets --all-features -- -D warnings` | PASS, zero warnings |
| Library unit tests | `cargo test --lib` | PASS, 79/79 |
| Complete test graph | `cargo test --all-targets` | PASS, 243/243 |
| Integration tests | `cargo test --tests` | PASS, 135 integration tests plus 108 unit tests |
| Published scientific regression | `cargo test --test screening_regression -- --nocapture` | PASS, 21/21 |
| Deterministic selection | `cargo test --test cli_screen deterministic` | PASS, 3/3 |
| Deterministic screen text | `cargo test --test cli_screen repeated_runs_are_byte_identical` | PASS, 1/1 |
| Deterministic deck | `cargo test --test cli_screen identical_inputs_produce_an_identical_deck` | PASS, 1/1 |
| Deterministic formats | `cargo test --test screening_regression every_output_format_is_byte_identical_across_runs` | PASS, 1/1 |
| Model format | `cargo test --test portable_model_format` | PASS, 17/17 |
| Model CLI | `cargo test --test cli_model` | PASS, 17/17 |
| Training API compatibility | `cargo test --test model_training_api` | PASS, 6/6 |
| Python lint | `uv run ruff check .` | PASS in clean release worktree |
| Python formatting | `uv run ruff format --check .` | PASS in clean release worktree |
| Quantum parser tests | `uv run python -m unittest tests/test_quantum_backend.py` | PASS, 11/11 |
| Retrospective-study integrity | `uv run --extra science python studies/study_011_ligand_ranking.py --verify-only` | PASS, 57 planned splits and both locked hashes verified |
| Citation schema | `uvx --from cffconvert cffconvert --validate -i CITATION.cff` | PASS against CFF 1.2.0 |

The Study 011 verifier reported:

```text
Study 011 verification passed
  design lock: efae084de2df171b971dc0308674b5e93bffeb6e8d91f8dd5682357ed130f920
  frozen predictions: b7f5f34fdfbf9389ac4bf7d1cb4afcab419f3c5b7314f9c23aa818ec81051cd1
  planned splits: 57
```

## CLI-help audit

The root help page and every command page must exist, exit successfully, and contain a
`Usage:` section. This exact loop passed for the root plus 15 command pages:

```bash
./target/release/stericx --help

set -e
help_count=0
for help_case in \
  'parse' 'buried-volume' 'predict' 'fit' 'evaluate' \
  'model' 'model inspect' 'model validate' 'simulate' 'descriptors' \
  'search' 'screen' 'compare' 'db' 'db build'; do
  read -r -a help_args <<< "$help_case"
  help_text="$(./target/release/stericx "${help_args[@]}" --help)"
  test -n "$help_text"
  [[ "$help_text" == *Usage:* ]]
  help_count=$((help_count + 1))
done
test "$help_count" -eq 15
```

The release binary also passed:

```bash
test "$(./target/release/stericx --version)" = "stericx 0.3.0"
```

## Model schema and backwards compatibility

[`MODEL_FORMAT.md`](MODEL_FORMAT.md) matches the implementation and tests:

- schema 2 remains the highest supported model schema and is a strict superset of schema 1;
- schema-2 documents round-trip byte-for-byte and exact `f64` values are preserved;
- every read validates dimensions, feature names, scaling, ranges, provenance, and repeated
  inference fields before returning a model;
- unknown future versions and malformed models are rejected rather than partially scored;
- the checked-in Study 001 schema-1 model still loads, inspects, and validates with an
  explicit `legacy_schema` warning;
- the schema-1 model still screens when the operator supplies `--descending`; missing
  response context, leverage geometry, and bootstrap ensemble are reported as unavailable,
  not guessed;
- the schema-1 fit report and frozen-prediction formats, sigpack v1/v2 readers, and v0.2 CLI
  aliases remain covered by regression tests.

The explicit legacy CLI check was:

```bash
./target/release/stericx model inspect docs/study_001/stericx_model.json
./target/release/stericx model validate docs/study_001/stericx_model.json
./target/release/stericx screen docs/study_001/stericx_model.json \
  --library data/reactions_raw.csv --descending --top 3 --format json
```

Result: schema 1, `portable=false`, validation `valid=true` with the documented legacy
warning, and three ranked candidates.

## Example-command audit

The README quick-start descriptor command and v0.3 screening commands passed unchanged:

```bash
./target/release/stericx descriptors data/xyz/SIG-NIHDA-401_9d42bff1.xyz
./target/release/stericx model inspect docs/study_001/stericx_portable_model.json
./target/release/stericx screen docs/study_001/stericx_portable_model.json \
  --library data/reactions_raw.csv --top 3
```

The complete commands in [`REACTION_SCREENING.md`](REACTION_SCREENING.md) were exercised by
the release smoke test below. All repository-relative files referenced by that tutorial,
the model-format guide, and Study 011 were checked to exist.

## End-to-end screening smoke test

The test starts with the checked-in experimental reaction CSV, rebuilds its packed records,
fits only the eight observations permitted by the predefined Study 011 `S001` split,
serializes and validates a portable model, screens after excluding those eight tested
ligands, retains applicability and both uncertainty signals, makes a diverse top-2
selection, and exports an auditable deck.

```bash
cargo build --release
mkdir -p .stericx/release-0.3.0-smoke

./target/release/stericx parse \
  --csv data/reactions_raw.csv \
  --xyz-dir data \
  --output .stericx/release-0.3.0-smoke/reactions.sigpack

cmp data/reactions.sigpack \
  .stericx/release-0.3.0-smoke/reactions.sigpack

./target/release/stericx fit \
  --data .stericx/release-0.3.0-smoke/reactions.sigpack \
  --metadata docs/examples/ni_hda_screening_split.csv \
  --output .stericx/release-0.3.0-smoke/fit-report.json \
  --predictions .stericx/release-0.3.0-smoke/frozen-predictions.csv \
  --portable-model .stericx/release-0.3.0-smoke/model.json \
  --model-id ni-hda-v0.3.0-release-smoke \
  --reaction-family "Ni-catalyzed homo-Diels-Alder" \
  --catalyst-metal Ni \
  --ligand-class "monodentate phosphorus(III)" \
  --source-url \
    "https://raw.githubusercontent.com/SigmanGroup/Ni-Catalyzed-hDA/main/data/kraken.csv" \
  --response-temp-k 298.15 \
  --optimize maximize \
  --bootstrap 1000 \
  --permutations 500 \
  --seed 20260904

./target/release/stericx model validate \
  .stericx/release-0.3.0-smoke/model.json --strict

./target/release/stericx model inspect \
  .stericx/release-0.3.0-smoke/model.json --format json \
  > .stericx/release-0.3.0-smoke/inspect.json

./target/release/stericx screen \
  .stericx/release-0.3.0-smoke/model.json \
  --library data/reactions_raw.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --top 2 \
  --diverse \
  --diversity-weight 0.5 \
  --domain-rule max-neighbor \
  --format json \
  > .stericx/release-0.3.0-smoke/screen.json

./target/release/stericx screen \
  .stericx/release-0.3.0-smoke/model.json \
  --library data/reactions_raw.csv \
  --exclude-tested docs/examples/ni_hda_screening_tested.csv \
  --top 2 \
  --diverse \
  --diversity-weight 0.5 \
  --domain-rule max-neighbor \
  --export-deck .stericx/release-0.3.0-smoke/candidate-deck.csv
```

The assertions applied to those outputs were:

```bash
jq -e '
  .schema_version == 2 and
  .portable == true and
  .training_observations == 8 and
  .stericx_version == "0.3.0"
' .stericx/release-0.3.0-smoke/inspect.json

jq -e '
  .library_size == 11 and .screened == 3 and .returned == 2 and
  .skipped == 0 and .ranking_order == "descending" and
  .model_optimization == "maximize" and .domain_rule == "max_neighbor" and
  .diversity_applied == true and .diversity_weight == 0.5 and
  .exclusion.excluded == 8 and
  ([.hits[].domain_verdict] | all(. == "interpolation")) and
  ([.hits[].prediction_interval_low] | all(. != null)) and
  ([.hits[].uncertainty.method] |
    all(. == "percentile_bootstrap_mean_response"))
' .stericx/release-0.3.0-smoke/screen.json

test "$(wc -l < .stericx/release-0.3.0-smoke/frozen-predictions.csv)" -eq 4
test "$(wc -l < .stericx/release-0.3.0-smoke/candidate-deck.csv)" -eq 3
```

### Smoke-test results

- Parsing processed 11 reaction records and wrote a 704-byte sigpack exactly equal to the
  committed `data/reactions.sigpack`.
- Fitting used 8 observations in 6 groups and froze 3 predictions. It selected
  `nbo_charge`; training R² was 0.6929, ordinary LOO Q² was 0.4625, and group-LOO Q² was
  **-192.5603**. The failed group-transfer diagnostic remains prominent and prevents this
  example model from being treated as deployment-ready.
- Schema-2 model validation reported 0 errors, 0 warnings, and `valid=true`.
- Exclusion matched all 8 identifiers, leaving 3 candidates; all 3 were screened, none was
  skipped, and the diverse selector returned 2.
- The domain rule was `max_neighbor` with threshold 0.0288334. All three eligible candidates
  were interpolation cases; this is only an applicability statement, not proof that the
  ranking is accurate.
- Diverse selection returned `SIG-NIHDA-724` first (ordinary rank 1, predicted ΔΔG‡
  1.5366101 kcal/mol, 95% prediction interval [0.3893550, 2.6838652], bootstrap mean-response
  interval [1.1682507, 1.9527320]) and `SIG-NIHDA-498` second (ordinary rank 3, predicted
  ΔΔG‡ 0.4727847 kcal/mol, 95% prediction interval [-0.6121128, 1.5576821], bootstrap
  mean-response interval [0.1883410, 0.8623360]).
- The candidate deck has one header plus 2 candidates and contains no experimental-response
  field. Its sidecar records the model, library, exclusion list, ranking direction, domain
  rule, uncertainty method, diversity objective, and cryptographic hashes.

The JSON screen was repeated and compared with `cmp`. The deck was exported twice to the
same path, with SHA-256 taken before and after. Both comparisons were byte-identical for the
fixed serialized model. The recorded smoke-run hashes were:

```text
candidate-deck.csv       0e0f420f37079ba2533d4f91674eb319a1d3bc9801d03d7832cd2afc5ea7352b
candidate-deck.meta.json 9b73d33fa21c07c13a35e0b47d6a85937fb4a0fdbe75907f9c7e60ee065ae2fd
```

These hashes are deliberately scoped to this serialized model. A new `fit` records a new
creation timestamp, so a separately fitted model and decks that embed its hash need not have
the same bytes. Screening determinism means identical model and library bytes plus identical
options produce identical output.

## README, changelog, and release metadata

- `Cargo.toml`, `Cargo.lock`, `CITATION.cff`, `.zenodo.json`, the README citation, and the
  binary all report `0.3.0`; the citation release date is 2026-09-04.
- `pyproject.toml` remains `stericx-tools` 0.1.0. It is the independently named,
  unpublished support-script environment rather than the StericX Rust package, and it was
  already 0.1.0 in the v0.2.0 release; changing it would misrepresent a separate package
  release.
- `CHANGELOG.md` leads with one v0.3.0 section and records additions, compatibility, and
  known scientific and data-licensing limitations.
- `RELEASING.md` includes the screening regression and Study 011 verification commands,
  requires a clean-worktree audit, and retains the tag → GitHub release → Zenodo order.
- `CITATION.cff` validates against schema 1.2.0. `.zenodo.json` is valid JSON and describes
  the v0.3.0 screening additions without changing the project identity or references.

Metadata checks:

```bash
cargo metadata --no-deps --format-version 1
./target/release/stericx --version
jq empty .zenodo.json
uvx --from cffconvert cffconvert --validate -i CITATION.cff
git diff --check
```

## Added-data licensing and provenance

No new third-party raw data were added between v0.2.0 and v0.3.0. One existing derived
binary representation is newly tracked so a clean clone contains every locked Study 011
input.

- Study 011's tables, JSON, and plots are StericX-generated results derived from the existing
  11-row Ni-hDA subset. Its `run_manifest.json` pins every source input and generated artifact.
- `data/reactions.sigpack` is a 704-byte StericX-generated packing of the already-committed
  `data/reactions_raw.csv` rows and geometry-derived descriptors. The release smoke test
  rebuilds it and requires byte equality before fitting.
- The tutorial split and exclusion CSV contain row identifiers, split roles, and existing
  scaffold-group labels needed to replay predefined split `S001`; they introduce no new
  experimental outcomes.
- `data/official/provenance.json` records the upstream raw URL, retrieval timestamp, row and
  response counts, and SHA-256
  `4af1a776378a1f1aa369e076c9b35de2afa8c60e30841ec63cd13d95ccb8f00d`.
- The official upstream GitHub API reports `SigmanGroup/Ni-Catalyzed-hDA` as public and
  unarchived, with no declared SPDX license. Public availability is therefore cited as
  provenance, not characterized as a license grant.
- Copyrighted AAAS supporting information used by earlier local studies remains excluded.
  The new screening regression suite pins only committed StericX outputs for Studies 007 and
  009 and explicitly states that it does not reproduce their copyrighted per-ligand inputs.
- StericX source and original generated outputs remain offered under `MIT OR Apache-2.0`;
  both license texts are present. The archive metadata selects MIT while `Cargo.toml`, the
  README, and `CITATION.cff` preserve the dual-license statement.

The provenance audit commands included:

```bash
git diff --name-status v0.2.0..HEAD -- data docs/examples docs/study_011 tests/data
gh api repos/SigmanGroup/Ni-Catalyzed-hDA \
  --jq '{html_url,visibility,license:.license.spdx_id,default_branch,archived}'
uv run --extra science python studies/study_011_ligand_ranking.py --verify-only
```

## Release conditions

The annotated `v0.3.0` tag and GitHub release may be created only when:

1. every command in the required-gates table passes from the release commit in a clean
   worktree;
2. the end-to-end smoke test and assertions pass there;
3. the push-triggered GitHub Actions workflow completes successfully for the release commit;
4. `v0.3.0` does not already exist locally, remotely, or as a GitHub release.

If any condition fails, the release stops. A scientific failure must be fixed in the method
or reported as a limitation; its threshold must not be loosened to force a pass.
