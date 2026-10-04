# Reproducing this investigation

The immutable cohort is the 1,541 IDs in `frozen/cohort.json`, with 31,611 conformers. No experiment chooses a constant, input, conformer or reference version by goodness of fit. Historical and current-API references are separate named comparisons.

Run scripts from the repository root using the existing `.venv`. The scripts use NumPy, pandas, SciPy, RDKit, matplotlib, requests and Morfeus 0.8.0. Rust dependency resolution is pinned by the root and adapter `Cargo.lock` files. The validated build records its compiler and source manifest. Binary copies are retained under `frozen/`, `bin/` and `candidates/angular_sampling/validated_build/bin/`.

Every experiment refuses to replace its existing output directory. To repeat a campaign, use an isolated checkout/output root and restore the frozen source/input archives as appropriate. Do not delete the evidence to make a script run. Several scripts intentionally use the repository's earlier scientific audit assets; their paths and hashes are recorded in manifests. The evidence bundle is reviewable locally, not a promise that the multi-gigabyte earlier audit is unnecessary for a clean-room rebuild.

The execution order was:

1. `scripts/freeze.py`, then `scripts/analyze.py`, and `scripts/convention_plan.py`.
2. `scripts/run_native.py --stage baseline_replay`, `--stage raw_center`, `--stage primary_radii`, `--stage primary_grid` (six workers). The baseline CSV-parser assertion failure is preserved; `campaigns/baseline_replay/replay_verification.json` checks the raw native JSON and IEEE values directly.
3. `scripts/reference_campaign.py --workers 6`. The interrupted original stream remains. The successful invocation adds `--name morfeus_primary_resumed --reuse-prefix docs/kraken_residual_forensics/campaigns/morfeus_primary/observations.jsonl.gz`; only complete ordered JSON records are reused.
4. `scripts/minimal_reference.py`, `scripts/provenance_analysis.py`, `scripts/graph_diagnosis.py`, `scripts/frame_diagnostics.py`, `scripts/geometry_bounds.py`, `scripts/convergence.py`, `scripts/summarize_convergence.py` and `scripts/check_available_sources.py`.
5. Build the optional production sampling API adapter with `cargo build --release --manifest-path docs/kraken_residual_forensics/adapter/Cargo.toml`. Run `scripts/sampling_campaign.py`. Its plan and independent three-sphere witness predate the production edit.
6. `scripts/recover_xyz.py`, `scripts/recovered_profile.py`, `scripts/extend_minimal.py` and `scripts/consolidate.py`. The 159 extra extremizers are listed in `minimal_reference/final_top100_additional_ids.json`; all 27 changed XYZ inputs are independently recalculated.
7. `scripts/historical_dependencies.py`, `scripts/selection_boundary.py`, `scripts/report_evidence.py`, `scripts/family_plots.py`, and `scripts/historical_target.py`. The latter explicitly scores the frozen historical CSV after detecting API reference-version changes for 821 and 1036. It does not change any predictions. The failed first boundary probe and corrected source-profile probe are both retained. The family plot v2 only adds title spacing; original plots remain.
8. `scripts/historical_coverage.py` adds 15 independent cases so both target versions' baseline/final top100 extrema are covered (5,491 total). `scripts/reference_precision.py` distinguishes actual reference revisions from decimal serialization changes using the source literals. `scripts/write_report.py` renders the human report; the final verification and manifest inventory seal the evidence.

The final report and requirement audit identify the authoritative output for each question. `scripts_history/` preserves executed script bytes under their SHA-256 names, including two exactly hash-matched reconstructed versions preceding CSV-parser and interrupted-stream recovery changes. These are execution-history records, not new scientific variants. `BUNDLE_MANIFEST.json` inventories final evidence except rebuildable compiler/cache files and itself; `verification.json` records consistency checks.

Validation after the production API edit:

```sh
cargo test --all-targets --all-features
.venv/bin/python -m unittest discover -s tests -p 'test_*.py'
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo doc --no-deps --all-features
```

The existing scientific replay controller also built and observed all 13 frozen oracle lanes. Its exact manifests and all raw observations are in `candidates/angular_sampling/validated_build/` and `full_observations/`; `scripts/compare_candidate_replay.py` checks 128,021 records including private frames/bins and IEEE values against the accepted C2 receipt. New option-specific tests are in `tests/sterimol_sampling.rs`.

All reported regressions use `R² = 1 − Σ(prediction − reference)² / Σ(reference − mean(reference))²`. Pearson r, slope and intercept are diagnostics only. They never change predictions. Read CSV floats with `float_precision="round_trip"` when checking exact numeric equality. Null historical energies, RMSD and membership mean unavailable data, never zero disagreement.
