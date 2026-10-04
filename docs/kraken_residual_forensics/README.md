# Kraken residual forensics

Read [KRAKEN_RESIDUAL_FORENSICS.md](KRAKEN_RESIDUAL_FORENSICS.md) for the findings, [REQUIREMENTS_AUDIT.md](REQUIREMENTS_AUDIT.md) for coverage, and [REPRODUCE.md](REPRODUCE.md) for the original experiment workflow.

The original 1,541 ligands and 31,611 conformers remain in every comparison. The investigation distinguishes the historical-reference and current-API results; no target fitting or outlier removal was used. Unavailable historical inputs and unresolved causal assignments are explicit.

## Restore the published evidence

Reports, key metrics, plots and authored scripts are directly readable in GitHub. Larger evidence is stored losslessly as ordered archive parts of at most 48 MiB, following this repository's existing scientific-audit publication convention. From the repository root:

```sh
python3 docs/kraken_residual_forensics/publication/manage.py restore
python3 docs/kraken_residual_forensics/publication/manage.py verify
```

Restoration verifies every part and the combined archive, rejects unsafe paths, preserves executable modes, and refuses to replace an existing file with different bytes. Verification checks every published file against both the publication manifest and the unchanged scientific `BUNDLE_MANIFEST.json`, plus the recorded repository source and accepted sampling change. It does not rerun scientific calculations.

Downloaded third-party source captures remain local under the existing repository publication policy. `publication/manifest.json` identifies every exclusion, its original path, size, hash and reason. Their acquisition metadata and source URLs are published. The verifier reports missing local-only captures explicitly; its success means published-evidence integrity, not that every original source capture or unavailable historical input is present. Numerical and geometry inputs, outputs, failed attempts and derived results are preserved in the publication.

All original scientific files, including the report, sealed manifests and audit-local `.gitignore`, retain their bytes. Root Git visibility rules, this README, the style configuration and publication tools are separate packaging additions. The sealed report's links into archived directories become available after restoration. Scratch checks may use `--audit-dir /path/to/copy` and `--repo-dir /path/to/repository`.
