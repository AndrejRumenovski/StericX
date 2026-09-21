# Independent review of the completed scientific audit

The renewed request was checked against the current StericX workspace on 21 September 2026. The current system is identical to the system examined in the completed [scientific accuracy audit](../../SCIENTIFIC_ACCURACY_AUDIT.md): commit `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7`, executable SHA256 `b577c49d3e98b4745213fce0e70e93a55fc61120c5abe6fdb0658547c43cab94`.

No production code, scientific input, frozen observation, existing claim classification or sealed report was changed. This review is additive. It found no material new gap in the investigated requirements; that does not convert existing incorrect, uncertain or out-of-scope claims into successful validation.

## Current identity and preserved evidence

[Current identity](current_identity.json) records the commit, executable hash, Rust/Cargo/Python versions, complete installed package versions and all 24 Morfeus source hashes. They match the initial snapshot. The tracked worktree has no changes.

[Package verification](sealed_package_verification.json) checked all 74,298 sealed artifacts and 37,593 original-evidence checks, without mismatches. The original final-manifest hash remains `c31b655641a352402b14dc2d4535261a9de9c9004979c5e571d78e45cac82f03`.

## Independent coverage review

Two additional read-only reviews inspected evidence and methods rather than accepting the coverage checklist as proof. Their bounded checks found:

- **Geometry:** 27,542 descriptor comparison rows retain absolute and relative errors. Recalculated base Sterimol MAE and maximum errors agree with the reported metrics. There are 58,176 individual quadrant/octant comparisons across orientations and both independent/reference calculations. Six representatives each have 100 random rotations plus X/Y/Z rotations, with translations. The historical reconstruction contains 20 actual conformers and explicit SDF connectivity.
- **Conformers and kinetics:** the saved campaign contains 195 kinetic requests, independent high-precision equations, actual frozen conformer filtering/weighting observations and native packed-record comparisons. Synthetic weighting checks are correctly separated from conformer-search completeness and physical population validity.
- **Kraken statistics:** all 448 requested statistic fields across 56 comparisons were independently recomputed from 86,435 saved rows. [The reproducible standard-library calculation](statistics_recheck.json) agrees with the report. A separate reviewer calculation agreed within 9.9×10⁻¹⁵. The complete dossier inventory has 1,120 entries, including both extrema for all 280 conformer-range entries.
- **Models and screening:** a separate calculation from frozen Study 011 rows confirms 92/141 interval coverage, 70/108 interpolation coverage, 47 successful splits, 9/57 top-one accuracy and 0.4649122807 top-two overlap. The saved native outer-fold observations contain seven failures. All 72 classifier-fit records retain the complete data and agree with the independent reference predictions.
- **Source attribution:** frozen primary Kraken code and SI support the raw-vector center, hydrogen-radius and 0.001 integration conventions. The official threshold notebook uses an inclusive activity boundary. The reports distinguish those convention differences from unresolved historical publication/ensemble discrepancies.

The reviews found the independent geometric/reference methods, conditional validation interpretations and separation of optimization, implementation, methodological and predictive validity appropriately scoped. They did not rerun all molecular reference kernels; the original experiments and exact native replays remain in the sealed package.

## Limits remain unchanged

The original Verloop/Radhakrishnan and some other historical full texts were not fully accessible. Current Kraken exports do not provide the historical per-conformer energies and retention records needed to reproduce every original ensemble reduction. Corrected experimental target conflicts and some publication-version discrepancies remain unresolved. Numerical agreement, synthetic weighting and geometric applicability checks do not establish prospective experimental prediction.

Use the [147-claim inventory](../../CLAIMS.md), [full report](../../SCIENTIFIC_ACCURACY_AUDIT.md) and [original reproduction instructions](../../REPRODUCE.md) for the scientific findings and proposed corrections. This review does not authorize or implement those corrections.

## Reproduce this verification

From the repository root, using the recorded environment:

```sh
python3 docs/scientific_accuracy_audit/scripts/audit_package.py verify
.venv/bin/python3 docs/scientific_accuracy_audit/rechecks/recheck_current_system.py
```

The second command creates a new timestamped directory and preserves earlier records. It checks current identity and independently recomputes all Kraken comparison statistics from frozen rows. [The review manifest](review_manifest.json) hashes the additive review artifacts and its script; it leaves the original sealed manifest intact.
