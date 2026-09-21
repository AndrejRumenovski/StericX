# Replaying independent Kraken comparisons

Use the frozen inputs and audit environment described in the master manifest. A native replay must finish and be frozen before interpreting references. No network requests are required for independent calculation against already frozen scientific inputs.

Scripts resolve `ROOT` as the parent of their `scripts/` directory. To replay safely, copy the audit tree to a separate destination preserving that layout, or archive generated phase directories first. Do not overwrite or delete the original frozen phases. The following phase directories are created with `exist_ok=False` and must be absent at their creation step:

1. `analysis`
2. `morfeus_outliers`
3. `full_analytic_reference`
4. `interpretation`
5. `focused_diagnosis`
6. `delta_outliers`
7. `delta_interpretation`

Keep `primary/`, `raw_sources/`, `primary_si/`, `prepared_all/` and the frozen `sut/` outputs available. `morfeus_outliers_attempt1/` is an explicitly invalid reference-adapter attempt and must not replace the corrected `morfeus_outliers/` data. Archive all seven generated directories together before starting, rather than letting later scripts accidentally read old analysis alongside newly computed references. Preserve final reports/claims too if comparing narrative regeneration.

From the repository root, after creating a separate replay copy or archiving the above phase directories:

```bash
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_analyze.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_morfeus_outliers.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_full_analytic_reference.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_interpret.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_matched_percent.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_focused_diagnosis.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_pyr_rounding_bound.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_delta_outliers.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_finish_delta.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_369_dossier.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_write_report.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_claims.py
uv run --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_finalize.py
```

The full analytic and Morfeus-outlier reference computations can run independently after `analysis`; the displayed sequential order is simplest. The focused rounding study runs 1,028,096 corner calculations. The delta campaign computes all missing reference records and deliberately reuses only records from the already frozen corrected Morfeus campaign. `kraken_finish_delta.py` creates the authoritative combined dossier under `delta_interpretation/`; it preserves the earlier `interpretation/` dossier unchanged. Always use the completed phase's manifest to establish that all records were written before interpreting its values.

Compare numerical rows, input IDs/hashes, error records and statistics. UTC timestamps and absolute acquisition/replay paths differ by design. Do not require timestamp-containing manifests to hash identically after a replay. The frozen reference software version is Morfeus 0.8.0; use the pinned environment and frozen source inventory. A different Morfeus version is a new scientific comparison, not an exact replay.

`kraken_energy_sources.py` documents optional fresh network checks. Its destination `energy_source_search/` must likewise be archived first: prior unsuccessful responses are immutable and intentionally cannot be overwritten in place. New source requests may give different HTTP statuses or data because the remote API/repositories can change. The audit conclusion is based on the original response bytes, not on assuming present-day requests reproduce their hashes.
