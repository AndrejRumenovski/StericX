# Reproduce the scientific accuracy audit

Start with the repository at commit `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7` and this complete audit directory under `docs/scientific_accuracy_audit/`. Retain the frozen captures. Run commands from the repository root. The audit intentionally makes no production scientific changes.

## 1. Verify the package and original evidence

```sh
python3 docs/scientific_accuracy_audit/scripts/audit_package.py verify
```

This checks the final inventory, original executable/source/reference snapshots, frozen scientific observations, acquisition hashes and unchanged production worktree. Software/tool versions are in `manifest_initial.json`; input/output hashes are also recorded in each component's original phase manifests. The final verification log is `reproducibility/package_verification.json`.

The initial source snapshot is immutable. The model acquisition journal was extended by two source downloads after its first manifest; `models/provenance_amendment_downloads.json` records this and points to the exact original journal bytes. Original scientific payloads were retained. Earlier audit-harness mistakes are preserved separately and discussed in component reports; use the corrected reference outputs, not an `attempt1` directory, for the final comparisons.

## 2. Re-run the entire frozen native observation streams

```sh
python3 docs/scientific_accuracy_audit/scripts/replay_frozen_observations.py
```

Each lane runs the frozen observer in a fresh process, enforces a ten-minute process limit, hashes the complete output stream, and compares it with the original capture. Repeated matching output files are removed after comparison; mismatches and stderr are preserved. New timestamped replay logs are stored under `reproducibility/`. The replay is a reproducibility check only; it is not independent scientific validation.

For a separate scratch replay of the final Kraken interpretation using the already frozen independent numerical references:

```sh
uv run --frozen --extra science python docs/scientific_accuracy_audit/scripts/replay_kraken_interpretation.py
```

This compares 19 generated artifacts, preserves a timestamped replay log, and leaves original outputs untouched. It does not rerun Morfeus kernels; section3 describes full reference recalculation.

To rebuild the observation adapter from the frozen source without overwriting its frozen executable:

```sh
cargo build --offline --release \
  --manifest-path docs/scientific_accuracy_audit/observer/Cargo.toml \
  --target-dir .stericx/audit-observer-rebuild
```

An offline build requires the dependencies already cached. Otherwise omit `--offline` while keeping the recorded lockfile. Compiler/platform changes may alter executable bytes; scientific output replay remains separately testable. The visibility wrapper is a byte-identical copy of frozen BV source followed by observation-only functions. Do not use it as an independent algorithm.

## 3. Recompute independent references

Use the repository's frozen `uv.lock` and science extra:

```sh
uv sync --frozen --extra science
```

To avoid changing the sealed result files, create a scratch checkout and copy the complete audit there. Run this block from the original repository root; subsequent regeneration commands run in the scratch checkout:

```sh
audit_replay_root=$(mktemp -d "${TMPDIR:-/tmp}/stericx-audit-replay.XXXXXX")
git worktree add --detach "$audit_replay_root" 6393aafe0d983e504baf8abc1e18dc2a0f0d40e7
cp -a docs/scientific_accuracy_audit "$audit_replay_root/docs/"
cd "$audit_replay_root"
uv sync --frozen --extra science
```

Geometry and model scripts write derived results in their component directory. Their guarded freeze routines do not authorize replacing the original captures. After regeneration, compare scientific values, not acquisition timestamps.

| Component | Reference regeneration entry points | Detailed instructions |
|---|---|---|
| Geometry | `scripts/geometry_run.py compare`, `geometry_summarize.py`, focused/numerical/alignment scripts | [geometry/REPORT.md](geometry/REPORT.md) and [METHODS.md](geometry/METHODS.md) |
| Kinetics / conformers | `uv run --frozen --extra science python docs/scientific_accuracy_audit/scripts/kinetics_audit.py analyze` | [kinetics/KINETICS_CONFORMERS.md](kinetics/KINETICS_CONFORMERS.md); checks frozen inputs/outputs before evaluating independent equations |
| Kraken | `kraken/scripts/kraken_analyze.py` and final interpretation/reference scripts | [kraken/REPORT.md](kraken/REPORT.md); uses complete downloaded geometries and fixed source conventions |
| Models / screening | `scripts/models_*.py` in the documented component order | [models/REPORT.md](models/REPORT.md); frozen SUT artifacts are observations and NumPy/SciPy/sklearn are independent references |

The Kraken scripts intentionally refuse an existing output directory. **Only in the scratch checkout just created**, archive those seven derived directories before regenerating them. Keep `primary/`, `raw_sources/`, `primary_si/`, `prepared_all/` and `sut/` intact: they contain the scientific inputs and frozen native observations.

```sh
python3 - <<'PY'
from pathlib import Path
import shutil
root = Path('docs/scientific_accuracy_audit/kraken')
archive = root / 'reproduction_archive'
archive.mkdir(exist_ok=False)
for name in ['analysis', 'morfeus_outliers', 'full_analytic_reference',
             'interpretation', 'focused_diagnosis', 'delta_outliers', 'delta_interpretation']:
    shutil.move(str(root / name), str(archive / name))
PY
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_analyze.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_morfeus_outliers.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_full_analytic_reference.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_interpret.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_matched_percent.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_focused_diagnosis.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_pyr_rounding_bound.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_delta_outliers.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_finish_delta.py
OPENBLAS_NUM_THREADS=1 uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_369_dossier.py
uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_write_report.py
uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_claims.py
uv run --frozen --extra science python docs/scientific_accuracy_audit/kraken/scripts/kraken_finalize.py
uv run --frozen --extra science python docs/scientific_accuracy_audit/scripts/plot_all_kraken_metrics.py
```

The first command reconstructs all 56 published-value comparisons from frozen observations. The next commands recompute independent Morfeus/analytic references, investigate alignment and coordinate rounding, expand to complete ensembles for every buried-volume range outlier, and map both extrema for each range. This takes longer than merely reinterpreting saved reference files. The original derived results remain in `reproduction_archive/` for comparison. [Kraken replay notes](kraken/REPLAY_NOTES.md) explain dependencies and the separate optional source-search regeneration. All 56 residual plots are indexed in [visuals/kraken_all_metrics/INDEX.md](visuals/kraken_all_metrics/INDEX.md).

The main science environment was Python3.12.13, NumPy2.5.1, SciPy1.18.0, scikit-learn1.9.0, RDKit2026.3.4 and Morfeus0.8.0. Exact installed distributions and Rust/compiler versions are recorded. CREST thermodynamics tests inject frozen realistic energies/tables into the frozen Python routine; they do not claim a new external CREST run or validate quantum chemistry physically.

## 4. Fresh acquisitions or new experiments

Use a new experiment directory and a new manifest. Do not overwrite the September19 input/source snapshot or raw observations. API responses may change; the audit records the captured content and its provenance rather than asserting future endpoints will return the same corpus. Failed HTTP responses and absent scientific data are part of the evidence.

The overview is [SCIENTIFIC_ACCURACY_AUDIT.md](SCIENTIFIC_ACCURACY_AUDIT.md), the individual classifications are [CLAIMS.md](CLAIMS.md), and coverage against the original24 requirements is [COMPLETION_CHECKLIST.md](COMPLETION_CHECKLIST.md).
