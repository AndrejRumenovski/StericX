# StericX demonstration readiness

**Reviewed locally: 2026-09-16.** Ready for a research discussion and demonstration
of descriptor calculation, library search, and model auditing. Native ligand
ranking remains experimental and has not shown useful predictive performance.

Start with the [5–7 minute walkthrough](DEMO.md). The
[visual results document](results.html) and [full report](REPRODUCTION_REPORT.md)
now distinguish descriptor agreement from reaction prediction. Read the
[scientific interpretation notes](SCIENTIFIC_NOTES.md) before discussing the
frozen forecast or the residual mechanism.

## Rehearsal

```bash
cargo build --release
bash scripts/demo.sh --pause
# Optional retrospective fit and candidate export:
bash scripts/demo.sh --screen --pause
```

The release binary has been rebuilt in this workspace. The offline demo uses
checked-in inputs and writes fresh `.stericx/demo-*` folders containing commands
and output logs. Use a successful rehearsal's logs as a meeting fallback.

## Verification

| Check | Result |
|---|---|
| Rust suite: `cargo test --all-targets --quiet` | **248 passed**, including 21 published screening regression tests |
| Rust lint: `cargo clippy --all-targets --all-features -- -D warnings` | Passed |
| Rust formatting: `cargo fmt --all --check` | Passed |
| Optimized binary: `cargo build --release` | Passed |
| Python quantum parsers: `python -m unittest tests/test_quantum_backend.py` in the project environment | **11 passed** |
| Ruff lint/format on 23 tracked Python files | Passed |
| Study 011 `--verify-only` | Passed: 57 planned splits, design lock and frozen-prediction hashes verified |
| CLI rehearsal | Descriptors, compare, constrained search, model inspection/validation, fit, screening, and candidate export passed |
| Output integrity | JSON and CSV parse correctly when combined with deck export |
| Invalid geometry | Non-finite XYZ/SDF coordinates are rejected; regression tests cover all three axes |
| Target convention | Explicit magnitude survives fit → save → inspect; fit reports and frozen prediction bytes are unchanged |
| Portable example correction | Every JSON field except `inference.response.sign_convention` is unchanged |
| Documentation | Local evidence links checked; no missing files in the revised entry points |

Repository-wide Ruff still reports existing findings in the user's untracked
`docs/media/geom.py` and `docs/media/render.py` (13 lint findings and two files
needing formatting). Those files were left intact. The tracked-script check above
does not claim the entire working directory passes Ruff.

No browser connection was available, so the HTML document's rendered layout was
not visually verified. Its existing styling was retained, and local links were
checked. This review did not rerun the expensive quantum studies, repeat the
historical speed benchmark, or test another operating system.

## Points to be ready to explain

- **Contribution:** fast Rust implementations of established descriptor definitions,
  with comparisons against reference calculations and published data.
- **Scope:** 1:1 R² = 0.9852 is for `vbur_max_delta_qvbur_min` on matched DFT ensembles;
  it is not a catalyst-performance score. Median error and RMSE are both useful.
- **Conventions:** bond-axis versus coordination-axis Sterimol, the 2.28 Å virtual
  metal, Bondi radii ×1.17, the 3.5 Å sphere, and conformer aggregation affect results.
- **Negative evidence:** native Ni-hDA Q² is approximately zero, and Study 011's
  top-1 recovery is 0.158 versus 0.333 random, with ten failed fits retained.
- **Next scientific question:** whether consistent training/candidate conformer
  aggregation and a chemically suitable held-out dataset improve ranking. Define
  and evaluate a new benchmark rather than replacing the historical result.
- **Forecast:** the response is selectivity magnitude, with no favored-enantiomer
  assignment and no measured prospective outcomes in the repository.

The emphasis on descriptor conventions and held-out reaction performance fits the
Sigman group's stated interest in linking structural descriptors to reaction
outcomes. This is a suggested discussion focus, not a prediction of Professor
Sigman's reaction. [Sigman Lab research](https://www.sigmanlab.com/research)

This report records local checks performed before the readiness changes were
committed.
