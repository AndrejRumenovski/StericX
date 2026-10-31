# Conformer population validation

The new Rust population API passes the independently specified
[numerical and provenance protocol](POPULATION_VALIDATION_PROTOCOL.md): **42
valid numerical cases and 11 invalid-input cases**, including all four energy
units and all three energy interpretations. Four supplemental Rust tests cover
closed-form ratios, explicit versus grouped degeneracy, temperature response,
nonfinite values and thermodynamic metadata. Six ensemble-statistics tests also
pass after the normalization change.

| Independent measure | Initial implementation | Compensated normalization |
| --- | ---: | ---: |
| Largest absolute weight error | 9.717701413156056e-17 | 6.049964658401754e-17 |
| Largest relative error (weights ≥1e-300) | 1.5888016734221925e-14 | 1.5888016734221925e-14 |
| Largest normalization error | 7.004859982806446e-15 | 7.632783294297951e-17 |
| Failed cases | 0 | 0 |

The initial 4,096-state normalization nearly consumed the predeclared absolute
error budget. Compensated summation reduces accumulation error without changing
the statistical-mechanical expression, constants, thresholds or legacy engine
calculations. Both runs are retained; the tolerance was not relaxed. This is
arithmetic validation, not evidence that electronic or force-field weights are
experimental equilibrium populations. No historical conformer energies were
created or inferred.

The reference uses 80-digit Decimal arithmetic, direct exponentiation and exact
SI constants independently of the Rust shifted-log implementation. Raw energy
bits accepted by Rust are recorded; no JSON energy-rounding differences occurred
in these cases. Free energies at a different temperature, missing provenance,
mixed energy models/kinds/units and invalid degeneracy are rejected. Equal
energies do not cause conformers to be removed or merged. Extreme arithmetic
cases include nonphysical temperatures solely to test numerical behavior.

Reproduction (choose fresh output directories):

```bash
cargo test --locked --test conformer_populations --test ensemble_statistics
cargo build --locked --example population_observer
.venv/bin/python scripts/validate_conformer_populations.py \
  --observer target/debug/examples/population_observer \
  --output docs/paton_recommendations/conformers/population_validation_v2
```

The [initial report](population_validation_v1/report.json) and
[final report](population_validation_v2/report.json) retain every comparison.
Their receipts include the pre-execution protocol hash, source and executable
hashes, Python version, command, raw requests, responses and stderr. The initial
Rust test compilation required an explicit f64 suffix on a test-only literal;
`step3_tests_initial.log` preserves that failed build. Subsequent test and build
logs retain successful execution. Cross-platform libm behavior is still a later
distribution gate; these results were obtained on the recorded Linux host.
