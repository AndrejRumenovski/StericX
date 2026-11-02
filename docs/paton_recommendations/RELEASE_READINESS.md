# Release readiness after the Paton implementation

The development distribution is built and tested. **No release candidate, tag,
GitHub release or PyPI publication has been created.** The engine remains 0.3.0;
the new Python distribution is 0.4.0.dev0. This preserves the version history
without implying that a gated release has been admitted.

| Required evidence | Result | Scope / receipt |
|---|---|---|
| Rust tests and quality checks | Pass | 348 local tests; Clippy, formatting, rustdoc; actual Linux/Windows/macOS Rust CI |
| Research-tool Python tests | Pass | 128 tests |
| Python installed-package tests | Pass | 71 tests per tested environment, including typed results, errors and optional adapters |
| Cross-platform wheels | Pass | 12 Linux/Windows/macOS × CPython 3.11–3.14 jobs; Linux/Windows x86_64 and macOS arm64 |
| Source installation | Pass | Clean local installation and actual CI source-install job |
| Legacy scientific preservation | Pass | 128,021 exact records in 13 lanes; 93 CLI cases at four thread counts |
| Independent legacy references | Preserved | Fresh Morfeus/Kraken replay, all regional bins and previous negative findings |
| Conformers and populations | Scoped pass | 56 conformers, 11 ligands; 42 valid and 11 invalid independent population cases |
| DBSTEP Sterimol | Scoped pass | 61 cases under explicit radii, axial and angular conventions |
| DBSTEP hard-sphere volume | **Strict gate false** | One of 69 matched-lattice cases fails at R=2.5 Å in three frames |
| Density volume | Scoped pass | 12 independent cases; explicit scalar density, threshold, units and supported grid geometry |
| Radius profiles | Partial strict agreement | All positions preserve individual calls; density passes; known hard-sphere witness remains |
| Numerical convergence | Scoped pass / retained limits | Analytic Gaussian convergence and existing integration/precision studies; no universal continuum claim |
| Performance | Preserved within stated protocol | 504 confirmed-window samples; no >10% median slowdown in 28 configurations |
| Documentation and tutorials | Prepared | Descriptor-focused README; executed notebook; scientific and performance evidence |

The hard-sphere witness is a demonstrated numerical integration difference:
native f32 sphere inclusion retains 5,565 nodes, while DBSTEP's f64 reclip retains
5,533. The discrepancy is not removed by a new tolerance, radius adjustment or
production arithmetic change. Its original strict comparison remains false.
Consequently this work does **not** issue the requested all-gates-passed release
candidate. The independently validated density mode is a separate descriptor and
cannot repair that hard-sphere gate.

Any later release decision must explicitly resolve the applicable comparison
contract with independent evidence while retaining this witness. A production
correction, if scientifically justified, requires a changed-output ledger and
fresh scientific, packaging and performance verification for the changed source.
Historical Kraken provenance gaps and experimental catalyst-prediction failures
also remain outside any descriptor-release claim.

The actual wheel matrix is bound to commit
`aa78ecb09d337488255937f6eb85e9027f0984db`. Subsequent commits add audit/reporting
adapters, documentation and evidence; they do not change the tested scientific
kernels, Python implementation or binding configuration.
[Distribution receipts](python/DISTRIBUTION.md) retain the exact
artifact hashes and initial failures. PyPI name availability was checked, but a
404 response reserves nothing; a future publication must recheck availability
and verify installation from the published registry in a clean environment.

See [scientific validation](SCIENTIFIC_VALIDATION.md),
[performance](performance/PERFORMANCE.md), and
[the implementation report](FINAL_REPORT.md) for the detailed scope and limits.
