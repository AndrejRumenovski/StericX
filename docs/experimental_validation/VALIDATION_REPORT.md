# Integrated scientific regression and independent validation

All accepted production changes are explicit APIs. The baseline public/default
behavior is reproduced exactly; independent analytical/reference checks separately
validate the new methods. Deterministic replay alone is not independent validation.

| Gate | Fresh result | Scope / limitation |
|---|---|---|
| Rust 1.97 all targets/features | 385 passed, zero failed/ignored | Includes all new public APIs and exact-addition oracle |
| Actual Rust 1.85 minimum | 385 passed, zero failed/ignored | Same frozen source and scientific test inputs |
| Format, Clippy, rustdoc | All pass | Clippy/rustdoc deny warnings |
| Python binding | 71 passed | Newly built wheel, isolated installed path |
| Python research | 156 passed, no skips | Includes conformer/thermodynamic and reference-controller safeguards |
| Upstream DBSTEP | 413 passed, three skips in isolated environment | Two RDKit-dependent skips subsequently exercised; all three API-module tests pass with RDKit available. PyMOL plugin unavailable; no GUI work attempted |
| Cross compilation | macOS AArch64 and Windows x86-64 library checks pass | Available Rust 1.98.1 targets; not physical hardware execution |
| Historical exact output | 128,021 / 128,021 unchanged, thirteen lanes | Includes errors, serialized numbers, ordering, all regional bins, model, aggregation and kinetics results |
| Original geometric analytical oracle | 167 / 167 valid enclosures | All eighteen historical failed targets now met by explicit geometric method |
| Exact-addition oracle | 6,005 rational witnesses pass | Independent Python Fraction enclosures, including subnormals/cancellation |
| Extended geometric campaign | 180 observations, zero invalid comparisons | 96 analytical unions and 84 molecular convergence requests; molecular consistency is not exact truth |
| New regional oracle | 1,248 / 1,248 integrals enclosed in 96 cases | Four intentionally stringent 1e-8 Å³ budget targets remain unmet with valid bounds |
| New radius profiles | 936 / 936 analytical integrals enclosed | 24 profiles, 120 requested entries, 96 distinct integrations; all requested ordinary targets met |
| Experimental molecular regional frames | 405 / 405 meet 0.05 Å³ | Represented local spheres; does not certify source geometry/populations |
| Fresh experimental Morfeus/DBSTEP comparison | 135 / 135 pass | Matched geometry/radii with explicit angular discretization bounds |
| Independent statistical fit check | Sixteen univariate outer folds agree within 1e-11 kcal/mol | SciPy centered covariance versus NumPy least-squares implementation |
| Fresh complete reproducibility replay | Exact numerical/model agreement | 167 analytical, 96 regional, 24 profile requests; published/generated native data; all locked model results |
| Memcheck | Zero errors, zero lost bytes in tested regional/profile processes | 544 reachable runtime bytes per process; rejected preload-path probe retained |

The historical geometry/Kraken reference controllers were freshly executed with both
inferred and recorded topology. They preserve all raw reference discrepancies and
unavailable results. The inferred-topology lane retains eight connectivity
ambiguities causing sixteen descriptor-operation refusals; these are historical
behavior, not newly deleted outliers. Matched same-center analytical L/B5 references
remain within maximum discrepancies 1.314e-6/6.984e-7 Å on 31,713 inferred-topology
successful conformers. All 31,721 conformers remain present. The explicit-topology
lane is separately retained rather than silently substituting its geometry policy.

The fresh geometry metrics, Kraken coverage, full analytical/pyramidalization metrics
and matched-default percentage rows agree with the previously admitted campaign.
Source-manifest hashes necessarily differ and are compared as provenance, not
scientific values. Raw default-versus-default library comparisons include known
radius/frame/angular conventions and are not advertised as zero residuals. The full independent finite-grid reference was also freshly recomputed for all
31,721 conformers, 95,163 frames and 1,141,956 bins, with zero reference errors;
all-bin residual metrics are compared with the previously admitted evidence. New
continuous regional integrals additionally have the independent rational oracle above.

## Failed attempts and corrections

- Missing system Valgrind during baseline inventory: record availability instead of
  treating an absent executable as a numerical failure. A previously retained local
  Valgrind distribution was later located and explicitly hashed.
- First generic replay adapter did not expose the new certified types; the second
  attempt encountered the production-private error constructor. A local, unused
  certified error bridge compiles the original private prefix without changing the
  historical observation path. Both failed builds remain archived.
- First engineering snapshot excluded test fixtures. Clippy/unit/MSRV builds failed
  to load them. A complete selected scientific/test snapshot then passes all checks.
  An overly broad 2.1-GiB archive extraction was stopped; its partial files/log remain.
- Initial regional nonnegative-contract and small-ball percentage test mistakes,
  and early syntax errors, remain in method-specific directories. Percentage test
  tolerance was tightened, not relaxed; original historical tolerances never changed.
- The first generated chemistry preparation lacked complete axial control; all
  candidates, energies and model outcomes remain, alongside the corrected lane.
- P-chirality exploratory interpretation omitted implicit-neighbor parity. Independent
  source embeddings resolved it; the final audit and pre-fit amendment are explicit.
- The first optimization-sensitivity summary compared unequal retained sets. The
  corrected summary uses matching identifiers; both summaries/controllers remain.
- A pandas read-only array stopped the first supplementary stereo fit before output;
  copied arrays fix the controller. The registered primary results are unchanged.
- Valgrind's initial space-containing preload path invalidated interception. The
  whitespace-free rerun passes without suppressions. The initial checker also missed
  one named timing field; numerical fields were never dropped to obtain agreement.
- Hardware `perf` remains unavailable. Callgrind counts are labeled simulated.

Evidence directories and controller commands are documented in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md). `campaign_manifest.json` seals required
reports, source, bundles, results and selected raw evidence. The final audit refuses
missing gates or altered baseline historical outputs. Unrelated `docs/media/` remains
untouched. No release, deployment, package publication or GUI change was made.
