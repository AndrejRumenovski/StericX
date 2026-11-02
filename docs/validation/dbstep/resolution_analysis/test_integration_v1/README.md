# Convergence regression discovery correction

The final project-wide `unittest` discovery initially failed because the new
convergence tests imported pytest. Their module-level functions also would not
have been discovered by unittest. That failed execution remains retained at
`docs/validation/dbstep/resolution_validation/final_python_tests.stderr.log`.

All seven tests now reside in `unittest.TestCase`; no pytest dependency is
introduced. The test names and every numerical assertion expression are
unchanged, verified by AST comparison. The exception assertion changes from
`pytest.raises(ValueError, match="do not cover")` to the equivalent
`self.assertRaisesRegex(ValueError, "do not cover")`.

The ordinary project environment lacks DBSTEP/Morfeus. It discovers all seven
tests and reports seven explicit reference-environment skips. The frozen
reference environment discovers and **passes all seven with zero skips**.
Reference-environment execution remains mandatory for scientific validation;
ordinary-environment skips are not reported as scientific passes.

[The receipt](receipt.json) binds both executions, the current test source and
its original source. The original test source is retained byte-for-byte here as
[historical_test_source.txt](historical_test_source.txt), in commit `156f0cd0`,
and in the already bound convergence campaign snapshot. All scientific study
controller, helper, proof, engine and reference source bindings remain unchanged.
The accepted numerical campaign was not rerun or overwritten for this test
discovery correction.

Historical scientific admission must verify the old test binding against its
immutable campaign snapshot, then separately verify this current-test receipt
and required reference-environment result. Every other original source binding
continues to be checked unchanged. Full project and reference-suite discovery
are rerun separately by the scientific replay reviewer.
