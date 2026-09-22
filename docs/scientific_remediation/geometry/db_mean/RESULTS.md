# Finite database means and standardized search distances

The retained [before database receipt](before_manifest.json) reproduces an
infinite group mean from four identical finite SDF observations. Accumulating
in f64 and checking the resulting f32 descriptor fixes this arithmetic failure.
The group mean is exactly the individual observation, approximately
1.9999593712379312e38 Å. Four focused tests cover ordinary means, large finite
values, empty input and nonfinite input ([test log](tests.log)).

The associated [search witness](search_after_manifest.json) exposed the same
intermediate-overflow problem in library standardization. Before correction,
identical large inputs and varying large inputs both produced successful JSON
with null distances. Standardization now uses f64 mean/variance intermediates
and rejects nonfinite distances before successful output. The population
standard deviation definition, f32 output and existing constant-feature cutoff
are unchanged. Eight focused tests include ordinary and large exact one-sigma
distances, constant large inputs, nonfinite descriptors and output overflow
([test log](search_tests.log)).

The final v4 native commands are retained in [final_v4/receipt.json](final_v4/receipt.json).
The database succeeds with the exact finite mean. Identical search inputs fail
with the existing all-constant-feature error. For a library containing
2^126 and 3·2^126, the analytic population standard deviation is 2^126; the
two distances exactly match the expected f32 values 0.6490590572357178 and
1.3509409427642822, in ascending order.

That receipt retains two initial checker failures: it parsed Rust's shortest
round-trip f32 JSON decimals as Python f64 and compared them directly to the
full f32 values. The [supplemental comparison](final_v4/serialized_f32_comparison.json)
decodes those same unchanged observations at the declared f32 boundary and
requires exact equality, with no tolerance. Both helper sources and all input,
binary, output and original receipt identities are preserved.

Repository Clippy and formatting checks passed after these source changes
([Clippy log](clippy.log), [format log](format.log)). A bounded scan of other
public aggregation callers found no additional demonstrated finite-sum overflow:
remaining matching regression sums were test code, model arithmetic was f64,
and the relevant reaction/descriptor means were handled by their owners. This
is a scoped review, not a proof of global numerical or experimental validity.
