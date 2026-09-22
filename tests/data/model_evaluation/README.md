These fixtures are byte-for-byte copies of the model evaluation failure
witnesses and their supporting synthetic dataset and fitted model from the
published scientific accuracy audit. Keeping them here lets the regression test
run from a fresh checkout without restoring the full evidence archive.

Source paths are relative to `docs/scientific_accuracy_audit/`. Sizes and SHA-256
hashes were verified against the audit's
[`publication/manifest.json`](../../../docs/scientific_accuracy_audit/publication/manifest.json).

| Fixture | Audit source | SHA-256 |
| --- | --- | --- |
| `evaluate_nan_predictions.csv` | `models/inputs/evaluate_nan_predictions.csv` | `dc329ab9d22a567a04ed2e3d60edff28d6e1c86fd02f8b183916b5c7550a18bd` |
| `evaluate_infinity_predictions.csv` | `models/inputs/evaluate_infinity_predictions.csv` | `e9ed836f09e6f5de6e599c8a71c04dcaae4fd3bd0ce0467e3bd52f000fc4acb0` |
| `evaluate_finite_mismatch_predictions.csv` | `models/inputs/evaluate_finite_mismatch_predictions.csv` | `49cf35860487cc81f836106b8a7ad31f61902bb94f6e4b7a888ae6d1731059e9` |
| `synthetic_linear.sigpack` | `models/inputs/synthetic_linear.sigpack` | `7330792a8662d7445a969fdf9759b757c6134bb0be4f5f9b4f3159fd0a6ad98d` |
| `synthetic_linear_labels.csv` | `models/inputs/synthetic_linear_labels.csv` | `c7c443db4a4d63f59847cf2a5ba2ac8aef44a0338235458bd803ee545fdc27db` |
| `report.json` | `models/raw/synthetic_linear/report.json` | `8f1aa7d7762a052dca2726fb25985255caba8f931efc687e78ea6c020933f6bf` |
