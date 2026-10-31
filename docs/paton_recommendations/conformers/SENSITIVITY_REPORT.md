# Conformer sensitivity on the retained MMFF dataset

The Rust ensemble API now produces sensitivity reports with each descriptor's
minimum, maximum, range, population standard deviation, arithmetic mean, median,
and explicitly justified weighted mean. Every exact tie at an extreme retains
its conformer index, original title and source. A report contains settings,
units, population provenance, missing-state policy and failed-state identities.
There is no combined flexibility score.

## Recorded results

The existing `data/reactions_raw.csv` supplies **11 ligand ensembles and 56
retained conformers**, with their recorded MMFF94 relative energies in kcal/mol.
All 154 scalar summary rows are available. The
[complete CSV](sensitivity_v1/sensitivity.csv) reports every ligand and all 14
scalar fields, including extrema identities. The
[manifest](sensitivity_v1/manifest.json) binds every source geometry, source CSV,
engine, command and per-ligand JSON/CSV report. No conformer was removed.

The following B5 values are in angstrom; rounding here is display-only. The
complete files preserve the full native and reduced precision.

| Ligand | Minimum | Maximum | Range | Standard deviation | MMFF-energy weighted mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| 401 | 6.59563 | 8.87194 | 2.27631 | 1.11823 | 7.08955 |
| 498 | 7.18615 | 7.18618 | 0.00003 | 0.00001 | 7.18616 |
| 723 | 7.85021 | 7.85021 | 0.00000 | 0.00000 | 7.85021 |
| 724 | 8.84945 | 8.84955 | 0.00010 | 0.00005 | 8.84950 |
| 785 | 9.46573 | 9.46579 | 0.00006 | 0.00003 | 9.46576 |
| 1057 | 7.85538 | 7.85538 | 0.00000 | 0.00000 | 7.85538 |
| 1058 | 7.84199 | 7.85529 | 0.01329 | 0.00665 | 7.85011 |
| 2062 | 9.49813 | 10.26925 | 0.77112 | 0.31320 | 10.26715 |
| 2063 | 9.44971 | 9.46918 | 0.01947 | 0.00948 | 9.46914 |
| 2064 | 9.44002 | 9.49534 | 0.05533 | 0.02185 | 9.47839 |
| 2067 | 8.75823 | 9.50603 | 0.74780 | 0.20980 | 9.33917 |

Descriptor sensitivities differ: ligand 401 spans 2.27631 Å in B5 and 2.06699 Å³
in buried volume; ligand 2067 spans 0.74780 Å and 20.74364 Å³ respectively. These
are properties of the retained geometries and stated descriptor settings. A
small observed range does not establish molecular rigidity or search coverage;
ligand 1057 has only one retained conformer.

## Conventions and limits

This report uses the current descriptor CLI defaults: coordination-axis
Sterimol with a +0.40 Å L correction, inferred center distance 2.28 Å, sphere
radius 3.5 Å, volume density 0.01 Å³, radii scale 1.17, and hydrogen exclusion
from hard-sphere volume. The former study's lower-level 2.1 Å center convention
is distinct; its frozen results have not been overwritten or redefined.

Weights are recalculated from the **printed existing MMFF energies** at the
source row's temperature. They approximate retained-state energy weighting, not
experimental equilibrium populations. Unit multiplicity is explicit. Distinct
thermodynamic state membership, symmetry, missing conformers and entropy have
not been independently established. These are not the unavailable historical
Kraken conformer energies or populations. The complete source XYZ comments and
energy metadata remain in each JSON report.

XYZ inputs have no supplied bond adjacency. Scalar atom correspondence is
checked, and per-conformer quadrant/octant frames remain available, but regional
ensemble statistics are explicitly unavailable without matching ordered
connectivity. The report does not silently average geometrically unrelated bins.
These generated descriptor values are not themselves an independent validation;
the next study compares reductions with independent references.

## Reproduction and API

```bash
cargo build --locked --example ensemble_report
.venv/bin/python scripts/report_conformer_sensitivity.py \
  --engine target/debug/examples/ensemble_report \
  --output docs/paton_recommendations/conformers/sensitivity_v1
```

Use a fresh output directory. The driver performs no descriptor chemistry: it
assembles unchanged source records, attaches explicit metadata, and calls Rust.
All failed commands remain in the manifest and their stderr files.

```rust
let descriptors = ensemble.descriptors(&options);
let statistics = descriptors.statistics(&Default::default(), population.as_ref())?;
let report = statistics.sensitivity("ligand identifier");
let csv_text = report.to_csv()?;
let json_text = serde_json::to_string_pretty(&report)?;
```

Archive JSON with the CSV projection to retain full settings and energy
provenance. Default incomplete-field statistics remain unavailable; explicitly
conditional reports retain the missing IDs and population mass. Seven statistics
tests cover closed-form values, failures, atom correspondence, tied extremes and
CSV identities containing commas, quotes and newlines. An initial test build
needed an explicit String type annotation in the new CSV formatter; its log is
retained alongside the passing run and clean Clippy check.
