# Reproducible chemical evidence bundles

These deterministic gzip/JSONL bundles preserve exact retained file bytes, original
workspace-relative paths, and SHA256 checksums. `manifest.json` hashes the bundles.
Restore/verify without overwriting existing evidence:

```sh
.venv/bin/python docs/experimental_validation/archive_chemical_data.py restore
```

- `published_inputs`: 135 original MolSSI Kraken SDF exports, library identities,
  author CSV, original preparation audit and source-P audit. The absent 906 record
  is retained as an HTTP failure. Source coordinates are never newly optimized.
- `generated_inputs`: the first ETKDG/MMFF experiment, including rejected candidates.
  Its incomplete axial checks are a disclosed preparation limitation.
- `corrected_generated_inputs`: the separate MDL-aromaticity preparation, all 218
  generated candidates before/after optimization, rejection reasons, 69 retained
  states, actual MMFF energies, and authoritative native observations.
- `native_results`: original published and first-generated native observations plus
  independent Morfeus/DBSTEP and statistical checks.
- `analysis_and_provenance`: both analysis versions, source/API receipts, missing
  energy responses, correction metadata/table comparison, and conformer sensitivity.

The SI PDFs remain in the local evidence directory and can be retrieved using their
recorded primary URLs. They are not bundled here. Numerical facts and original
source attributions are retained; third-party files are not relicensed under
StericX's software license. The corrected SI metadata specifies CC BY-NC 4.0;
MolSSI/API and source-repository material has its own provenance and terms, with no
blanket redistribution license inferred. This local research archive is not a release.

All molecules are free ligands. No source conformer energies were recovered; no
published Boltzmann population is reconstructed. Generated force-field energies
belong only to the explicitly separate preparation experiment. Experimental sample
size is eleven ligands, never the number of conformers or the 1,566-row library.
