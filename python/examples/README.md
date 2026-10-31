# Research workflows

`stericx_research.ipynb` is a complete, executed notebook. It covers single
geometries, ordered conformers, arithmetic and explicitly justified energy
weighting, sensitivity, CSV/JSON/NumPy/pandas output, RDKit, batch errors,
hard-sphere radius profiles and independently checked density-grid profiles.

From the repository root, install the development package and optional tools:

```bash
python -m pip install './python[notebook]'
python python/examples/research_workflow.py --output workflow-output
cd python/examples
jupyter lab stericx_research.ipynb
```

Use a fresh output directory for the script. Run the notebook with this directory
as its working directory; generated exports go into ignored `notebook_output/`.
The package is not published on PyPI. See the distribution evidence in
[`DISTRIBUTION.md`](../../docs/paton_recommendations/python/DISTRIBUTION.md).

The six ligand-401 geometries and their recorded MMFF94 energies come from the
repository's retained conformer study. Energy weighting is conditional on these
states, unit multiplicity and the stated temperature; it is not a historical
Kraken population or a measured free-energy distribution. The mapped SDF is the
existing topology fixture. `data/manifest.json` records each source path and
SHA-256 hash. The Ne cube is from the pinned official DBSTEP repository; its
MIT license is included. Its quantum-chemical method is unspecified in the
fixture, and the tutorial does not invent one.

Indices are zero-based; lengths are Å, volumes Å³, and buried fractions are
reported in percent. The scan varies the sphere radius at a fixed center.
Density occupation uses an explicit threshold in electrons/bohr³ and a separate
grid-node volume definition. The notebook deliberately retains a missing batch
file and invalid/out-of-domain scan positions as visible errors.

The script and all 11 notebook code cells were executed against the tested
Linux release wheel. Validation logs and the artifact receipt are in
[`docs/paton_recommendations/python/`](../../docs/paton_recommendations/python/).
