# StericX Buried-Volume Descriptor Specification

## Frozen reference

- Kraken source commit: `7b5f182fdc77334b713729a1f99ae25eaedbce69`
- Reference implementation: [morfeus_properties.py](https://github.com/SigmanGroup/kraken/blob/7b5f182fdc77334b713729a1f99ae25eaedbce69/kraken/morfeus_properties.py)
- Reference source SHA-256: `93d3ad5486e226dd4b49a8953797e471b19862e60b8da617d8afdc9377e8ca27`
- Kraken-pinned Morfeus version: 0.7.2
- Validation runtime Morfeus version: 0.8.0
- Integration sphere radius: 3.5 Å
- Filled-grid density: 0.01 Å³ per point
- Bondi radius scale: 1.17
- Hydrogens: excluded
- Virtual metal distance from phosphorus: 2.1 Å

For each conformer, phosphorus defines the Z axis. Each of the three nearest
heavy phosphorus substituents defines the XZ plane once. `qvbur_min` and
`qvbur_max` are extrema over all twelve resulting quadrant volumes.
`max_delta_qvbur` is the largest absolute difference between cyclically
adjacent quadrants over those three orientations.

Across conformers, StericX records the Boltzmann average, minimum, maximum,
range, and the property value from the conformer having minimum total buried
volume. The version-two binary schema stores the stable 64-byte v1 reaction
record followed by one 64-byte buried-volume descriptor block.

## Approximation boundary

Official free-ligand Kraken calculations choose a phosphorus lone-pair
direction from xTB localized-molecular-orbital centres. Plain XYZ files do not
contain those centres. StericX therefore places its virtual centre opposite the
sum of the three normalized P-substituent vectors, with a deterministic
maximum-clearance normal for planar geometries. Reference parity tests use that
same centre so they isolate the Rust geometry implementation. Comparison with
official Kraken values separately measures the combined effect of approximate
centres, RDKit/MMFF conformers, and the absence of CREST/xTB/DFT geometries.
