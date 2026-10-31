"""Descriptor workflow using retained geometries and recorded MMFF energies."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import stericx as sx
from rdkit import Chem

DATA = Path(__file__).parent / "data"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    settings = json.loads((DATA / "ligand_401_options.json").read_text())
    settings["volume"] = sx.VolumeOptions(**settings.get("volume", {}))
    options = sx.DescriptorOptions(**settings)

    molecule = sx.Molecule.from_file(DATA / "ligand_401_conformer_0.xyz")
    single = molecule.descriptors(options)
    print("Sterimol (Å):", single.sterimol.require())
    print("Hard-sphere volume (Å³ and %):", single.buried_volume.require().summary)
    (args.output / "single.json").write_text(json.dumps(single.to_dict(), indent=2))

    ensemble = sx.Ensemble.from_file(DATA / "ligand_401_ensemble.xyz")
    unweighted = ensemble.descriptors(options)
    assert unweighted.sterimol_b5.boltzmann_mean is None
    energies = [
        sx.Energy(**row)
        for row in json.loads((DATA / "ligand_401_energies.json").read_text())
    ]
    ensemble.set_energies(energies)
    result = ensemble.descriptors(
        options,
        population=sx.Boltzmann(
            298.15,
            "Retained MMFF94 states at unit multiplicity; "
            "incomplete search and entropy not established.",
        ),
        ligand_id="SIG-NIHDA-401",
    )
    print("B5 ensemble (Å):", result.sterimol_b5)
    print("Conformer sensitivity:", result.sensitivity["ranges"])
    print("Regional ensemble limitation:", result.regional_unavailable_reason)
    result.to_pandas().to_csv(args.output / "conformers.csv", index=False)
    (args.output / "ensemble.json").write_text(json.dumps(result.to_dict(), indent=2))
    (args.output / "sensitivity.csv").write_text(result.sensitivity_csv)

    # SDF preserves the graph. This fixture has no assigned energy.
    mapped = sx.Ensemble.from_sdf(DATA / "mapped_ligand.sdf")
    mapped_options = sx.DescriptorOptions(axis="coordination", donor_index=1)
    mapped_result = mapped.descriptors(mapped_options)
    print("Mapped molecular-frame region fields:", len(mapped_result.regions))
    rdkit_molecule = Chem.MolFromMolFile(
        str(DATA / "mapped_ligand.sdf"), removeHs=False
    )
    rdkit_result = sx.Molecule.from_rdkit(rdkit_molecule).descriptors(mapped_options)
    assert (
        rdkit_result.sterimol.require()
        == mapped_result.conformers[0].sterimol.require()
    )

    batch = sx.calculate_files(
        [
            DATA / "ligand_401_conformer_0.xyz",
            DATA / "ligand_401_conformer_1.xyz",
            DATA / "missing.xyz",
        ],
        options,
        threads=2,
    )
    (args.output / "batch.json").write_text(
        json.dumps([asdict(item) for item in batch], indent=2)
    )
    assert batch[-1].error is not None  # Missing inputs remain records.
    profile = molecule.buried_volume_scan([2.5, 3.5, 4.5], options=options)
    (args.output / "hard_sphere_profile.json").write_text(
        json.dumps(asdict(profile), indent=2)
    )

    density = sx.DensityGrid.from_cube(
        DATA / "Ne_medium.cube",
        field="total_electron_density",
        density_unit="electrons_per_bohr_cubed",
        provenance="DBSTEP Ne total-SCF-density fixture; method unspecified",
    )
    density_profile = density.buried_volume_scan(
        [1.0, 1.5, 2.0],
        center_angstrom=(0.0, 0.0, 0.0),
        threshold_electrons_per_bohr3=0.0016,
    )
    (args.output / "density_profile.json").write_text(
        json.dumps(asdict(density_profile), indent=2)
    )
    print(
        "Density node counts:",
        [
            (s.outcome.require().occupied_nodes, s.outcome.require().sphere_nodes)
            for s in density_profile.samples
        ],
    )


if __name__ == "__main__":
    main()
