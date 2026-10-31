"""StericX: native steric descriptors, explicit ensembles, and density-grid volume.

Coordinates and radii use Å; atom indices are zero-based. Numerical chemistry
and all population/statistics calculations run in the shared Rust engine.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import asdict
from os import PathLike
from pathlib import Path

from . import _native
from ._native import CalculationError, DensityError, InputError, PopulationError
from ._types import (
    BatchItem,
    Boltzmann,
    BuriedVolume,
    ConformerResult,
    DensityVolume,
    DescriptorOptions,
    DescriptorStatistics,
    Energy,
    EnsembleResult,
    Failure,
    Geometry,
    Outcome,
    Population,
    Pyramidalization,
    RadiusSample,
    RegionalVolume,
    Sterimol,
    SuppliedWeights,
    VolumeOptions,
    VolumeProfile,
    conformer,
    density_volume,
    ensemble_result,
    profile,
    volume,
)

__version__ = "0.4.0.dev0"


def _json(value) -> str:
    try:
        return json.dumps(value, allow_nan=False)
    except (TypeError, ValueError) as error:
        raise InputError(str(error)) from error


def _options(options: DescriptorOptions | None) -> str:
    if options is not None and not isinstance(options, DescriptorOptions):
        raise InputError("options must be DescriptorOptions or None")
    return _json(asdict(options or DescriptorOptions()))


class Ensemble:
    """Ordered conformers. Missing energies never imply thermodynamic populations."""

    def __init__(self, engine):
        self._engine = engine

    @classmethod
    def from_file(cls, path: str | PathLike[str]) -> Ensemble:
        return cls(_native.EnsembleEngine.from_file(str(path)))

    @classmethod
    def from_sdf(cls, path: str | PathLike[str]) -> Ensemble:
        try:
            text = Path(path).read_text(encoding="utf-8")
        except OSError as error:
            raise InputError(str(error)) from error
        return cls.from_sdf_string(text, source=str(path))

    @classmethod
    def from_sdf_string(cls, text: str, *, source: str = "memory.sdf") -> Ensemble:
        return cls(_native.EnsembleEngine.from_text(text, "sdf", source))

    @classmethod
    def from_xyz_string(cls, text: str, *, source: str = "memory.xyz") -> Ensemble:
        return cls(_native.EnsembleEngine.from_text(text, "xyz", source))

    @classmethod
    def from_rdkit(cls, molecule) -> Ensemble:
        """Use existing RDKit coordinates and atom order without optimization."""
        conformers = list(molecule.GetConformers())
        if not conformers:
            raise InputError("RDKit molecule has no conformer coordinates")
        records = []
        for conf in conformers:
            identifier = (
                molecule.GetProp("_Name")
                if molecule.HasProp("_Name")
                else f"rdkit:{conf.GetId()}"
            )
            lines = [str(molecule.GetNumAtoms()), f"RDKit conformer {conf.GetId()}"]
            for atom in molecule.GetAtoms():
                point = conf.GetAtomPosition(atom.GetIdx())
                lines.append(
                    atom.GetSymbol()
                    + " "
                    + " ".join(format(float(x), ".17g") for x in point)
                )
            records.append(
                {
                    "xyz": "\n".join(lines) + "\n",
                    "identifier": identifier,
                    "conformer_id": conf.GetId(),
                }
            )
        bonds = [
            [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
            for bond in molecule.GetBonds()
        ]
        return cls(
            _native.EnsembleEngine.from_coordinate_records(_json(records), _json(bonds))
        )

    def __len__(self) -> int:
        return len(self._engine)

    @property
    def records(self) -> tuple[dict, ...]:
        return tuple(json.loads(self._engine.records()))

    def set_energies(self, energies: Sequence[Energy | None]) -> None:
        """Assign explicit energies transactionally in source order."""
        self._engine.set_energies(
            _json([asdict(e) if e is not None else None for e in energies])
        )

    def descriptors(
        self,
        options: DescriptorOptions | None = None,
        *,
        population: Boltzmann | SuppliedWeights | None = None,
        missing: str = "require_complete",
        ligand_id: str | None = None,
    ) -> EnsembleResult:
        if population is not None and not isinstance(
            population, (Boltzmann, SuppliedWeights)
        ):
            raise PopulationError(
                "population must be Boltzmann, SuppliedWeights or None"
            )
        raw = self._engine.analyze(
            _options(options),
            _json({"missing": missing}),
            _json(asdict(population) if population else None),
            ligand_id,
        )
        return ensemble_result(json.loads(raw))

    def buried_volume_scan(
        self,
        radii_angstrom: Sequence[float],
        *,
        conformer_index: int = 0,
        options: DescriptorOptions | None = None,
    ) -> VolumeProfile[BuriedVolume]:
        return profile(
            json.loads(
                self._engine.volume_scan(
                    conformer_index, list(radii_angstrom), _options(options)
                )
            ),
            volume,
        )


class Molecule:
    """One source conformer. Use Ensemble explicitly for multirecord input."""

    def __init__(self, ensemble: Ensemble):
        if len(ensemble) != 1:
            raise InputError(
                f"Molecule requires exactly one conformer, received {len(ensemble)}; "
                "use Ensemble"
            )
        self._ensemble = ensemble

    @classmethod
    def from_file(cls, path: str | PathLike[str]) -> Molecule:
        return cls(Ensemble.from_file(path))

    @classmethod
    def from_xyz_string(cls, text: str, *, source: str = "memory.xyz") -> Molecule:
        return cls(Ensemble.from_xyz_string(text, source=source))

    @classmethod
    def from_sdf_string(cls, text: str, *, source: str = "memory.sdf") -> Molecule:
        return cls(Ensemble.from_sdf_string(text, source=source))

    @classmethod
    def from_rdkit(cls, molecule, *, conformer_id: int | None = None) -> Molecule:
        from rdkit import Chem

        if conformer_id is None:
            return cls(Ensemble.from_rdkit(molecule))
        copy = Chem.Mol(molecule)
        selected = Chem.Conformer(molecule.GetConformer(conformer_id))
        copy.RemoveAllConformers()
        copy.AddConformer(selected, assignId=False)
        return cls(Ensemble.from_rdkit(copy))

    def descriptors(self, options: DescriptorOptions | None = None) -> ConformerResult:
        raw = json.loads(self._ensemble._engine.descriptors(_options(options)))
        return conformer(raw["conformers"][0])

    def buried_volume_scan(
        self,
        radii_angstrom: Sequence[float],
        *,
        options: DescriptorOptions | None = None,
    ) -> VolumeProfile[BuriedVolume]:
        return self._ensemble.buried_volume_scan(radii_angstrom, options=options)


class DensityGrid:
    """Declared total electron density; separate from hard-sphere molecular volume."""

    def __init__(self, engine):
        self._engine = engine

    @classmethod
    def from_cube(
        cls,
        path: str | PathLike[str],
        *,
        field: str,
        density_unit: str,
        provenance: str,
    ) -> DensityGrid:
        declaration = {"field": field, "unit": density_unit, "provenance": provenance}
        return cls(_native.DensityEngine.from_file(str(path), _json(declaration)))

    @property
    def metadata(self) -> dict:
        return json.loads(self._engine.metadata())

    def buried_volume(
        self,
        *,
        center_angstrom: Sequence[float],
        radius_angstrom: float,
        threshold_electrons_per_bohr3: float,
    ) -> DensityVolume:
        return density_volume(
            json.loads(
                self._engine.volume(
                    tuple(center_angstrom),
                    radius_angstrom,
                    threshold_electrons_per_bohr3,
                )
            )
        )

    def buried_volume_scan(
        self,
        radii_angstrom: Sequence[float],
        *,
        center_angstrom: Sequence[float],
        threshold_electrons_per_bohr3: float,
    ) -> VolumeProfile[DensityVolume]:
        return profile(
            json.loads(
                self._engine.scan(
                    tuple(center_angstrom),
                    list(radii_angstrom),
                    threshold_electrons_per_bohr3,
                )
            ),
            density_volume,
        )


def calculate_files(
    paths: Sequence[str | PathLike[str]],
    options: DescriptorOptions | None = None,
    *,
    threads: int | None = None,
) -> tuple[BatchItem, ...]:
    """Ordered native batch; unreadable files and failed conformers remain records."""
    rows = json.loads(
        _native.batch([str(p) for p in paths], _options(options), threads)
    )
    results = []
    for row in rows:
        result = row["result"]
        if result["status"] == "error":
            results.append(
                BatchItem(
                    row["path"], (), Failure("input_error", result["message"]), None
                )
            )
        else:
            data = result["value"]
            results.append(
                BatchItem(
                    row["path"],
                    tuple(conformer(c) for c in data["conformers"]),
                    None,
                    data["options"],
                )
            )
    return tuple(results)


__all__ = [
    "BatchItem",
    "Boltzmann",
    "BuriedVolume",
    "CalculationError",
    "ConformerResult",
    "DensityError",
    "DensityGrid",
    "DensityVolume",
    "DescriptorOptions",
    "DescriptorStatistics",
    "Energy",
    "Ensemble",
    "EnsembleResult",
    "Failure",
    "Geometry",
    "InputError",
    "Molecule",
    "Outcome",
    "Population",
    "PopulationError",
    "Pyramidalization",
    "RadiusSample",
    "RegionalVolume",
    "Sterimol",
    "SuppliedWeights",
    "VolumeOptions",
    "VolumeProfile",
    "calculate_files",
]
