"""Typed transport objects. Equations and scientific reductions live in Rust."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Generic, Literal, TypeVar

T = TypeVar("T")


@dataclass(frozen=True)
class VolumeOptions:
    sphere_radius: float = 3.5
    density: float = 0.01
    center_distance: float = 2.28
    radii_scale: float = 1.17
    include_hydrogens: bool = False


@dataclass(frozen=True)
class DescriptorOptions:
    donor_element: str = "P"
    donor_index: int | None = None
    reference_index: int | None = None
    axis: Literal["bond", "coordination"] = "bond"
    volume: VolumeOptions = field(default_factory=VolumeOptions)
    center: tuple[float, float, float] | None = None
    angular_samples: int | None = None
    l_correction: float | None = None


@dataclass(frozen=True)
class Energy:
    value: float
    unit: Literal["kcal_mol", "kj_mol", "hartree", "electron_volt"]
    kind: Literal["electronic", "force_field", "gibbs_free_energy"]
    model: str
    provenance: str
    temperature_k: float | None = None
    degeneracy: int = 1


@dataclass(frozen=True)
class Boltzmann:
    temperature_k: float
    assumptions: str
    kind: Literal["boltzmann"] = field(default="boltzmann", init=False)


@dataclass(frozen=True)
class SuppliedWeights:
    weights: tuple[float, ...]
    provenance: str
    kind: Literal["supplied"] = field(default="supplied", init=False)


@dataclass(frozen=True)
class Failure:
    code: str
    message: str


@dataclass(frozen=True)
class Outcome(Generic[T]):
    value: T | None
    error: Failure | None

    @property
    def successful(self) -> bool:
        return self.error is None

    def require(self) -> T:
        from ._native import CalculationError

        if self.error is not None:
            raise CalculationError(f"{self.error.code}: {self.error.message}")
        assert self.value is not None
        return self.value


@dataclass(frozen=True)
class Sterimol:
    l: float  # noqa: E741 — conventional Sterimol L field
    b1: float
    b5: float
    unit: str = "angstrom"


@dataclass(frozen=True)
class Pyramidalization:
    pyr_p: float
    pyr_alpha: float


@dataclass(frozen=True)
class VolumeSummary:
    buried_volume: float
    percent_buried_volume: float
    qvbur_min: float
    qvbur_max: float
    max_delta_qvbur: float
    ovbur_min: float
    ovbur_max: float
    near_vbur: float
    far_vbur: float


@dataclass(frozen=True)
class RegionalVolume:
    plane_atom_index: int
    axes: tuple[tuple[float, ...], ...]
    buried_volume: float
    quadrants: tuple[float, ...]
    octants: tuple[float, ...]
    near_vbur: float
    far_vbur: float


@dataclass(frozen=True)
class BuriedVolume:
    summary: VolumeSummary
    center: tuple[float, float, float]
    orientations: tuple[RegionalVolume, ...]
    volume_unit: str = "angstrom^3"
    percent_unit: str = "percent"


@dataclass(frozen=True)
class Geometry:
    elements: tuple[str, ...]
    coordinates_angstrom: tuple[tuple[float, ...], ...]
    radii_angstrom: tuple[float, ...]
    bonds: tuple[tuple[int, int], ...] | None


@dataclass(frozen=True)
class ConformerResult:
    index: int
    identifier: str
    source: str
    source_format: str
    source_geometry: str
    properties: tuple[dict[str, str], ...]
    geometry: Outcome[Geometry]
    energy: Energy | None
    energy_error: str | None
    sterimol: Outcome[Sterimol]
    buried_volume: Outcome[BuriedVolume]
    pyramidalization: Outcome[Pyramidalization]

    def descriptor_values(self) -> dict[str, float | None]:
        values: dict[str, float | None] = {}
        for key in ("l", "b1", "b5"):
            values[f"sterimol_{key}"] = (
                getattr(self.sterimol.value, key) if self.sterimol.successful else None
            )
        for key in VolumeSummary.__dataclass_fields__:
            values[key] = (
                getattr(self.buried_volume.value.summary, key)
                if self.buried_volume.successful
                else None
            )
        for key in ("pyr_p", "pyr_alpha"):
            values[key] = (
                getattr(self.pyramidalization.value, key)
                if self.pyramidalization.successful
                else None
            )
        return values

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ConformerIdentity:
    index: int
    identifier: str
    source: str


@dataclass(frozen=True)
class DescriptorStatistics:
    unit: str
    total_count: int
    successful_count: int
    missing_conformers: tuple[ConformerIdentity, ...]
    complete: bool
    unavailable_reason: str | None
    mean: float | None
    minimum: float | None
    maximum: float | None
    standard_deviation: float | None
    median: float | None
    weighted_mean: float | None
    boltzmann_mean: float | None
    retained_population_mass: float | None
    minimum_conformers: tuple[ConformerIdentity, ...]
    maximum_conformers: tuple[ConformerIdentity, ...]


@dataclass(frozen=True)
class Population:
    weights: tuple[float, ...]
    origin: dict[str, Any]


@dataclass(frozen=True)
class EnsembleResult:
    conformers: tuple[ConformerResult, ...]
    descriptor_options: dict[str, Any]
    statistics: dict[str, DescriptorStatistics]
    regions: dict[str, DescriptorStatistics]
    statistics_error: Failure | None
    population: Population | None
    missing_policy: str
    structural_correspondence: str | None
    regional_unavailable_reason: str | None
    sensitivity: dict[str, Any] | None
    sensitivity_csv: str | None
    # The unmodified complete Rust payload preserves all settings/provenance.
    _raw: dict[str, Any] = field(repr=False)

    def __getattr__(self, name: str) -> DescriptorStatistics:
        if name in self.statistics:
            return self.statistics[name]
        raise AttributeError(name)

    @property
    def sterimol_b5(self) -> DescriptorStatistics:
        return self.statistics["sterimol_b5"]

    @property
    def sterimol_b1(self) -> DescriptorStatistics:
        return self.statistics["sterimol_b1"]

    @property
    def sterimol_l(self) -> DescriptorStatistics:
        return self.statistics["sterimol_l"]

    def to_dict(self) -> dict[str, Any]:
        import copy

        return copy.deepcopy(self._raw)

    def to_records(self) -> list[dict[str, Any]]:
        return [
            {
                "conformer_index": c.index,
                "identifier": c.identifier,
                "source": c.source,
                **c.descriptor_values(),
                **{
                    f"{family}_error": getattr(c, family).error.message
                    if getattr(c, family).error
                    else None
                    for family in (
                        "geometry",
                        "sterimol",
                        "buried_volume",
                        "pyramidalization",
                    )
                },
                "energy_error": c.energy_error,
            }
            for c in self.conformers
        ]

    def to_pandas(self):
        """One row per original conformer, including failures; requires pandas."""
        import pandas as pd

        table = pd.DataFrame.from_records(self.to_records())
        table.attrs.update(
            units={k: s.unit for k, s in self.statistics.items()},
            descriptor_options=self.descriptor_options,
            missing_policy=self.missing_policy,
        )
        return table

    def to_numpy(self):
        """Structured array; failed numeric descriptors are NaN with error columns."""
        import numpy as np

        rows = self.to_records()
        numeric = (
            list(self.conformers[0].descriptor_values()) if self.conformers else []
        )
        keys = list(rows[0]) if rows else ["conformer_index"]
        dtype = [
            (
                key,
                "i8"
                if key == "conformer_index"
                else "f8"
                if key in numeric
                else "U"
                + str(
                    max(1, max((len(str(row[key] or "")) for row in rows), default=1))
                ),
            )
            for key in keys
        ]
        array = np.empty(len(rows), dtype=dtype)
        for i, row in enumerate(rows):
            for key in keys:
                array[key][i] = (
                    np.nan
                    if key in numeric and row[key] is None
                    else ""
                    if row[key] is None
                    else row[key]
                )
        return array


@dataclass(frozen=True)
class DensityVolume:
    method: str
    center_angstrom: tuple[float, float, float]
    radius_angstrom: float
    threshold_electrons_per_bohr3: float
    occupied_nodes: int
    sphere_nodes: int
    cell_volume_angstrom3: float
    buried_volume_angstrom3: float
    sampled_sphere_volume_angstrom3: float
    analytic_sphere_volume_angstrom3: float
    sphere_sampling_error_percent: float
    percent_buried_volume: float
    grid: dict[str, Any]


@dataclass(frozen=True)
class RadiusSample(Generic[T]):
    requested_radius_angstrom: float
    evaluated_radius_angstrom: float | None
    outcome: Outcome[T]


@dataclass(frozen=True)
class VolumeProfile(Generic[T]):
    definition: str
    center_angstrom: tuple[float, ...] | None
    samples: tuple[RadiusSample[T], ...]
    metadata: dict[str, Any]


@dataclass(frozen=True)
class BatchItem:
    path: str
    conformers: tuple[ConformerResult, ...]
    error: Failure | None
    descriptor_options: dict[str, Any] | None


def outcome(row: dict[str, Any], convert, code: str = "descriptor_error"):
    if row["status"] == "error":
        return Outcome(None, Failure(code, row["message"]))
    return Outcome(convert(row["value"]), None)


def volume(row):
    frames = tuple(
        RegionalVolume(
            **{
                **f,
                "axes": tuple(map(tuple, f["axes"])),
                "quadrants": tuple(f["quadrants"]),
                "octants": tuple(f["octants"]),
            }
        )
        for f in row["orientations"]
    )
    return BuriedVolume(VolumeSummary(**row["summary"]), tuple(row["center"]), frames)


def conformer(row):
    def geometry(g):
        return Geometry(
            tuple(g["elements"]),
            tuple(map(tuple, g["coordinates_angstrom"])),
            tuple(g["radii_angstrom"]),
            tuple(map(tuple, g["bonds"])) if g["bonds"] is not None else None,
        )

    return ConformerResult(
        **{
            **row,
            "properties": tuple(row["properties"]),
            "geometry": outcome(row["geometry"], geometry, "geometry_error"),
            "energy": Energy(**row["energy"]) if row["energy"] else None,
            "sterimol": outcome(row["sterimol"], lambda s: Sterimol(**s)),
            "buried_volume": outcome(row["buried_volume"], volume),
            "pyramidalization": outcome(
                row["pyramidalization"], lambda s: Pyramidalization(**s)
            ),
        }
    )


def ensemble_result(row):
    stat = row["statistics"]
    error = (
        Failure("ensemble_correspondence_error", stat["message"])
        if stat["status"] == "error"
        else None
    )
    data = stat.get("value", {})

    def statistics(s):
        return DescriptorStatistics(
            **{
                **s,
                **{
                    k: tuple(ConformerIdentity(**v) for v in s[k])
                    for k in (
                        "missing_conformers",
                        "minimum_conformers",
                        "maximum_conformers",
                    )
                },
            }
        )

    pop = data.get("population")
    return EnsembleResult(
        tuple(conformer(c) for c in row["descriptors"]["conformers"]),
        row["descriptors"]["options"],
        {k: statistics(s) for k, s in data.get("descriptors", {}).items()},
        {k: statistics(s) for k, s in data.get("regions", {}).items()},
        error,
        Population(tuple(pop["weights"]), pop["origin"]) if pop else None,
        data.get("options", {}).get("missing", "unavailable"),
        data.get("structural_correspondence"),
        data.get("regional_unavailable_reason"),
        row.get("sensitivity"),
        row.get("sensitivity_csv"),
        row,
    )


def density_volume(row):
    return DensityVolume(**{**row, "center_angstrom": tuple(row["center_angstrom"])})


def profile(row, convert):
    return VolumeProfile(
        row["definition"],
        tuple(row["center_angstrom"]) if row["center_angstrom"] else None,
        tuple(
            RadiusSample(
                s["requested_radius_angstrom"],
                s["evaluated_radius_angstrom"],
                outcome(s["outcome"], convert),
            )
            for s in row["samples"]
        ),
        {k: v for k, v in row.items() if k != "samples"},
    )
