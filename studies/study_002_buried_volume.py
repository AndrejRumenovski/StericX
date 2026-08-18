#!/usr/bin/env python3
"""Validate native coordination-aware buried volumes and rerun Ni-hDA.

The reference path reproduces the public Kraken/Morfeus protocol on the same
RDKit conformers used by StericX. It then validates the version-two binary
matrix, compares the approximate geometries with the official Kraken DFT
descriptor, and evaluates a preregistered one-feature Ni-hDA model.

Dependencies:
    morfeus, numpy, matplotlib, scipy, pandas

Run:
    uv run --extra science python studies/study_002_buried_volume.py
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import struct
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from morfeus import BuriedVolume
    from scipy.stats import linregress
except ImportError as exc:  # pragma: no cover - environment dependent
    raise SystemExit(
        "Missing study dependency. Run with "
        "`uv run --extra science python studies/study_002_buried_volume.py`."
    ) from exc


KRAKEN_COMMIT: Final[str] = "7b5f182fdc77334b713729a1f99ae25eaedbce69"
KRAKEN_SOURCE_URL: Final[str] = (
    "https://github.com/SigmanGroup/kraken/blob/"
    f"{KRAKEN_COMMIT}/kraken/morfeus_properties.py"
)
KRAKEN_SOURCE_SHA256: Final[str] = (
    "93d3ad5486e226dd4b49a8953797e471b19862e60b8da617d8afdc9377e8ca27"
)
KRAKEN_MORFEUS_VERSION: Final[str] = "0.7.2"
VALIDATION_MORFEUS_VERSION: Final[str] = importlib.metadata.version("morfeus-ml")
TRAIN_IDS: Final[tuple[int, ...]] = (
    401,
    498,
    724,
    785,
    1057,
    1058,
    2062,
    2063,
    2064,
    2067,
)
BLIND_ID: Final[int] = 723
OFFICIAL_FEATURE: Final[str] = "vbur_max_delta_qvbur_min"
SPHERE_RADIUS: Final[float] = 3.5
DENSITY: Final[float] = 0.01
CENTER_DISTANCE: Final[float] = 2.1
RADII_SCALE: Final[float] = 1.17
SIGPACK_V2_MAGIC: Final[bytes] = b"SIGPKV2\0"
SIGPACK_V2_RECORD_BYTES: Final[int] = 128
VDW_RADII: Final[dict[str, float]] = {
    "H": 1.20,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "SI": 2.10,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "NI": 1.63,
}
BV_FIELDS: Final[tuple[str, ...]] = (
    "vbur_boltz",
    "vbur_min",
    "vbur_max",
    "vbur_delta",
    "qvbur_min_boltz",
    "qvbur_max_boltz",
    "max_delta_qvbur_boltz",
    "max_delta_qvbur_min",
    "max_delta_qvbur_max",
    "max_delta_qvbur_delta",
    "max_delta_qvbur_vburminconf",
    "near_vbur_boltz",
    "far_vbur_boltz",
    "conformer_count",
    "sphere_radius",
    "grid_density",
)
CORE_PARITY_FIELDS: Final[tuple[str, ...]] = (
    "vbur",
    "qvbur_min",
    "qvbur_max",
    "max_delta_qvbur",
    "near_vbur",
    "far_vbur",
)


@dataclass(frozen=True)
class ComparisonMetrics:
    """Direct parity and linear-regression statistics."""

    parameter: str
    count: int
    r2: float | None
    slope: float | None
    intercept: float | None
    rmse: float
    mae: float
    mean_relative_error_percent: float | None
    max_absolute_error: float


@dataclass(frozen=True)
class ModelMetrics:
    """Regression metrics for one fixed partition."""

    count: int
    r2: float | None
    mae: float
    rmse: float


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description="Run StericX Study 002 buried-volume fidelity validation."
    )
    parser.add_argument(
        "--reactions-csv",
        type=Path,
        default=root / "data" / "reactions_raw.csv",
    )
    parser.add_argument("--xyz-dir", type=Path, default=root / "data")
    parser.add_argument(
        "--catalog",
        type=Path,
        default=root / "data" / "official" / "ni_hda_kraken.csv",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=root / "target" / "release" / "stericx",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "docs" / "study_002",
    )
    parser.add_argument(
        "--sigpack-output",
        type=Path,
        default=root / "data" / "reactions_v2.sigpack",
    )
    parser.add_argument(
        "--native-conformers-output",
        type=Path,
        default=root / "data" / "buried_volume_conformers.csv",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="Require an existing release binary.",
    )
    return parser.parse_args(argv)


def atomic_write_text(path: Path, content: str) -> None:
    """Atomically replace one UTF-8 artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    """Atomically replace one CSV artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False, float_format="%.10g")
    os.replace(temporary, path)


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    """Read one standard XYZ frame."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path} has no XYZ header")
    atom_count = int(lines[0].strip())
    if len(lines) < atom_count + 2:
        raise ValueError(f"{path} is truncated")
    elements: list[str] = []
    coordinates: list[list[float]] = []
    for line in lines[2 : atom_count + 2]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"{path} contains an incomplete atom line")
        elements.append(fields[0])
        coordinates.append([float(value) for value in fields[1:4]])
    array = np.asarray(coordinates, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError(f"{path} contains non-finite coordinates")
    return elements, array


def donor_neighbors(
    elements: list[str],
    coordinates: np.ndarray,
    donor_idx: int,
    reference_idx: int,
) -> list[int]:
    """Match StericX's deterministic nearest-heavy-atom donor topology."""
    candidates = sorted(
        (
            (float(np.sum((coordinates[index] - coordinates[donor_idx]) ** 2)), index)
            for index, element in enumerate(elements)
            if index != donor_idx and element.upper() != "H"
        ),
        key=lambda item: (item[0], item[1]),
    )
    neighbors = [reference_idx]
    neighbors.extend(index for _, index in candidates if index != reference_idx)
    neighbors = neighbors[:3]
    if len(neighbors) != 3:
        raise ValueError("three heavy donor substituents are required")
    return neighbors


def inferred_virtual_center(
    elements: list[str],
    coordinates: np.ndarray,
    donor_idx: int,
    neighbors: list[int],
) -> np.ndarray:
    """Infer the lone-pair direction used when xTB LMO centres are unavailable."""
    vectors = np.asarray(
        [
            (coordinates[index] - coordinates[donor_idx])
            / np.linalg.norm(coordinates[index] - coordinates[donor_idx])
            for index in neighbors
        ]
    )
    direction = -vectors.sum(axis=0)
    norm = float(np.linalg.norm(direction))
    if norm > 0.01:
        direction /= norm
    else:
        direction = np.cross(vectors[0], vectors[1])
        direction /= np.linalg.norm(direction)
        positive = minimum_heavy_clearance(
            elements,
            coordinates,
            donor_idx,
            coordinates[donor_idx] + direction,
        )
        negative = minimum_heavy_clearance(
            elements,
            coordinates,
            donor_idx,
            coordinates[donor_idx] - direction,
        )
        if negative > positive:
            direction *= -1.0
    return coordinates[donor_idx] + CENTER_DISTANCE * direction


def minimum_heavy_clearance(
    elements: list[str],
    coordinates: np.ndarray,
    donor_idx: int,
    point: np.ndarray,
) -> float:
    """Smallest squared distance from a candidate centre to a heavy ligand atom."""
    return min(
        float(np.sum((coordinate - point) ** 2))
        for index, (element, coordinate) in enumerate(
            zip(elements, coordinates, strict=True)
        )
        if index != donor_idx and element.upper() != "H"
    )


def morfeus_reference(
    elements: list[str],
    coordinates: np.ndarray,
    donor_idx: int,
    reference_idx: int,
) -> dict[str, float]:
    """Reproduce Kraken's three-orientation Morfeus calculation."""
    neighbors = donor_neighbors(elements, coordinates, donor_idx, reference_idx)
    center = inferred_virtual_center(elements, coordinates, donor_idx, neighbors)
    extended_elements = [*elements, "Ni"]
    extended_coordinates = np.vstack([coordinates, center])
    center_index = len(extended_elements)
    radii = np.asarray(
        [
            VDW_RADII.get(element.upper(), 1.80) * RADII_SCALE
            for element in extended_elements
        ]
    )

    qvbur_all: list[float] = []
    ovbur_all: list[float] = []
    max_delta_all: list[float] = []
    first: dict[str, float] | None = None
    for neighbor_idx in neighbors:
        buried = BuriedVolume(
            extended_elements,
            extended_coordinates.copy(),
            center_index,
            excluded_atoms=[center_index],
            radii=radii,
            include_hs=False,
            radius=SPHERE_RADIUS,
            density=DENSITY,
            z_axis_atoms=[donor_idx + 1],
            xz_plane_atoms=[neighbor_idx + 1],
        )
        buried.octant_analysis()
        quadrants = np.asarray(
            list(buried.quadrants["buried_volume"].values()),
            dtype=float,
        )
        octants = np.asarray(
            list(buried.octants["buried_volume"].values()),
            dtype=float,
        )
        qvbur_all.extend(quadrants.tolist())
        ovbur_all.extend(octants.tolist())
        max_delta_all.append(
            max(abs(quadrants[index] - quadrants[index - 1]) for index in range(4))
        )
        if first is None:
            first = {
                "vbur": float(buried.buried_volume),
                "percent_vbur": float(buried.fraction_buried_volume * 100.0),
                "near_vbur": float(octants[4:].sum()),
                "far_vbur": float(octants[:4].sum()),
            }
    assert first is not None
    return {
        **first,
        "qvbur_min": min(qvbur_all),
        "qvbur_max": max(qvbur_all),
        "max_delta_qvbur": max(max_delta_all),
        "ovbur_min": min(ovbur_all),
        "ovbur_max": max(ovbur_all),
    }


def build_reference(reactions: pd.DataFrame, xyz_dir: Path) -> pd.DataFrame:
    """Calculate reference descriptors for every retained conformer."""
    rows: list[dict[str, object]] = []
    total = int(reactions["Conformer_XYZ_Paths"].str.split(";").map(len).sum())
    completed = 0
    for _, reaction in reactions.iterrows():
        paths = str(reaction["Conformer_XYZ_Paths"]).split(";")
        weights = [
            float(value)
            for value in str(reaction["Conformer_Boltzmann_Weights"]).split(";")
        ]
        for conformer_index, (relative_path, weight) in enumerate(
            zip(paths, weights, strict=True)
        ):
            elements, coordinates = read_xyz(xyz_dir / relative_path)
            values = morfeus_reference(
                elements,
                coordinates,
                int(reaction["Attach_Atom_Idx"]),
                int(reaction["Primary_Bond_Vector_Idx"]),
            )
            rows.append(
                {
                    "Reaction_ID": str(reaction["Reaction_ID"]),
                    "Conformer_Index": conformer_index,
                    "Conformer_XYZ_Path": relative_path,
                    "Boltzmann_Weight": weight,
                    **values,
                }
            )
            completed += 1
            print(f"[reference {completed:03d}/{total:03d}] {relative_path}")
    return pd.DataFrame.from_records(rows)


def ensure_binary(binary: Path, no_build: bool) -> None:
    """Build the optimized Rust executable when necessary."""
    if binary.is_file() and no_build:
        return
    if no_build:
        raise FileNotFoundError(f"release binary not found: {binary}")
    root = Path(__file__).resolve().parent.parent
    completed = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())
    if not binary.is_file():
        raise FileNotFoundError(f"release build did not produce {binary}")


def run_native(
    binary: Path,
    reactions_csv: Path,
    xyz_dir: Path,
    sigpack_output: Path,
    conformer_output: Path,
) -> str:
    """Execute the native descriptor pipeline and return its console log."""
    sigpack_output.parent.mkdir(parents=True, exist_ok=True)
    conformer_output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(binary),
        "buried-volume",
        "--csv",
        str(reactions_csv),
        "--xyz-dir",
        str(xyz_dir),
        "--output",
        str(sigpack_output),
        "--per-conformer-output",
        str(conformer_output),
        "--sphere-radius",
        str(SPHERE_RADIUS),
        "--density",
        str(DENSITY),
        "--center-distance",
        str(CENTER_DISTANCE),
        "--radii-scale",
        str(RADII_SCALE),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or "native command failed")
    print(completed.stdout.rstrip())
    return completed.stdout


def decode_sigpack_v2(path: Path, reaction_ids: list[str]) -> pd.DataFrame:
    """Decode and validate the version-two file emitted by Rust."""
    raw = path.read_bytes()
    if len(raw) < 64 or raw[:8] != SIGPACK_V2_MAGIC:
        raise ValueError("invalid sigpack v2 magic or truncated header")
    version, endian = struct.unpack_from("=II", raw, 8)
    record_count = struct.unpack_from("=Q", raw, 16)[0]
    record_size, descriptor_count = struct.unpack_from("=II", raw, 24)
    if (version, endian, record_size, descriptor_count) != (
        2,
        0x01020304,
        SIGPACK_V2_RECORD_BYTES,
        32,
    ):
        raise ValueError("sigpack v2 header does not match the frozen schema")
    if record_count != len(reaction_ids):
        raise ValueError("sigpack v2 record count does not match reactions CSV")
    expected_bytes = 64 + record_count * record_size
    if len(raw) != expected_bytes:
        raise ValueError(
            f"sigpack v2 declares {expected_bytes} bytes, found {len(raw)}"
        )
    dtype = np.dtype(
        [
            ("reaction", "=f4", (16,)),
            ("buried_volume", "=f4", (16,)),
        ]
    )
    records = np.frombuffer(raw, dtype=dtype, count=record_count, offset=64)
    rows = []
    for reaction_id, record in zip(reaction_ids, records, strict=True):
        values = {
            field: float(record["buried_volume"][index])
            for index, field in enumerate(BV_FIELDS)
        }
        rows.append({"Reaction_ID": reaction_id, **values})
    return pd.DataFrame.from_records(rows)


def aggregate_reference(reference: pd.DataFrame) -> pd.DataFrame:
    """Build the exact v2 ensemble columns from reference conformer rows."""
    rows: list[dict[str, object]] = []
    for reaction_id, group in reference.groupby("Reaction_ID", sort=False):
        weights = group["Boltzmann_Weight"].to_numpy(dtype=float, copy=True)
        weights /= weights.sum()
        vbur = group["vbur"].to_numpy(dtype=float)
        delta = group["max_delta_qvbur"].to_numpy(dtype=float)
        minimum_vbur_index = int(np.argmin(vbur))
        rows.append(
            {
                "Reaction_ID": reaction_id,
                "vbur_boltz": float(np.dot(vbur, weights)),
                "vbur_min": float(vbur.min()),
                "vbur_max": float(vbur.max()),
                "vbur_delta": float(np.ptp(vbur)),
                "qvbur_min_boltz": float(
                    np.dot(group["qvbur_min"].to_numpy(dtype=float), weights)
                ),
                "qvbur_max_boltz": float(
                    np.dot(group["qvbur_max"].to_numpy(dtype=float), weights)
                ),
                "max_delta_qvbur_boltz": float(np.dot(delta, weights)),
                "max_delta_qvbur_min": float(delta.min()),
                "max_delta_qvbur_max": float(delta.max()),
                "max_delta_qvbur_delta": float(np.ptp(delta)),
                "max_delta_qvbur_vburminconf": float(delta[minimum_vbur_index]),
                "near_vbur_boltz": float(
                    np.dot(group["near_vbur"].to_numpy(dtype=float), weights)
                ),
                "far_vbur_boltz": float(
                    np.dot(group["far_vbur"].to_numpy(dtype=float), weights)
                ),
                "conformer_count": float(len(group)),
                "sphere_radius": SPHERE_RADIUS,
                "grid_density": DENSITY,
            }
        )
    return pd.DataFrame.from_records(rows)


def comparison_metrics(
    parameter: str,
    reference: np.ndarray,
    candidate: np.ndarray,
) -> ComparisonMetrics:
    """Calculate direct and trendline parity statistics."""
    residual = candidate - reference
    varying = float(np.ptp(reference)) > np.finfo(float).eps
    regression = linregress(reference, candidate) if varying else None
    nonzero = np.abs(reference) > 1.0
    relative = (
        float(np.mean(np.abs(residual[nonzero] / reference[nonzero])) * 100.0)
        if nonzero.any()
        else None
    )
    return ComparisonMetrics(
        parameter=parameter,
        count=int(reference.size),
        r2=float(regression.rvalue**2) if regression else None,
        slope=float(regression.slope) if regression else None,
        intercept=float(regression.intercept) if regression else None,
        rmse=float(np.sqrt(np.mean(residual**2))),
        mae=float(np.mean(np.abs(residual))),
        mean_relative_error_percent=relative,
        max_absolute_error=float(np.max(np.abs(residual))),
    )


def model_metrics(actual: np.ndarray, predicted: np.ndarray) -> ModelMetrics:
    """Calculate fixed-partition model metrics."""
    residual = predicted - actual
    total = float(np.sum((actual - actual.mean()) ** 2))
    r2 = (
        float(1.0 - np.sum(residual**2) / total)
        if actual.size > 1 and total > np.finfo(float).eps
        else None
    )
    return ModelMetrics(
        count=int(actual.size),
        r2=r2,
        mae=float(np.mean(np.abs(residual))),
        rmse=float(np.sqrt(np.mean(residual**2))),
    )


def fit_line(x_values: np.ndarray, y_values: np.ndarray) -> tuple[float, float]:
    """Fit y = intercept + slope*x by ordinary least squares."""
    slope, intercept = np.polyfit(x_values, y_values, 1)
    return float(intercept), float(slope)


def fixed_feature_loo(x_values: np.ndarray, y_values: np.ndarray) -> np.ndarray:
    """Leave each row out without changing the preregistered descriptor."""
    predictions = np.empty_like(y_values)
    for held_out in range(y_values.size):
        mask = np.arange(y_values.size) != held_out
        intercept, slope = fit_line(x_values[mask], y_values[mask])
        predictions[held_out] = intercept + slope * x_values[held_out]
    return predictions


def plot_parity(
    output: Path,
    reference: np.ndarray,
    candidate: np.ndarray,
    label: str,
    metrics: ComparisonMetrics,
    reference_name: str = "Morfeus reference",
) -> None:
    """Save one high-resolution descriptor parity plot."""
    figure, axis = plt.subplots(figsize=(6.2, 5.8))
    axis.scatter(reference, candidate, color="#176B87", s=42, alpha=0.82)
    limits = [
        min(float(reference.min()), float(candidate.min())),
        max(float(reference.max()), float(candidate.max())),
    ]
    padding = max((limits[1] - limits[0]) * 0.06, 0.05)
    limits = [limits[0] - padding, limits[1] + padding]
    axis.plot(limits, limits, "--", color="#333333", linewidth=1.1)
    axis.set(xlim=limits, ylim=limits)
    axis.set_xlabel(f"{reference_name} {label} (Å³)")
    axis.set_ylabel(f"StericX {label} (Å³)")
    r2 = "unavailable" if metrics.r2 is None else f"{metrics.r2:.6f}"
    axis.set_title(f"{label}: native parity ($R^2$={r2})")
    figure.tight_layout()
    figure.savefig(output, dpi=400)
    plt.close(figure)


def write_descriptor_spec(output_dir: Path) -> None:
    """Freeze the recovered public Kraken convention and approximation boundary."""
    text = f"""# StericX Buried-Volume Descriptor Specification

## Frozen reference

- Kraken source commit: `{KRAKEN_COMMIT}`
- Reference implementation: [morfeus_properties.py]({KRAKEN_SOURCE_URL})
- Reference source SHA-256: `{KRAKEN_SOURCE_SHA256}`
- Kraken-pinned Morfeus version: {KRAKEN_MORFEUS_VERSION}
- Validation runtime Morfeus version: {VALIDATION_MORFEUS_VERSION}
- Integration sphere radius: {SPHERE_RADIUS:.1f} Å
- Filled-grid density: {DENSITY:.2f} Å³ per point
- Bondi radius scale: {RADII_SCALE:.2f}
- Hydrogens: excluded
- Virtual metal distance from phosphorus: {CENTER_DISTANCE:.1f} Å

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
"""
    atomic_write_text(output_dir / "DESCRIPTOR_SPEC.md", text)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        reactions = pd.read_csv(args.reactions_csv)
        required = {
            "Reaction_ID",
            "Source_ID",
            "Conformer_XYZ_Paths",
            "Conformer_Boltzmann_Weights",
            "Attach_Atom_Idx",
            "Primary_Bond_Vector_Idx",
        }
        missing = sorted(required.difference(reactions.columns))
        if missing:
            raise ValueError(f"reactions CSV lacks: {', '.join(missing)}")

        reference = build_reference(reactions, args.xyz_dir)
        atomic_write_csv(args.output_dir / "morfeus_conformer_reference.csv", reference)
        ensure_binary(args.binary, args.no_build)
        native_log = run_native(
            args.binary,
            args.reactions_csv,
            args.xyz_dir,
            args.sigpack_output,
            args.native_conformers_output,
        )
        native = pd.read_csv(args.native_conformers_output)
        merged = reference.merge(
            native,
            on=["Reaction_ID", "Conformer_Index", "Conformer_XYZ_Path"],
            suffixes=("_reference", "_native"),
            validate="one_to_one",
        )
        parity: list[ComparisonMetrics] = []
        for field in CORE_PARITY_FIELDS:
            metrics = comparison_metrics(
                field,
                merged[f"{field}_reference"].to_numpy(dtype=float),
                merged[f"{field}_native"].to_numpy(dtype=float),
            )
            parity.append(metrics)
            plot_parity(
                args.output_dir / f"{field}_parity.png",
                merged[f"{field}_reference"].to_numpy(dtype=float),
                merged[f"{field}_native"].to_numpy(dtype=float),
                field,
                metrics,
            )

        reference_ensemble = aggregate_reference(reference)
        native_ensemble = decode_sigpack_v2(
            args.sigpack_output,
            reactions["Reaction_ID"].astype(str).tolist(),
        )
        ensemble_merged = reference_ensemble.merge(
            native_ensemble,
            on="Reaction_ID",
            suffixes=("_reference", "_native"),
            validate="one_to_one",
        )
        ensemble_parity = [
            comparison_metrics(
                field,
                ensemble_merged[f"{field}_reference"].to_numpy(dtype=float),
                ensemble_merged[f"{field}_native"].to_numpy(dtype=float),
            )
            for field in BV_FIELDS
        ]
        atomic_write_csv(args.output_dir / "ensemble_reference.csv", reference_ensemble)
        atomic_write_csv(args.output_dir / "native_v2_descriptors.csv", native_ensemble)

        catalog_features = pd.read_csv(
            args.catalog,
            index_col=0,
            usecols=lambda column: (
                column == OFFICIAL_FEATURE or column.startswith("Unnamed:")
            ),
        )
        catalog_features.index = pd.to_numeric(catalog_features.index).astype(int)
        source_map = reactions[["Reaction_ID", "Source_ID"]].copy()
        source_map["Source_ID"] = pd.to_numeric(source_map["Source_ID"]).astype(int)
        official_comparison = native_ensemble.merge(
            source_map,
            on="Reaction_ID",
            validate="one_to_one",
        )
        official_comparison["official_kraken_max_delta_qvbur_min"] = [
            float(catalog_features.at[source_id, OFFICIAL_FEATURE])
            for source_id in official_comparison["Source_ID"]
        ]
        official_metrics = comparison_metrics(
            "max_delta_qvbur_min_official",
            official_comparison["official_kraken_max_delta_qvbur_min"].to_numpy(
                dtype=float
            ),
            official_comparison["max_delta_qvbur_min"].to_numpy(dtype=float),
        )
        atomic_write_csv(
            args.output_dir / "official_kraken_comparison.csv",
            official_comparison,
        )
        plot_parity(
            args.output_dir / "official_kraken_descriptor_comparison.png",
            official_comparison["official_kraken_max_delta_qvbur_min"].to_numpy(
                dtype=float
            ),
            official_comparison["max_delta_qvbur_min"].to_numpy(dtype=float),
            "max_delta_qvbur_min",
            official_metrics,
            reference_name="Official Kraken",
        )

        training_target_table = pd.read_csv(
            args.catalog,
            index_col=0,
            usecols=lambda column: column == "ddG_abs" or column.startswith("Unnamed:"),
        )
        training_target_table.index = pd.to_numeric(training_target_table.index).astype(
            int
        )
        feature_by_source = official_comparison.set_index("Source_ID")[
            "max_delta_qvbur_min"
        ]
        x_train = feature_by_source.loc[list(TRAIN_IDS)].to_numpy(dtype=float)
        y_train = training_target_table.loc[list(TRAIN_IDS), "ddG_abs"].to_numpy(
            dtype=float
        )
        del training_target_table
        intercept, slope = fit_line(x_train, y_train)
        training_prediction = intercept + slope * x_train
        loo_prediction = fixed_feature_loo(x_train, y_train)
        training_metrics = model_metrics(y_train, training_prediction)
        loo_metrics = model_metrics(y_train, loo_prediction)
        model_artifact = {
            "schema_version": 2,
            "model": "preregistered_univariate_ols",
            "feature": "stericx_max_delta_qvbur_min",
            "intercept": intercept,
            "slope": slope,
            "training_source_ids": list(TRAIN_IDS),
            "blind_source_ids": [BLIND_ID],
            "training_feature_minimum": float(x_train.min()),
            "training_feature_maximum": float(x_train.max()),
            "sigpack_v2_sha256": hashlib.sha256(
                args.sigpack_output.read_bytes()
            ).hexdigest(),
            "kraken_source_commit": KRAKEN_COMMIT,
            "virtual_center_method": "opposed normalized donor-substituent vectors",
            "created_at_utc": datetime.now(UTC).isoformat(),
        }
        atomic_write_text(
            args.output_dir / "native_model.json",
            json.dumps(model_artifact, indent=2, sort_keys=True) + "\n",
        )

        blind_feature = float(feature_by_source.at[BLIND_ID])
        blind_prediction = intercept + slope * blind_feature
        frozen = pd.DataFrame(
            [
                {
                    "Source_ID": BLIND_ID,
                    "Feature": "stericx_max_delta_qvbur_min",
                    "Feature_Value": blind_feature,
                    "Predicted_ddG_kcal_mol": blind_prediction,
                    "Target_Accessed_During_Fit": False,
                }
            ]
        )
        frozen_path = args.output_dir / "frozen_predictions.csv"
        atomic_write_csv(frozen_path, frozen)
        frozen_sha256 = hashlib.sha256(frozen_path.read_bytes()).hexdigest()
        revealed_target_table = pd.read_csv(
            args.catalog,
            index_col=0,
            usecols=lambda column: column == "ddG_abs" or column.startswith("Unnamed:"),
        )
        revealed_target_table.index = pd.to_numeric(revealed_target_table.index).astype(
            int
        )
        blind_target = float(revealed_target_table.at[BLIND_ID, "ddG_abs"])
        blind_error = blind_prediction - blind_target
        scored = frozen.assign(
            Experimental_ddG_kcal_mol=blind_target,
            Residual_kcal_mol=blind_error,
            Absolute_Error_kcal_mol=abs(blind_error),
        )
        atomic_write_csv(args.output_dir / "scored_blind_predictions.csv", scored)

        figure, axis = plt.subplots(figsize=(6.2, 5.8))
        axis.scatter(y_train, loo_prediction, color="#176B87", s=60, label="LOO train")
        axis.scatter(
            [blind_target],
            [blind_prediction],
            color="#D95F02",
            marker="*",
            s=180,
            label="Historical blind 723",
        )
        limits = [min(y_train.min(), loo_prediction.min()) - 0.1, 2.2]
        axis.plot(limits, limits, "--", color="#333333")
        axis.set(xlim=limits, ylim=limits)
        axis.set_xlabel(r"Experimental $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
        axis.set_ylabel(r"Predicted $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
        axis.legend(frameon=False)
        axis.set_title("Study 002: Native Buried-Volume Model")
        figure.tight_layout()
        figure.savefig(
            args.output_dir / "ni_hda_native_buried_volume_parity.png",
            dpi=400,
        )
        plt.close(figure)

        core_errors = [
            metric.mean_relative_error_percent
            for metric in parity
            if metric.mean_relative_error_percent is not None
        ]
        geometry_gate = max(core_errors, default=math.inf) < 1.0
        official_gate = (official_metrics.r2 or 0.0) > 0.99
        model_gate = (loo_metrics.r2 or -math.inf) >= 0.752
        holdout_gate = abs(blind_error) <= 0.373
        results = {
            "schema_version": 2,
            "generated_at_utc": datetime.now(UTC).isoformat(),
            "kraken_source_commit": KRAKEN_COMMIT,
            "kraken_source_url": KRAKEN_SOURCE_URL,
            "kraken_source_sha256": KRAKEN_SOURCE_SHA256,
            "kraken_morfeus_version": KRAKEN_MORFEUS_VERSION,
            "validation_morfeus_version": VALIDATION_MORFEUS_VERSION,
            "geometry": {
                "sphere_radius_angstrom": SPHERE_RADIUS,
                "density_angstrom3_per_point": DENSITY,
                "center_distance_angstrom": CENTER_DISTANCE,
                "radii_scale": RADII_SCALE,
                "include_hydrogens": False,
                "official_center_method": "xTB localized molecular orbital centre",
                "stericx_center_method": "opposed normalized donor-substituent vectors",
            },
            "records": {
                "reactions": len(reactions),
                "conformers": len(reference),
                "sigpack_v2_bytes": args.sigpack_output.stat().st_size,
            },
            "per_conformer_parity": [asdict(metric) for metric in parity],
            "ensemble_v2_parity": [asdict(metric) for metric in ensemble_parity],
            "official_kraken_comparison": asdict(official_metrics),
            "native_model": {
                "feature": "stericx_max_delta_qvbur_min",
                "intercept": intercept,
                "slope": slope,
                "training": asdict(training_metrics),
                "fixed_feature_loo": asdict(loo_metrics),
                "historical_blind": {
                    "source_id": BLIND_ID,
                    "feature_value": blind_feature,
                    "predicted_ddg_kcal_mol": blind_prediction,
                    "experimental_ddg_kcal_mol": blind_target,
                    "absolute_error_kcal_mol": abs(blind_error),
                },
                "frozen_prediction_sha256": frozen_sha256,
            },
            "success_gates": {
                "morfeus_mean_relative_error_below_1_percent": geometry_gate,
                "official_kraken_descriptor_r2_above_0_99": official_gate,
                "native_fixed_feature_loo_q2_at_least_0_752": model_gate,
                "historical_blind_error_at_most_0_373_kcal_mol": holdout_gate,
            },
            "native_console_log": native_log,
        }
        atomic_write_text(
            args.output_dir / "study_results.json",
            json.dumps(results, indent=2, sort_keys=True) + "\n",
        )
        write_descriptor_spec(args.output_dir)
        write_report(args.output_dir, results)

        print("\nStericX Study 002 complete")
        print(
            "  Morfeus parity gate: "
            f"{'PASS' if geometry_gate else 'FAIL'} "
            f"(worst mean relative error={max(core_errors):.6f}%)"
        )
        print(
            "  Official Kraken descriptor: "
            f"R²={official_metrics.r2:.4f}, RMSE={official_metrics.rmse:.4f} Å³"
        )
        print(
            "  Native Ni-hDA model: "
            f"train R²={training_metrics.r2:.4f}, "
            f"LOO Q²={loo_metrics.r2:.4f}"
        )
        print(
            f"  Historical blind 723: predicted={blind_prediction:.4f}, "
            f"experimental={blind_target:.4f}, error={abs(blind_error):.4f} kcal/mol"
        )
        print(f"  Report: {args.output_dir / 'STUDY_002.md'}")
        return 0
    except (
        FileNotFoundError,
        KeyError,
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
        ValueError,
    ) as exc:
        print(f"Study 002 failed: {exc}", file=sys.stderr)
        return 1


def write_report(output_dir: Path, results: dict[str, object]) -> None:
    """Write the measured Study 002 model card without hiding failed gates."""
    native = results["native_model"]
    official = results["official_kraken_comparison"]
    gates = results["success_gates"]
    parity = results["per_conformer_parity"]
    records = results["records"]
    assert isinstance(native, dict)
    assert isinstance(official, dict)
    assert isinstance(gates, dict)
    assert isinstance(parity, list)
    assert isinstance(records, dict)
    training = native["training"]
    loo = native["fixed_feature_loo"]
    blind = native["historical_blind"]
    assert isinstance(training, dict)
    assert isinstance(loo, dict)
    assert isinstance(blind, dict)
    parity_rows = "\n".join(
        f"| {metric['parameter']} | {metric['r2']:.6f} | "
        f"{metric['rmse']:.6g} | {metric['mean_relative_error_percent']:.6f}% |"
        for metric in parity
        if isinstance(metric, dict)
        and metric["r2"] is not None
        and metric["mean_relative_error_percent"] is not None
    )
    gate_rows = "\n".join(
        f"| `{name}` | {'PASS' if passed else 'FAIL'} |"
        for name, passed in gates.items()
    )
    report = f"""# StericX Study 002

## Coordination-aware buried volume

This study implements the public Kraken three-orientation quadrant protocol and
validates the native Rust voxel engine against Morfeus on exactly the same
geometries. The complete convention and its approximation boundary are frozen
in [DESCRIPTOR_SPEC.md](DESCRIPTOR_SPEC.md).

The version-two matrix contains {records["reactions"]} reaction records and
aggregates {records["conformers"]} retained conformers. Version one remains
readable and unchanged.

| Per-conformer descriptor | R² vs Morfeus | RMSE (Å³) | Mean relative error |
|---|---:|---:|---:|
{parity_rows}

![Total buried-volume parity](vbur_parity.png)

![Quadrant-anisotropy parity](max_delta_qvbur_parity.png)

## Official Kraken comparison

The native geometry uses ETKDGv3/MMFF94 conformers and an inferred lone-pair
direction. Official Kraken uses CREST/xTB/DFT ensembles and xTB localized
molecular-orbital centres. This comparison therefore tests the complete
approximate workflow, not the Rust voxel arithmetic alone.

| Quantity | Value |
|---|---:|
| R² against official `{OFFICIAL_FEATURE}` | {official["r2"]:.4f} |
| RMSE | {official["rmse"]:.4f} Å³ |
| Slope | {official["slope"]:.4f} |
| Intercept | {official["intercept"]:.4f} Å³ |

![Official Kraken comparison](official_kraken_descriptor_comparison.png)

## Locked Ni-hDA rerun

The model uses only the native ensemble minimum of `max_delta_qvbur`. The ten
published training IDs and historical blind ligand 723 are unchanged.

| Quantity | Value |
|---|---:|
| Training R² | {training["r2"]:.4f} |
| Training RMSE | {training["rmse"]:.4f} kcal/mol |
| Fixed-feature LOO Q² | {loo["r2"]:.4f} |
| Fixed-feature LOO RMSE | {loo["rmse"]:.4f} kcal/mol |
| Blind prediction | {blind["predicted_ddg_kcal_mol"]:.4f} kcal/mol |
| Blind experimental | {blind["experimental_ddg_kcal_mol"]:.4f} kcal/mol |
| Blind absolute error | {blind["absolute_error_kcal_mol"]:.4f} kcal/mol |

![Native model parity](ni_hda_native_buried_volume_parity.png)

## Preregistered success gates

| Gate | Result |
|---|---|
{gate_rows}

Failed gates are retained as scientific results. Passing the same-geometry
Morfeus gate establishes implementation fidelity; it does not imply that
approximate RDKit geometries reproduce the official quantum-chemical
descriptor.

## Next experimental boundary

This remains a historical reproduction. A prospective claim requires a frozen
ranked ligand deck followed by new measurements performed without refitting to
those outcomes.
"""
    atomic_write_text(output_dir / "STUDY_002.md", report)


if __name__ == "__main__":
    raise SystemExit(main())
