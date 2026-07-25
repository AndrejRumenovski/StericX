#!/usr/bin/env python3
"""Benchmark StericX Sterimol parameters against morfeus-fsu.

The script reads every XYZ structure in ``data/xyz``, obtains the attachment
axis from ``data/reactions_raw.csv`` (or the XYZ comment as a fallback), and
calculates reference L, B1, and B5 values with ``morfeus.Sterimol``. It then
invokes the release StericX CLI once over the same ordered structures, decodes
the resulting 64-byte records, reports regressions, and writes correlation
plots under ``docs/``.

Dependencies:
    morfeus-fsu (import name: morfeus), numpy, matplotlib, scipy, pandas

Run:
    python validate_stericx.py
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from morfeus import Sterimol
    from scipy.stats import linregress
except ImportError as exc:  # pragma: no cover - depends on the host environment.
    missing = getattr(exc, "name", "required package")
    raise SystemExit(
        f"Missing dependency `{missing}`. Install morfeus-fsu, numpy, "
        "matplotlib, scipy, and pandas before running validate_stericx.py."
    ) from exc


RECORD_DTYPE = np.dtype(
    [
        ("l", "=f4"),
        ("b1", "=f4"),
        ("b5", "=f4"),
        ("nbo_charge", "=f4"),
        ("ir_freq", "=f4"),
        ("temp_k", "=f4"),
        ("exp_ddg", "=f4"),
        ("reserved", "=f4", (9,)),
    ]
)
EXPECTED_RECORD_BYTES = 64
AXIS_PATTERN = re.compile(r"attach_idx=(\d+)\s+axis_idx=(\d+)")
VDW_RADII_ANGSTROM = {
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
}


@dataclass(frozen=True)
class AxisMetadata:
    """Attachment-axis and physical-organic values for one XYZ file."""

    reaction_id: str
    attach_idx: int
    neighbor_idx: int
    nbo_charge: float = 0.0
    ir_freq: float = 0.0
    temp_k: float = 298.15
    exp_ddg: float = 0.0


@dataclass(frozen=True)
class RegressionMetrics:
    """Linear comparison statistics for one Sterimol parameter."""

    parameter: str
    sample_count: int
    r_squared: float
    slope: float
    intercept: float
    rmse: float


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Compare StericX Sterimol descriptors with morfeus-fsu."
    )
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        default=project_root / "data" / "xyz",
        help="Directory containing XYZ files.",
    )
    parser.add_argument(
        "--reactions-csv",
        type=Path,
        default=project_root / "data" / "reactions_raw.csv",
        help="CSV containing attachment and primary-vector indices.",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=project_root / "target" / "release" / "stericx",
        help="Compiled StericX release binary.",
    )
    parser.add_argument(
        "--docs-dir",
        type=Path,
        default=project_root / "docs",
        help="Destination directory for correlation plots.",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="Do not build the release binary when it is missing.",
    )
    return parser.parse_args()


def read_xyz(path: Path) -> tuple[list[str], np.ndarray, str]:
    """Read elements, coordinates, and comment from one standard XYZ file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError("file has fewer than two XYZ header lines")
    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"invalid atom count `{lines[0].strip()}`") from exc
    if atom_count < 2:
        raise ValueError("at least two atoms are required for Sterimol")
    if len(lines) < atom_count + 2:
        raise ValueError(
            f"declares {atom_count} atoms but contains only {len(lines) - 2} atom lines"
        )

    elements: list[str] = []
    coordinates: list[list[float]] = []
    for atom_offset, line in enumerate(lines[2 : atom_count + 2], start=1):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"atom {atom_offset} has fewer than four fields")
        try:
            xyz = [float(value) for value in fields[1:4]]
        except ValueError as exc:
            raise ValueError(f"atom {atom_offset} has invalid coordinates") from exc
        if not np.isfinite(xyz).all():
            raise ValueError(f"atom {atom_offset} has non-finite coordinates")
        elements.append(fields[0])
        coordinates.append(xyz)
    return elements, np.asarray(coordinates, dtype=float), lines[1]


def load_axis_metadata(csv_path: Path) -> dict[str, AxisMetadata]:
    """Load CSV metadata keyed by XYZ basename."""
    if not csv_path.is_file():
        print(
            f"[warning] Metadata CSV not found: {csv_path}; "
            "falling back to XYZ comments."
        )
        return {}

    frame = pd.read_csv(csv_path)
    required = {
        "Reaction_ID",
        "Ligand_XYZ_Path",
        "Attach_Atom_Idx",
        "Primary_Bond_Vector_Idx",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        print(
            "[warning] Metadata CSV lacks "
            + ", ".join(missing)
            + "; falling back to XYZ comments."
        )
        return {}

    metadata: dict[str, AxisMetadata] = {}
    for row_number, row in frame.iterrows():
        filename = Path(str(row["Ligand_XYZ_Path"])).name
        if filename in metadata:
            raise ValueError(
                f"metadata CSV contains duplicate XYZ basename `{filename}`"
            )
        metadata[filename] = AxisMetadata(
            reaction_id=str(row["Reaction_ID"]),
            attach_idx=int(row["Attach_Atom_Idx"]),
            neighbor_idx=int(row["Primary_Bond_Vector_Idx"]),
            nbo_charge=float(row.get("NBO_Charge", 0.0)),
            ir_freq=float(row.get("IR_Frequency", 0.0)),
            temp_k=float(row.get("Temp_K", 298.15)),
            exp_ddg=float(row.get("Exp_ddG_kcal_mol", 0.0)),
        )
        values = metadata[filename]
        if not np.isfinite(
            [values.nbo_charge, values.ir_freq, values.temp_k, values.exp_ddg]
        ).all():
            raise ValueError(
                f"metadata row {row_number + 2} contains non-finite values"
            )
    print(f"[info] Loaded attachment metadata for {len(metadata)} structures.")
    return metadata


def metadata_from_comment(path: Path, comment: str) -> AxisMetadata:
    """Extract zero-based StericX axis indices from the XYZ comment."""
    match = AXIS_PATTERN.search(comment)
    if match is None:
        raise ValueError(
            "attachment indices are absent from both reactions CSV and XYZ comment"
        )
    return AxisMetadata(
        reaction_id=path.stem,
        attach_idx=int(match.group(1)),
        neighbor_idx=int(match.group(2)),
    )


def morfeus_sterimol(
    elements: list[str],
    coordinates: np.ndarray,
    metadata: AxisMetadata,
) -> tuple[float, float, float]:
    """Calculate a Morfeus reference using StericX-compatible atomic radii."""
    atom_count = len(elements)
    if not 0 <= metadata.attach_idx < atom_count:
        raise IndexError(f"attachment index {metadata.attach_idx} is out of bounds")
    if not 0 <= metadata.neighbor_idx < atom_count:
        raise IndexError(f"neighbor index {metadata.neighbor_idx} is out of bounds")
    if metadata.attach_idx == metadata.neighbor_idx:
        raise ValueError("attachment and neighbor indices are identical")

    radii = np.asarray(
        [VDW_RADII_ANGSTROM.get(element.upper(), 1.80) for element in elements],
        dtype=float,
    )
    # Morfeus is 1-indexed and names the axis origin "dummy". It excludes that
    # atom from surface projection, matching the conventional Sterimol setup.
    reference = Sterimol(
        elements,
        coordinates.copy(),
        dummy_index=metadata.attach_idx + 1,
        attached_index=metadata.neighbor_idx + 1,
        radii=radii,
        n_rot_vectors=3_600,
    )
    return (
        # StericX implements the raw geometric extent max(z + r). Morfeus also
        # exposes that value before applying its historical +0.40 Å correction.
        float(reference.L_value_uncorrected),
        float(reference.B_1_value),
        float(reference.B_5_value),
    )


def ensure_release_binary(binary: Path, no_build: bool) -> None:
    """Build the release binary when needed and allowed."""
    if binary.is_file():
        return
    if no_build:
        raise FileNotFoundError(f"release binary not found: {binary}")

    project_root = Path(__file__).resolve().parent
    print("[info] Release binary is missing; running `cargo build --release`.")
    completed = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=project_root,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        details = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"release build failed:\n{details}")
    if not binary.is_file():
        raise FileNotFoundError(f"release build did not produce {binary}")


def run_stericx(
    binary: Path,
    xyz_dir: Path,
    validation_rows: pd.DataFrame,
    work_dir: Path,
) -> np.ndarray:
    """Run the Rust parser and decode its flat 64-byte output records."""
    csv_path = work_dir / "validation_reactions.csv"
    sigpack_path = work_dir / "validation.sigpack"
    validation_rows.to_csv(csv_path, index=False)
    command = [
        str(binary),
        "parse",
        "--csv",
        str(csv_path),
        "--xyz-dir",
        str(xyz_dir),
        "--output",
        str(sigpack_path),
    ]
    print("[info] Executing:", " ".join(command))
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.stdout.strip():
        print("[stericx stdout]")
        print(completed.stdout.rstrip())
    if completed.returncode != 0:
        details = completed.stderr.strip() or "no stderr output"
        raise RuntimeError(
            f"StericX parse exited with status {completed.returncode}:\n{details}"
        )
    if not sigpack_path.is_file():
        raise FileNotFoundError("StericX completed without producing a sigpack file")
    byte_count = sigpack_path.stat().st_size
    if byte_count % EXPECTED_RECORD_BYTES != 0:
        raise ValueError(
            f"sigpack size {byte_count} is not divisible by {EXPECTED_RECORD_BYTES}"
        )

    records = np.fromfile(sigpack_path, dtype=RECORD_DTYPE)
    if len(records) != len(validation_rows):
        raise ValueError(
            f"StericX returned {len(records)} records for "
            f"{len(validation_rows)} validation rows"
        )
    return records


def regression_metrics(
    parameter: str,
    reference: np.ndarray,
    stericx: np.ndarray,
) -> RegressionMetrics:
    """Calculate direct RMSE and least-squares trendline statistics."""
    if len(reference) < 2:
        raise ValueError(f"{parameter} regression requires at least two samples")
    if np.ptp(reference) <= np.finfo(float).eps:
        slope = float("nan")
        intercept = float("nan")
        r_squared = float("nan")
    else:
        fit = linregress(reference, stericx)
        slope = float(fit.slope)
        intercept = float(fit.intercept)
        r_squared = float(fit.rvalue**2)
    rmse = float(np.sqrt(np.mean(np.square(stericx - reference))))
    return RegressionMetrics(
        parameter=parameter,
        sample_count=len(reference),
        r_squared=r_squared,
        slope=slope,
        intercept=intercept,
        rmse=rmse,
    )


def plot_correlation(
    parameter: str,
    reference: np.ndarray,
    stericx: np.ndarray,
    metrics: RegressionMetrics,
    output_path: Path,
) -> None:
    """Write a publication-quality correlation plot."""
    figure, axis = plt.subplots(figsize=(6.5, 6.0), constrained_layout=True)
    axis.scatter(
        reference,
        stericx,
        s=48,
        alpha=0.82,
        color="#176B87",
        edgecolors="white",
        linewidths=0.7,
        label=f"Structures (n={len(reference)})",
        zorder=3,
    )

    lower = float(min(reference.min(), stericx.min()))
    upper = float(max(reference.max(), stericx.max()))
    padding = max((upper - lower) * 0.08, 0.05)
    x_values = np.linspace(lower - padding, upper + padding, 200)
    axis.plot(
        x_values,
        x_values,
        linestyle="--",
        linewidth=1.2,
        color="#7A7A7A",
        label="Identity",
        zorder=1,
    )
    if np.isfinite(metrics.slope):
        axis.plot(
            x_values,
            metrics.slope * x_values + metrics.intercept,
            linewidth=2.0,
            color="#D95F02",
            label=(
                f"Fit: y={metrics.slope:.3f}x{metrics.intercept:+.3f}\n"
                f"$R^2$={metrics.r_squared:.4f}, RMSE={metrics.rmse:.3f} Å"
            ),
            zorder=2,
        )

    axis.set_xlim(x_values[0], x_values[-1])
    axis.set_ylim(x_values[0], x_values[-1])
    axis.set_xlabel(f"Morfeus {parameter} (Å)", fontsize=12)
    axis.set_ylabel(f"StericX {parameter} (Å)", fontsize=12)
    axis.set_title(f"Sterimol {parameter}: StericX vs morfeus-fsu", fontsize=13)
    axis.grid(alpha=0.2, linewidth=0.7)
    axis.legend(frameon=True, fontsize=9, loc="best")
    figure.savefig(output_path, dpi=400, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    args = parse_args()
    try:
        if not args.xyz_dir.is_dir():
            raise FileNotFoundError(f"XYZ directory not found: {args.xyz_dir}")
        xyz_files = sorted(args.xyz_dir.glob("*.xyz"))
        if not xyz_files:
            raise FileNotFoundError(f"no .xyz files found in {args.xyz_dir}")
        ensure_release_binary(args.binary, args.no_build)
        metadata_by_filename = load_axis_metadata(args.reactions_csv)

        reference_rows: list[dict[str, object]] = []
        validation_rows: list[dict[str, object]] = []
        failures: list[tuple[str, str]] = []
        print(f"[info] Evaluating {len(xyz_files)} XYZ files with Morfeus.")
        for index, xyz_path in enumerate(xyz_files, start=1):
            try:
                elements, coordinates, comment = read_xyz(xyz_path)
                metadata = metadata_by_filename.get(xyz_path.name)
                if metadata is None:
                    metadata = metadata_from_comment(xyz_path, comment)
                l_value, b1_value, b5_value = morfeus_sterimol(
                    elements,
                    coordinates,
                    metadata,
                )
                reference_rows.append(
                    {
                        "Reaction_ID": metadata.reaction_id,
                        "xyz_filename": xyz_path.name,
                        "morfeus_l": l_value,
                        "morfeus_b1": b1_value,
                        "morfeus_b5": b5_value,
                    }
                )
                validation_rows.append(
                    {
                        "Reaction_ID": metadata.reaction_id,
                        "Ligand_XYZ_Path": xyz_path.name,
                        "Attach_Atom_Idx": metadata.attach_idx,
                        "Primary_Bond_Vector_Idx": metadata.neighbor_idx,
                        "NBO_Charge": metadata.nbo_charge,
                        "IR_Frequency": metadata.ir_freq,
                        "Temp_K": metadata.temp_k,
                        "Exp_ddG_kcal_mol": metadata.exp_ddg,
                    }
                )
                print(
                    f"[{index:04d}/{len(xyz_files):04d}] {xyz_path.name}: "
                    f"L={l_value:.4f}, B1={b1_value:.4f}, B5={b5_value:.4f}"
                )
            except Exception as exc:
                failures.append((xyz_path.name, str(exc)))
                print(
                    f"[{index:04d}/{len(xyz_files):04d}] "
                    f"{xyz_path.name}: skipped - {exc}",
                    file=sys.stderr,
                )

        if len(reference_rows) < 2:
            raise RuntimeError(
                "fewer than two structures produced valid Morfeus references"
            )
        validation_frame = pd.DataFrame(validation_rows)
        with tempfile.TemporaryDirectory(prefix="stericx_validation_") as temporary:
            records = run_stericx(
                args.binary,
                args.xyz_dir,
                validation_frame,
                Path(temporary),
            )

        reference_frame = pd.DataFrame(reference_rows)
        comparisons = {
            "L": (
                reference_frame["morfeus_l"].to_numpy(dtype=float),
                records["l"].astype(float),
                "sterimol_l_corr.png",
            ),
            "B1": (
                reference_frame["morfeus_b1"].to_numpy(dtype=float),
                records["b1"].astype(float),
                "sterimol_b1_corr.png",
            ),
            "B5": (
                reference_frame["morfeus_b5"].to_numpy(dtype=float),
                records["b5"].astype(float),
                "sterimol_b5_corr.png",
            ),
        }
        args.docs_dir.mkdir(parents=True, exist_ok=True)
        metrics: list[RegressionMetrics] = []
        for parameter, (reference, stericx, filename) in comparisons.items():
            result = regression_metrics(parameter, reference, stericx)
            metrics.append(result)
            output_path = args.docs_dir / filename
            plot_correlation(parameter, reference, stericx, result, output_path)
            print(f"[info] Wrote {output_path}")

        summary = pd.DataFrame(
            [
                {
                    "Parameter": metric.parameter,
                    "N": metric.sample_count,
                    "R^2": metric.r_squared,
                    "Slope": metric.slope,
                    "Intercept": metric.intercept,
                    "RMSE_A": metric.rmse,
                }
                for metric in metrics
            ]
        )
        print("\nSterimol validation summary")
        print(summary.to_string(index=False, float_format=lambda value: f"{value:.6f}"))
        print(f"\nCompleted: {len(reference_rows)} compared, {len(failures)} skipped.")
        if failures:
            print("Skipped files:")
            for filename, reason in failures:
                print(f"  - {filename}: {reason}")
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"[fatal] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
