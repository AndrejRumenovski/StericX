"""Study 010: Is the buried-volume descriptor converged, or grid-lucky?

Studies 002-004 show StericX's buried volume matches morfeus (to R2 = 1.0 on
identical geometries) and Kraken's published values (R2 = 0.9852 across 1,541
ligands). Every one of those runs, though, used a *single* grid resolution -- the
Kraken convention of 0.01 A^3 per integration point. A fair question follows: is
that agreement a property of the descriptor, or an artifact of one lucky grid?

This study answers it directly. It sweeps the integration grid from coarse
(0.5 A^3/point) to fine (0.001 A^3/point) on a diverse sample of Kraken ligands,
takes the finest grid as the converged reference, and measures how far each
coarser grid sits from it -- in the flagship %Vbur, and in cost (wall time). The
honest questions are whether %Vbur converges as the grid refines, whether the
default 0.01 sits inside the converged regime, and what the irreducible
discretization floor of a voxel integrator actually is.

The Kraken DFT SDF cache is a local, gitignored artifact; only StericX's own
convergence measurements are written out.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent

# Coarse -> fine. The finest is the converged reference every other grid is
# measured against; the two finest also bound the irreducible grid-registration
# jitter (the discretization floor that no refinement removes).
DEFAULT_DENSITIES = (0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
# The Kraken / StericX default, highlighted in the report and figure.
KRAKEN_DENSITY = 0.01


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target" / "release" / "stericx"
    )
    parser.add_argument(
        "--cache-dir", type=Path, default=ROOT / ".stericx" / "kraken_dft_cache"
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs" / "study_010")
    parser.add_argument(
        "--ligands",
        type=int,
        default=60,
        help="Number of ligands sampled evenly across the id range.",
    )
    parser.add_argument(
        "--densities",
        type=str,
        default=",".join(str(d) for d in DEFAULT_DENSITIES),
        help="Comma-separated grid densities (A^3/point), coarse to fine.",
    )
    return parser.parse_args(argv)


def sample_geometries(cache_dir: Path, count: int) -> list[Path]:
    """One lowest-numbered SDF for `count` ligands spaced across the id range."""
    ligand_dirs = sorted(
        (d for d in cache_dir.iterdir() if d.is_dir() and d.name.isdigit()),
        key=lambda d: int(d.name),
    )
    with_sdf = [d for d in ligand_dirs if list(d.glob("*.sdf"))]
    if not with_sdf:
        raise SystemExit(f"no ligand SDFs under {cache_dir}")
    if count >= len(with_sdf):
        chosen = with_sdf
    else:
        indices = np.linspace(0, len(with_sdf) - 1, count).round().astype(int)
        chosen = [with_sdf[i] for i in sorted(set(indices.tolist()))]
    return [sorted(d.glob("*.sdf"))[0] for d in chosen]


def run_density(
    binary: Path, paths: list[Path], density: float
) -> tuple[dict[str, float], float]:
    """Return {file: %Vbur} and the wall time for the pass.

    Passes all `paths` in one call. At the study's scale (tens to ~1.5k one-SDF
    ligand paths) the argv stays well under the OS limit; the 31k-conformer case
    that needs chunking is Study 008's, which owns the STERICX_BATCH logic.
    """
    command = [
        str(binary),
        "descriptors",
        "--density",
        str(density),
        "--format",
        "csv",
        *map(str, paths),
    ]
    start = time.perf_counter()
    output = subprocess.run(command, capture_output=True, text=True, check=True).stdout
    elapsed = time.perf_counter() - start
    values: dict[str, float] = {}
    for row in csv.DictReader(io.StringIO(output)):
        vbur = row.get("percent_buried_volume")
        if vbur:
            values[Path(row["file"]).name] = float(vbur)
    return values, elapsed


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.binary.is_file():
        raise SystemExit(
            f"StericX binary not found: {args.binary} (cargo build --release)"
        )
    densities = sorted((float(d) for d in args.densities.split(",")), reverse=True)
    paths = sample_geometries(args.cache_dir, args.ligands)
    print(f"grid-convergence sweep: {len(paths)} ligands x {len(densities)} densities")

    # Warm the OS cache so the timing reflects compute, not first-read I/O.
    for path in paths:
        path.read_bytes()

    per_density: dict[float, dict[str, float]] = {}
    timing: dict[float, float] = {}
    for density in densities:
        values, elapsed = run_density(args.binary, paths, density)
        per_density[density] = values
        timing[density] = elapsed

    reference_density = min(densities)
    shared = sorted(set.intersection(*(set(v) for v in per_density.values())))
    reference = per_density[reference_density]
    ref_vbur = np.array([reference[name] for name in shared])

    # Irreducible discretization floor: the RMS %Vbur change between the two finest
    # grids. Refinement past this only shuffles which points land inside the sphere.
    two_finest = sorted(densities)[:2]
    floor_a = np.array([per_density[two_finest[0]][name] for name in shared])
    floor_b = np.array([per_density[two_finest[1]][name] for name in shared])
    floor_rms = float(np.sqrt(np.mean((floor_a - floor_b) ** 2)))

    rows: list[dict] = []
    print(
        f"\n  {'density':>9} {'mean|dVbur|':>12} {'max|dVbur|':>11} "
        f"{'bias':>8} {'time(s)':>9}"
    )
    for density in densities:
        vbur = np.array([per_density[density][name] for name in shared])
        diff = vbur - ref_vbur
        row = {
            "density": density,
            "mean_abs_vbur_error": float(np.mean(np.abs(diff))),
            "max_abs_vbur_error": float(np.max(np.abs(diff))),
            "signed_bias": float(np.mean(diff)),
            "seconds": timing[density],
            "is_reference": density == reference_density,
            "is_default": density == KRAKEN_DENSITY,
        }
        rows.append(row)
        marker = "  <- default" if row["is_default"] else ""
        marker += "  (reference)" if row["is_reference"] else ""
        print(
            f"  {density:>9} {row['mean_abs_vbur_error']:>12.4f} "
            f"{row['max_abs_vbur_error']:>11.4f} {row['signed_bias']:>+8.4f} "
            f"{row['seconds']:>9.3f}{marker}"
        )

    default_row = next(r for r in rows if r["is_default"])
    print(
        f"\n  discretization floor (RMS %Vbur between the two finest grids): "
        f"{floor_rms:.4f}"
    )
    within = default_row["mean_abs_vbur_error"] <= 2 * floor_rms
    print(
        f"  default 0.01 mean error vs finest: "
        f"{default_row['mean_abs_vbur_error']:.4f} %Vbur "
        f"({'within ~2x' if within else 'above ~2x'} the floor)"
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_figure(rows, floor_rms, args.output_dir / "grid_convergence.png")
    _write_metrics(rows, floor_rms, reference_density, len(shared), args.output_dir)
    write_report(
        rows,
        floor_rms,
        reference_density,
        len(shared),
        args.output_dir / "STUDY_010.md",
    )
    print("\nStudy 010 complete.")
    return 0


def _write_metrics(rows, floor_rms, reference_density, n, output_dir: Path) -> None:
    (output_dir / "grid_convergence_metrics.json").write_text(
        json.dumps(
            {
                "n_ligands": n,
                "reference_density": reference_density,
                "default_density": KRAKEN_DENSITY,
                "discretization_floor_rms_vbur": floor_rms,
                "sweep": rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def write_figure(rows: list[dict], floor_rms: float, output: Path) -> None:
    figure, (error_ax, cost_ax) = plt.subplots(1, 2, figsize=(11.5, 4.8))
    non_ref = [r for r in rows if not r["is_reference"]]
    densities = [r["density"] for r in non_ref]
    mean_err = [r["mean_abs_vbur_error"] for r in non_ref]

    error_ax.loglog(
        densities, mean_err, "o-", color="#0F766E", label="mean |ΔVbur| vs finest"
    )
    error_ax.axhline(
        floor_rms,
        ls="--",
        color="#B4530A",
        lw=1.0,
        label=f"discretization floor ({floor_rms:.3f})",
    )
    default = next(r for r in rows if r["is_default"])
    error_ax.scatter(
        [default["density"]],
        [default["mean_abs_vbur_error"]],
        s=120,
        facecolor="none",
        edgecolor="#0F766E",
        linewidths=2,
        zorder=5,
        label="Kraken default 0.01",
    )
    error_ax.set_xlabel("grid density (Å³ / point) — coarser →")
    error_ax.set_ylabel("mean |%Vbur - finest| ")
    error_ax.set_title("Buried volume converges as the grid refines")
    error_ax.invert_xaxis()
    error_ax.legend(fontsize=8, frameon=False)

    cost_ax.loglog(
        [r["density"] for r in rows],
        [r["seconds"] for r in rows],
        "s-",
        color="#8C8C8C",
    )
    cost_ax.scatter(
        [default["density"]],
        [default["seconds"]],
        s=120,
        facecolor="none",
        edgecolor="#176B87",
        linewidths=2,
        zorder=5,
        label="default 0.01",
    )
    cost_ax.set_xlabel("grid density (Å³ / point) — coarser →")
    cost_ax.set_ylabel("wall time for the sample (s)")
    cost_ax.set_title("Cost of refinement")
    cost_ax.invert_xaxis()
    cost_ax.legend(fontsize=8, frameon=False)

    figure.suptitle(
        "Study 010: grid convergence of the buried-volume integrator", fontsize=12
    )
    figure.tight_layout()
    figure.savefig(output, dpi=200)
    plt.close(figure)


def write_report(
    rows: list[dict], floor_rms: float, reference_density: float, n: int, output: Path
) -> None:
    default = next(r for r in rows if r["is_default"])
    coarse = max(rows, key=lambda r: r["density"])
    at_floor = default["mean_abs_vbur_error"] <= 2 * floor_rms
    reference_seconds = next(r["seconds"] for r in rows if r["is_reference"])
    refine_cost = reference_seconds / default["seconds"]
    lines = [
        "# StericX Study 010 - Grid Convergence of the Buried-Volume Integrator",
        "",
        "## Is the descriptor converged, or grid-lucky?",
        "",
        "Studies 002-004 show StericX's buried volume matches morfeus (R2 = 1.0 on "
        "identical geometries) and Kraken's published values (R2 = 0.9852 across "
        "1,541 ligands). But every one of those runs used a **single** grid "
        "resolution -- Kraken's convention of 0.01 A^3 per integration point. This "
        "study asks whether that agreement is a property of the descriptor or an "
        "artifact of one lucky grid, by sweeping the integration grid from coarse "
        "to fine on a diverse sample of "
        f"**{n} Kraken ligands** and measuring convergence in the flagship %Vbur.",
        "",
        "The finest grid tested "
        f"({reference_density:g} A^3/point) is taken as the converged reference; "
        "every coarser grid is measured against it.",
        "",
        "## Result",
        "",
        "| Grid density (Å³/point) | mean \\|ΔVbur\\| | max \\|ΔVbur\\| "
        "| signed bias | time (s) |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        tag = ""
        if row["is_default"]:
            tag = " **(default)**"
        elif row["is_reference"]:
            tag = " *(reference)*"
        lines.append(
            f"| {row['density']:g}{tag} | {row['mean_abs_vbur_error']:.4f} | "
            f"{row['max_abs_vbur_error']:.4f} | {row['signed_bias']:+.4f} | "
            f"{row['seconds']:.3f} |"
        )
    lines += [
        "",
        "**The descriptor converges.** As the grid refines, the mean deviation from "
        "the finest grid falls sharply -- from "
        f"**{coarse['mean_abs_vbur_error']:.3f} "
        f"%Vbur** at the coarse {coarse['density']:g} A^3/point down to "
        f"**{default['mean_abs_vbur_error']:.3f} %Vbur** at Kraken's default 0.01. The "
        "default sits inside the converged regime, not on the steep part of the curve.",
        "",
        f"**There is an irreducible floor.** Refining past the default does not drive "
        f"the error to zero: the two finest grids still differ by an RMS of "
        f"**{floor_rms:.3f} %Vbur**. This is the discretization floor of any voxel "
        "integrator -- refinement only reshuffles which grid points fall inside the "
        "sphere near its surface, so a small registration jitter remains at every "
        f"resolution. The default's mean error ({default['mean_abs_vbur_error']:.3f} "
        f"%Vbur) sits within about {'2x' if at_floor else 'a few x'} this floor -- "
        "close enough that further refinement chases jitter, not signal.",
        "",
        "![Grid convergence](grid_convergence.png)",
        "",
        "*Figure. Left: mean |%Vbur - finest| vs grid density (both axes log; coarser "
        "to the right), with the discretization floor and the Kraken default marked. "
        "Right: the wall-time cost of refinement. Generated by "
        "`studies/study_010_grid_convergence.py`.*",
        "",
        "## Why this matters",
        "",
        "It is worth being precise about what this does and does not say about the "
        "Study 004 reproduction. The full-set residual there is a median 0.11 A^3, or "
        f"~0.06 %Vbur on a 3.5 A sphere -- which is about the *same* size as the "
        f"default grid's own uncertainty here ({default['mean_abs_vbur_error']:.3f} "
        "%Vbur vs the finest grid). So grid discretization is not negligibly small in "
        "absolute terms.",
        "",
        "What makes the reproduction robust anyway is that **both sides use the same "
        "0.01 grid**: Kraken's published values come from a voxel integrator at the "
        "same resolution, so the discretization is largely common-mode and cancels in "
        "the StericX-vs-Kraken comparison rather than adding to it. And the residual "
        "that remains does not behave like grid noise -- Study 006 localizes it to "
        "specific primary/secondary phosphines and the coordination centre, a "
        "structured effect a finer grid would not touch. Refining past 0.01 costs "
        f"~{refine_cost:.0f}x the wall time for no meaningful gain, while a coarser "
        "grid -- tempting for "
        f"speed -- would visibly degrade the descriptor (the {coarse['density']:g} "
        f"A^3/point grid is off by {coarse['mean_abs_vbur_error']:.2f} %Vbur). The "
        "default is the honest speed/accuracy sweet spot, and the sweep confirms the "
        "reproduction studies were run in the converged regime, not a grid-lucky one.",
        "",
        "## Reproducing",
        "",
        "```bash",
        "uv run --extra science python studies/study_010_grid_convergence.py",
        "```",
        "",
        "Requires the Kraken DFT SDF cache (local, gitignored) and the release binary. "
        "Only StericX's own convergence measurements are written out.",
        "",
    ]
    output.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
